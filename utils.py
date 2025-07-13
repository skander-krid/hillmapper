import duckdb
import os
import json

redset = "'serverless_full.parquet'"
db = duckdb.connect("../../v0.1.1/redbench/tpcds/10gb.duckdb", read_only=True)

def count_tables(scansets):
    return max([max(scanset) for scanset in scansets if scanset])

def get_redset_scansets(
        instance_id=None,
        user_id=None,
        database_id=None,
        query_type="select",
        limit=1000,
        min_num_tables=2
):
    redset_scansets = duckdb.sql(f"""
        select read_table_ids
        from {redset}
        where true
            and read_table_ids is not null
            {f"and instance_id = {instance_id}" if instance_id is not None else ""}
            {f"and user_id = {user_id}" if user_id is not None else ""}
            {f"and database_id = {database_id}" if database_id is not None else ""}
            {f"and query_type = '{query_type}'" if query_type is not None else ""}
            {f"and len(regexp_extract_all(read_table_ids, ',')) >= " + str(min_num_tables-1) if min_num_tables else ""}
        limit {limit}
    """).fetchall()
    redset_scansets = map(lambda x: x[0], redset_scansets)
    redset_scansets = map(lambda x: x.split(','), redset_scansets)
    redset_scansets = map(lambda x: [int(i) for i in x], redset_scansets)
    redset_scansets = map(sorted, redset_scansets)
    redset_scansets = map(tuple, redset_scansets)
    redset_scansets = filter(lambda x: len(x) > 0, redset_scansets)
    redset_scansets = list(redset_scansets)
    redset_scansets, _ = translate_scansets(redset_scansets)
    validate_scansets(redset_scansets)
    return redset_scansets, count_tables(redset_scansets)

def get_tpcds_table_ids():
    table_to_id = dict()
    with open("tpcds_schema.sql", "r") as file:
        schema = [line for line in file.readlines() if line.startswith("CREATE TABLE")]
    for idx, line in enumerate(schema):
        table_name = line.split("CREATE TABLE ")[1].split("(")[0].strip()
        table_to_id[table_name] = idx + 1
    return table_to_id

def extract_scanset_from_profile(node, table_to_id=None):
    if table_to_id is None:
        table_to_id = get_tpcds_table_ids()
    scanset = []
    if node.get("operator_type") == "TABLE_SCAN":
        table = node.get("extra_info", {}).get("Table")
        if table:
            scanset.append(table_to_id[table.lower()])
    for child in node.get("children", []):
        scanset += extract_scanset_from_profile(child, table_to_id)
    return scanset

def extract_scanset_from_query(template):
    filepath = os.path.join("../redbench_v2/matching/tpcds/queries", str(template), "0.sql")
    with open(filepath, "r") as file:
        sql = " ".join(file.readlines())
    db.execute("PRAGMA enable_profiling='json'")
    db.execute("PRAGMA profiling_output = '/tmp/profile.json';")
    try:
        db.execute(sql)
    except Exception as e:
        return None
    with open("/tmp/profile.json", "r") as file:
        profile = file.read()
    return extract_scanset_from_profile(json.loads(profile))

def get_tpcds_scansets():
    tpcds_scansets = []
    for template in range(1, 100):
        scanset = extract_scanset_from_query(template)
        if scanset is not None:
            tpcds_scansets.append(scanset)
    tpcds_scansets = map(set, tpcds_scansets)
    tpcds_scansets = filter(lambda x: len(x) > 0, tpcds_scansets)
    tpcds_scansets = map(sorted, tpcds_scansets)
    tpcds_scansets = map(tuple, tpcds_scansets)
    tpcds_scansets = list(tpcds_scansets)
    validate_scansets(tpcds_scansets)
    return tpcds_scansets, count_tables(tpcds_scansets)

# Remove non-existing tables.
def translate_scansets(scansets):
    tables = sorted(list(set([table for scanset in scansets for table in scanset])))
    conversion = dict()
    for idx, table in enumerate(tables):
        conversion[table] = idx + 1
    return [tuple(sorted(list({conversion[table] for table in scanset}))) for scanset in scansets], conversion

def validate_scansets(scansets):
    # Make sure there are no duplicates in the scansets and that they are sorted.
    for scanset in scansets:
        if not isinstance(scanset, tuple):
            raise ValueError(f"Scanset {scanset} is not a tuple.")
        if len(scanset) == 0:
            raise ValueError("Empty scanset found.")
        if not all(isinstance(table, int) for table in scanset):
            raise ValueError(f"Scanset {scanset} contains non-integer table IDs.")
        if not all(table > 0 for table in scanset):
            raise ValueError(f"Scanset {scanset} contains non-positive table IDs.")
        if not all(scanset[i] < scanset[i + 1] for i in range(len(scanset) - 1)):
            raise ValueError(f"Scanset {scanset} is not sorted.")
        if len(scanset) != len(set(scanset)):
            raise ValueError(f"Scanset {scanset} contains duplicate table IDs.")

    # Make sure that all tables between min table and max table are present.
    min_table = min(table for scanset in scansets for table in scanset)
    max_table = max(table for scanset in scansets for table in scanset)
    for table in range(min_table, max_table + 1):
        if not any(table in scanset for scanset in scansets):
            raise ValueError(f"Table {table} is missing from all scansets.")