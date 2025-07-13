import duckdb
import os
import json

redset = "'serverless_full.parquet'"
db = duckdb.connect("../../v0.1.1/redbench/tpcds/10gb.duckdb", read_only=True)

def get_redset_scansets():
    redset_scansets = duckdb.sql(f"""
        select
            read_table_ids
        from 'serverless_full.parquet'
        where true
            and database_id = 0
            and instance_id = 126
            and read_table_ids is not null
            and query_type = 'select'
    """).fetchall()
    redset_scansets = map(lambda x: x[0], redset_scansets)
    redset_scansets = map(lambda x: x.split(','), redset_scansets)
    redset_scansets = map(lambda x: [int(i) for i in x], redset_scansets)
    redset_scansets = map(sorted, redset_scansets)
    redset_scansets = map(tuple, redset_scansets)
    return list(redset_scansets)

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
        with open("profile.json", "w") as outfile:
            json.dump(json.loads(profile), outfile, indent=2)
    return extract_scanset_from_profile(json.loads(profile))

def get_tpcds_scansets():
    tpcds_scansets = []
    for template in range(1, 100):
        scanset = extract_scanset_from_query(template)
        if scanset is not None:
            tpcds_scansets.append(scanset)
    return list(map(tuple, tpcds_scansets))
