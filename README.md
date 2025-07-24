# HILLMAPPER

**HILLMAPPER** is a high-performance **C++ extension** (via [PyBind11](https://pybind11.readthedocs.io)) that provides a fast, parallel **hill-climbing algorithm** implementation.  

## Features

- **Blazing fast** C++ core compiled with [PyBind11](https://pybind11.readthedocs.io).  
- Provides `find_optimal_bijection` — a **parallel hill-climbing optimizer**.  
- Precompiled for **Ubuntu (manylinux)** and **macOS (universal2)**.  
- **Windows support coming soon.**

## Requirements

- **Python 3.11+** 
- Runs on **Ubuntu** and **macOS** (x86_64 & Apple Silicon)

## Installation

Install from PyPI
```bash
pip install hillmapper
```

## Usage example
```python
#!/usr/bin/env python3
# import it first!
from hillmapper import find_optimal_bijection


def main():
    # Frequency counters for Redset scansets
    redset_scanset_counters = [4, 2, 3]

    # Redset scansets (1-based table IDs, since C++ uses 1..N internally)
    redset_scansets = [
        [1, 2],   # Scanset 0
        [2, 3],   # Scanset 1
        [1, 3],   # Scanset 2
    ]

    # Benchmark scansets (also 1-based)
    tpcds_scansets = [
        [1, 2],   # Benchmark scanset 0
        [2, 3],   # Benchmark scanset 1
    ]

    # Map each Redset table (1-based) to the indices of the scansets it affects
    r_table_to_scansets = [
        [],       # index 0 unused (placeholder so tables match 1..N)
        [0, 2],   # Table 1 in scansets 0, 2
        [0, 1],   # Table 2 in scansets 0, 1
        [1, 2],   # Table 3 in scansets 1, 2
    ]

    n_redset_tables = 3  # Tables 1..3
    n_tpcds_tables = 2   # Tables 1..2

    best_distance, best_assignment = find_optimal_bijection(
        2,                          # n_threads
        50,                         # n_iterations_per_thread
        redset_scanset_counters,
        redset_scansets,
        tpcds_scansets,
        r_table_to_scansets,
        n_redset_tables,
        n_tpcds_tables
    )

    print("=== Hillmapper Optimization Demo ===")
    print(f"Best Manhattan distance: {best_distance}")
    print("Optimal table assignment (redset_table → tpcds_table):")
    for redset_table, tpcds_table in best_assignment.items():
        print(f"  Redset table {redset_table} → TPCDS table {tpcds_table}")


if __name__ == "__main__":
    main()
```
