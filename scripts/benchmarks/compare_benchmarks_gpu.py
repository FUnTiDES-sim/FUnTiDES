#!/usr/bin/env python3

"""
Benchmark Comparison Tool

This script analyzes and compares the performance results of two different runs
of the finite element (FE) Cartesian solver benchmark suite.

It is designed to be used in a CI/CD pipeline or locally to detect performance
regressions or improvements between two codebase states (e.g., comparing a new
'feature' branch against the 'main' baseline).

Key Features:
- Parses Google Benchmark JSON outputs from two specified directories.
- Compares `cpu_time` to isolate algorithmic changes from background OS noise.
- Applies a ±2% threshold filter to ignore standard measurement fluctuations.
- Outputs a color-coded terminal report: Green for improvements, Red for regressions.

Usage:
    ./compare_benchmarks.py <path_to_reference_dir> <path_to_new_dir>
"""

import json
import argparse
from pathlib import Path
import sys

# The 4 files we expect to find in the benchmark directories
EXPECTED_FILES = [
    "bench_solver_fe_cartesian_struct_acoustic.json",
    "bench_solver_fe_cartesian_unstruct_acoustic.json",
    "bench_solver_fe_cartesian_struct_elastic.json",
    "bench_solver_fe_cartesian_unstruct_elastic.json",
]

# ANSI color codes for terminal output
COLOR_RED = "\033[91m"
COLOR_GREEN = "\033[92m"
COLOR_RESET = "\033[0m"


def load_benchmark_data(filepath: Path):
    """Loads a Google Benchmark JSON file and returns a dictionary {name: data}."""
    if not filepath.exists():
        return None

    with open(filepath, "r") as f:
        try:
            data = json.load(f)
            # Index each run by its name for easy comparison O(1) lookups
            return {b["name"]: b for b in data.get("benchmarks", [])}
        except json.JSONDecodeError:
            print(f"Error: The file {filepath} is not a valid JSON.", file=sys.stderr)
            return None


def format_diff(old_val, new_val):
    """Calculates and formats the percentage difference with color coding."""
    if old_val == 0:
        return "N/A"

    diff_pct = ((new_val - old_val) / old_val) * 100

    # A shorter time (negative) is a performance improvement (green)
    # A longer time (positive) is a performance degradation (red)
    if diff_pct < -2.0:  # 2% threshold to filter out standard measurement noise
        return f"{COLOR_GREEN}{diff_pct:+.2f}%{COLOR_RESET}"
    elif diff_pct > 2.0:
        return f"{COLOR_RED}{diff_pct:+.2f}%{COLOR_RESET}"
    else:
        return f"{diff_pct:+.2f}%"


def compare_directories(dir_ref: Path, dir_new: Path):
    print(f"Performance Comparison:")
    print(f"Reference (Base) : {dir_ref}")
    print(f"New (Target)     : {dir_new}")
    print("=" * 80)

    for filename in EXPECTED_FILES:
        ref_file = dir_ref / filename
        new_file = dir_new / filename

        print(f"\n📊 File: {filename}")

        ref_data = load_benchmark_data(ref_file)
        new_data = load_benchmark_data(new_file)

        if not ref_data or not new_data:
            print("  -> File missing in one of the directories, skipped.")
            continue

        # Print table header
        print(
            f"{'Benchmark Name':<65} | {'Ref (ms)':<10} | {'New (ms)':<10} | {'Diff (%)'}"
        )
        print("-" * 105)

        # Iterate through the keys of the reference file
        for name, ref_bench in ref_data.items():
            if name not in new_data:
                continue  # The benchmark was removed or renamed in the new version

            new_bench = new_data[name]

            # Google benchmark provides real_time and cpu_time.
            # We use cpu_time by default as it is generally more stable.
            ref_time = ref_bench.get("cpu_time", 0)
            new_time = new_bench.get("cpu_time", 0)

            diff_str = format_diff(ref_time, new_time)

            # Truncate the name if it's too long to maintain table alignment
            display_name = name if len(name) < 63 else name[:60] + "..."

            print(
                f"{display_name:<65} | {ref_time:<10.2f} | {new_time:<10.2f} | {diff_str}"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compares the JSON results of two Google Benchmark directories."
    )
    parser.add_argument(
        "dir_ref",
        type=str,
        help="Directory containing the reference results (e.g., main branch)",
    )
    parser.add_argument(
        "dir_new",
        type=str,
        help="Directory containing the new results (e.g., feature branch)",
    )

    args = parser.parse_args()

    dir_ref = Path(args.dir_ref)
    dir_new = Path(args.dir_new)

    if not dir_ref.is_dir() or not dir_new.is_dir():
        print("Error: Both arguments must point to valid directories.", file=sys.stderr)
        sys.exit(1)

    compare_directories(dir_ref, dir_new)
