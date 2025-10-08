#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
from pathlib import Path


def run_once(threads, extra_pytest_args, bench_root, marker):
    env = os.environ.copy()
    for var in ("KOKKOS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT"):
        env[var] = str(threads)

    out_dir = bench_root
    out_dir.mkdir(parents=True, exist_ok=True)
    outfile = out_dir / f"python_{marker}_t{threads}_latest.json"

    cmd = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "tests/benchmarks/python",
        "--benchmark-only",
        f"--benchmark-json={outfile}",
    ]
    if marker:
        cmd.extend(("-m", marker))
    cmd += extra_pytest_args

    print(f"[RUN] threads={threads} -> {outfile}")
    subprocess.check_call(cmd, env=env)


def parse_args():
    p = argparse.ArgumentParser(
        description="Run python benchmarks for multiple thread counts."
    )
    p.add_argument(
        "--threads",
        type=str,
        default="2,4,8,16,32,64",
        help="Comma-separated list of thread counts (default: 2,4,8,16,32,64)",
    )
    p.add_argument(
        "--marker",
        type=str,
        default="",
        help="Optional pytest marker expression to filter benchmarks",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build/benchmarks/python"),
        help="Output directory for JSON results",
    )
    p.add_argument(
        "pytest_args",
        nargs=argparse.REMAINDER,
        help="Additional arguments passed through to pytest after '--'",
    )
    return p.parse_args()


def main():
    args = parse_args()
    thread_list = [int(t) for t in args.threads.split(",") if t.strip()]
    for t in thread_list:
        run_once(t, args.pytest_args, args.output_dir, args.marker)
    print("All Python benchmarks completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
