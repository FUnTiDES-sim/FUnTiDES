#!/usr/bin/env python3
"""
BP5 File Inspector - similar to bpls
Shows structure and metadata of ADIOS2 BP5 files
"""

import sys
import os
import struct
import json
from pathlib import Path


def inspect_bp5(bp_dir):
    """Inspect BP5 directory structure"""

    if not os.path.isdir(bp_dir):
        print(f"Error: {bp_dir} is not a directory")
        return False

    print(f"Inspecting BP5 file: {bp_dir}")
    print("=" * 70)

    # List files
    files = sorted(os.listdir(bp_dir))
    print("\nFiles in directory:")
    for f in files:
        filepath = os.path.join(bp_dir, f)
        size = os.path.getsize(filepath)
        size_mb = size / (1024 * 1024)
        print(f"  {f:20s} {size_mb:10.1f} MB  ({size:,} bytes)")

    total_size = sum(os.path.getsize(os.path.join(bp_dir, f)) for f in files)
    print(f"  {'Total':20s} {total_size/(1024*1024):10.1f} MB")

    # Try to read metadata
    md_file = os.path.join(bp_dir, "md.0")
    if os.path.exists(md_file):
        print("\n" + "=" * 70)
        print("Metadata (md.0):")
        print("=" * 70)
        try:
            with open(md_file, 'rb') as f:
                # BP5 metadata format is complex, try to extract some text
                data = f.read()
                # Look for printable strings
                strings = []
                current = []
                for byte in data:
                    if 32 <= byte < 127:  # Printable ASCII
                        current.append(chr(byte))
                    else:
                        if len(current) > 4:
                            strings.append(''.join(current))
                        current = []
                if current:
                    strings.append(''.join(current))

                # Print interesting strings
                for s in strings[:20]:  # First 20 strings
                    if len(s) > 3:
                        print(f"  {s}")
        except Exception as e:
            print(f"  Could not parse: {e}")

    # Try to read profiling data
    prof_file = os.path.join(bp_dir, "profiling.json")
    if os.path.exists(prof_file):
        print("\n" + "=" * 70)
        print("Profiling Information (profiling.json):")
        print("=" * 70)
        try:
            with open(prof_file, 'r') as f:
                prof_data = json.load(f)
                print(json.dumps(prof_data, indent=2)[:1000])  # First 1000 chars
        except Exception as e:
            print(f"  Could not parse: {e}")

    # Analyze data.0 file
    data_file = os.path.join(bp_dir, "data.0")
    if os.path.exists(data_file):
        print("\n" + "=" * 70)
        print("Data Analysis (data.0):")
        print("=" * 70)
        size = os.path.getsize(data_file)
        num_floats = size // 4
        num_doubles = size // 8

        print(f"  Total size: {size:,} bytes ({size/(1024*1024):.1f} MB)")
        print(f"  As float32: {num_floats:,} values")
        print(f"  As float64: {num_doubles:,} values")

        # If you know grid dimensions
        print("\n  Possible grid dimensions (as float32):")

        # Test common sizes
        test_sizes = [
            (201, 201, 201),
            (256, 256, 256),
            (100, 100, 100),
            (512, 512, 512),
            (128, 128, 128),
        ]

        for nx, ny, nz in test_sizes:
            grid_size = nx * ny * nz
            if num_floats % grid_size == 0:
                num_steps = num_floats // grid_size
                remainder = num_floats % grid_size
                print(f"    {nx}³: {num_steps} timesteps (remainder: {remainder})")
            elif (num_floats // grid_size) > 0:
                num_steps = num_floats // grid_size
                remainder = num_floats % grid_size
                print(f"    {nx}³: {num_steps} complete timesteps + {remainder} values")

    print("\n" + "=" * 70)
    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 inspect_bp5.py <bp_file_or_directory>")
        print("\nExample:")
        print("  python3 inspect_bp5.py snapshots.bp")
        return 1

    bp_path = sys.argv[1]

    if not os.path.exists(bp_path):
        print(f"Error: {bp_path} not found")
        return 1

    if os.path.isdir(bp_path):
        success = inspect_bp5(bp_path)
    else:
        print(f"Error: {bp_path} is not a directory (BP5 is stored as a directory)")
        return 1

    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
