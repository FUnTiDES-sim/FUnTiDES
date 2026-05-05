#!/usr/bin/env python3
"""
BP5 File Inspector - Parse metadata to show variables and timesteps
"""

import sys
import os
import struct
import json


def read_uint64_le(data, offset):
    """Read little-endian uint64"""
    return struct.unpack('<Q', data[offset:offset+8])[0]


def read_uint32_le(data, offset):
    """Read little-endian uint32"""
    return struct.unpack('<I', data[offset:offset+4])[0]


def inspect_bp5_detailed(bp_dir):
    """Detailed inspection of BP5 metadata"""

    print(f"Inspecting BP5: {bp_dir}")
    print("=" * 80)

    # Read metadata file
    md_file = os.path.join(bp_dir, "md.0")
    if not os.path.exists(md_file):
        print("Error: md.0 not found")
        return False

    with open(md_file, 'rb') as f:
        md_data = f.read()

    print(f"\nMetadata file size: {len(md_data):,} bytes")

    # Try to parse BP5 metadata structure
    # BP5 format: Complex binary format, but we can extract some info

    # Look for variable names (typically stored as strings)
    print("\nSearching for variable names and metadata...")

    # Extract all printable strings
    strings = []
    current = []
    for i, byte in enumerate(md_data):
        if 32 <= byte < 127:  # Printable ASCII
            current.append((i, chr(byte)))
        else:
            if len(current) > 2:
                start_pos = current[0][0]
                text = ''.join(c[1] for c in current)
                strings.append((start_pos, text))
            current = []

    if current:
        start_pos = current[0][0]
        text = ''.join(c[1] for c in current)
        strings.append((start_pos, text))

    # Filter for interesting strings
    print("\nStrings found in metadata:")
    for pos, s in strings:
        # Look for variable names and important keywords
        if any(keyword in s.lower() for keyword in
               ['pressure', 'field', 'time', 'step', 'array', 'variable',
                'acoustic', 'receiver', 'coords', 'iteration']):
            print(f"  [{pos:6d}] {s}")

    # Try to read data file info
    data_file = os.path.join(bp_dir, "data.0")
    if os.path.exists(data_file):
        data_size = os.path.getsize(data_file)
        print(f"\nData file (data.0): {data_size:,} bytes")

        # Read first few bytes to see if there's a pattern
        with open(data_file, 'rb') as f:
            first_bytes = f.read(min(256, data_size))
            print(f"First bytes (hex):")
            for i in range(0, min(len(first_bytes), 64), 16):
                hex_str = ' '.join(f'{b:02x}' for b in first_bytes[i:i+16])
                ascii_str = ''.join(chr(b) if 32 <= b < 127 else '.' for b in first_bytes[i:i+16])
                print(f"  {i:04x}: {hex_str:48s} {ascii_str}")

    # Check md.idx file
    md_idx_file = os.path.join(bp_dir, "md.idx")
    if os.path.exists(md_idx_file):
        print(f"\nMetadata index file (md.idx): {os.path.getsize(md_idx_file):,} bytes")

    # Check mmd.0 file
    mmd_file = os.path.join(bp_dir, "mmd.0")
    if os.path.exists(mmd_file):
        with open(mmd_file, 'rb') as f:
            mmd_data = f.read()
        print(f"\nMeta-metadata file (mmd.0): {len(mmd_data):,} bytes")
        # Try to find strings
        strings = []
        current = []
        for byte in mmd_data:
            if 32 <= byte < 127:
                current.append(chr(byte))
            else:
                if len(current) > 2:
                    strings.append(''.join(current))
                current = []
        if current:
            strings.append(''.join(current))

        if strings:
            print("  Strings in mmd.0:")
            for s in strings[:10]:
                print(f"    {s}")

    print("\n" + "=" * 80)
    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 inspect_bp5_detailed.py <bp_directory>")
        return 1

    bp_path = sys.argv[1]

    if not os.path.isdir(bp_path):
        print(f"Error: {bp_path} is not a directory")
        return 1

    inspect_bp5_detailed(bp_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
