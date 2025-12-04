#!/usr/bin/env python3
"""Quick verification of generated test data."""
import numpy as np

def read_sep_binary(header_file):
    """Read SEP header and binary data."""
    params = {}
    with open(header_file) as f:
        for line in f:
            if '=' in line and not line.startswith('#'):
                key, val = line.strip().split('=', 1)
                params[key.strip()] = val.strip()

    n1, n2, n3 = int(params['n1']), int(params['n2']), int(params['n3'])
    esize = int(params.get('esize', 4))
    data_format = params.get('data_format', 'native_float')

    import os
    bin_file = params['in']
    if not os.path.isabs(bin_file):
        bin_file = os.path.join(os.path.dirname(header_file), bin_file)

    dtype = np.float32 if esize == 4 else np.float64
    if 'xdr' in data_format:
        dtype = '>f4' if esize == 4 else '>f8'

    data = np.fromfile(bin_file, dtype=dtype)
    return data.reshape((n1, n2, n3), order='F')

# Test index data
print("Verifying test_index.H...")
data = read_sep_binary("test_index.H")
assert data[0,0,0] == 0, f"Expected 0, got {data[0,0,0]}"
assert data[1,0,0] == 1, f"Expected 1, got {data[1,0,0]}"
assert data[0,1,0] == 1000, f"Expected 1000, got {data[0,1,0]}"
assert data[0,0,1] == 1000000, f"Expected 1000000, got {data[0,0,1]}"
assert data[5,3,2] == 5 + 3*1000 + 2*1000000, f"Mismatch at (5,3,2)"
print("  PASSED")

print("\nAll verifications passed!")
