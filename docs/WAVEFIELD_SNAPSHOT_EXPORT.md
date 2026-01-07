# Wavefield Snapshot Export Documentation

This document describes how to save and export wavefield snapshots from the FD solver Python bindings.

## Overview

The FD solver Python bindings now support saving wavefield snapshots during simulation. This allows you to:
- Visualize wave propagation over time
- Analyze wavefield characteristics
- Debug simulation behavior
- Create animations of wave propagation

## Implementation

### C++ Function: `get_wavefield_snapshot()`

Location: `src/solver/fd/pywrap/src/bindings.cpp:391-429`

```cpp
py::array_t<float> get_wavefield_snapshot(
    fdtd::kernel::FdtdKernels& kernels,
    int nx, int ny, int nz,
    int lx, int ly, int lz,
    int buffer_idx
);
```

**Parameters:**
- `kernels`: FdtdKernels object containing the wavefield data
- `nx, ny, nz`: Interior grid dimensions (without ghost cells)
- `lx, ly, lz`: Stencil half-widths
- `buffer_idx`: Current buffer index (0 or 1)

**Returns:**
- NumPy array of shape `(nx, ny, nz)` containing the pressure field

**Implementation Details:**
- Extracts data from the Kokkos `pnGlobal` array
- Removes ghost cells, returning only the interior domain
- Copies data to a NumPy array for easy manipulation in Python
- Uses row-major (C-style) ordering for NumPy compatibility

### Python Helper: `save_snapshot()`

Location: `tests/benchmarks/python/test_fd_solver_example.py:23-44`

```python
def save_snapshot(kernels, nx, ny, nz, lx, ly, lz, itime, buffer_idx,
                  output_dir="snapshots"):
    """Save wavefield snapshot to disk."""
```

**Parameters:**
- `kernels`: FdtdKernels object
- `nx, ny, nz`: Grid dimensions
- `lx, ly, lz`: Stencil half-widths
- `itime`: Current time step number
- `buffer_idx`: Current buffer index (0 or 1)
- `output_dir`: Directory to save snapshots (default: "snapshots")

**Returns:**
- String: Path to the saved snapshot file

**File Format:**
- Files are saved as NumPy binary format (`.npy`)
- Naming convention: `wavefield_NNNNNN.npy` where NNNNNN is the zero-padded time step
- Each file contains a 3D float32 array

## Usage Example

### Basic Usage

```python
from pyfuntides import fd_solver
import numpy as np
import os

# Initialize solver (see test_fd_solver_example.py for full setup)
# ... initialization code ...

# Time-stepping loop with snapshot saving
i1, i2 = 0, 1
snapshot_interval = 10  # Save every 10 time steps

for itime in range(n_steps):
    # Compute one time step
    solver.compute_one_stepPML(itime, i1, i2)
    i1, i2 = i2, i1  # Swap buffers

    # Save snapshot at specified intervals
    if itime % snapshot_interval == 0:
        wavefield = fd_solver.get_wavefield_snapshot(
            kernels, nx, ny, nz,
            stencils.lx, stencils.ly, stencils.lz,
            i2  # Current buffer index
        )

        # Save to file
        os.makedirs("snapshots", exist_ok=True)
        filename = f"snapshots/wavefield_{itime:06d}.npy"
        np.save(filename, wavefield)
        print(f"Saved: {filename}")
```

### Loading and Visualizing Snapshots

```python
import numpy as np
import matplotlib.pyplot as plt

# Load a snapshot
wavefield = np.load('snapshots/wavefield_000100.npy')

# Print basic information
print(f"Shape: {wavefield.shape}")
print(f"Data type: {wavefield.dtype}")
print(f"Min: {wavefield.min():.6e}, Max: {wavefield.max():.6e}")
print(f"Mean: {wavefield.mean():.6e}")

# Visualize a 2D slice (XY plane at z=100)
plt.figure(figsize=(10, 8))
plt.imshow(wavefield[:, :, 100].T, cmap='seismic', aspect='auto')
plt.colorbar(label='Pressure')
plt.title('Wavefield at z=100')
plt.xlabel('X index')
plt.ylabel('Y index')
plt.savefig('wavefield_slice.png', dpi=150)
plt.close()

# Visualize a 2D slice (XZ plane at y=100)
plt.figure(figsize=(10, 8))
plt.imshow(wavefield[:, 100, :].T, cmap='seismic', aspect='auto')
plt.colorbar(label='Pressure')
plt.title('Wavefield at y=100')
plt.xlabel('X index')
plt.ylabel('Z index')
plt.savefig('wavefield_xz_slice.png', dpi=150)
plt.close()
```

### Creating an Animation

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from glob import glob

# Load all snapshots
snapshot_files = sorted(glob('snapshots/wavefield_*.npy'))
print(f"Found {len(snapshot_files)} snapshots")

# Create animation
fig, ax = plt.subplots(figsize=(10, 8))

def update(frame):
    wavefield = np.load(snapshot_files[frame])
    # Show XY slice at z=100
    ax.clear()
    im = ax.imshow(wavefield[:, :, 100].T, cmap='seismic',
                   aspect='auto', vmin=-2e4, vmax=2e4)
    ax.set_title(f'Wavefield at time step {frame * 10}')
    ax.set_xlabel('X index')
    ax.set_ylabel('Y index')
    return [im]

anim = animation.FuncAnimation(fig, update, frames=len(snapshot_files),
                              interval=100, blit=True)
anim.save('wavefield_animation.mp4', writer='ffmpeg', fps=10, dpi=150)
plt.close()
```

## Configuration

### Enabling Snapshot Saving

Snapshots are controlled via the `FdtdOptions.output` configuration:

```python
options = fd_solver.FdtdOptions()

# Enable snapshot saving
options.output.save_snapshots = True
options.output.snapshot_interval = 10  # Save every 10 time steps
```

### Performance Considerations

**File Size:**
- Each snapshot file is approximately 32 MB for a 200×200×200 grid (float32)
- Formula: `4 bytes × nx × ny × nz`
- Example: For 200³ grid: `4 × 200³ = 32,000,000 bytes ≈ 30.5 MB`

**Storage Requirements:**
- 900 time steps with interval=10 → 90 snapshots → ~2.7 GB total
- Consider using compression or lower sampling rates for long simulations

**Performance Impact:**
- I/O operations can slow down the simulation
- Larger intervals reduce impact but provide less temporal resolution
- Consider using asynchronous I/O for production runs

### Recommended Intervals

| Simulation Length | Recommended Interval | Total Snapshots | Storage (200³) |
|-------------------|----------------------|-----------------|----------------|
| 1,000 steps       | 10                   | 100             | ~3 GB          |
| 10,000 steps      | 50                   | 200             | ~6 GB          |
| 100,000 steps     | 100-500              | 200-1000        | 6-30 GB        |

## File Format Details

### NumPy .npy Format

The snapshots are saved in NumPy's native binary format (`.npy`):
- **Advantages:**
  - Fast read/write
  - Preserves array metadata (shape, dtype)
  - Easy to load in Python
  - Compact storage

- **Disadvantages:**
  - Python/NumPy specific
  - Not self-documenting
  - No compression by default

### Alternative Formats

You can modify `save_snapshot()` to use other formats:

**HDF5 (for larger datasets):**
```python
import h5py

def save_snapshot_hdf5(wavefield, filename, itime):
    with h5py.File(filename, 'a') as f:
        f.create_dataset(f'wavefield_{itime:06d}', data=wavefield,
                        compression='gzip', compression_opts=9)
```

**VTK (for ParaView visualization):**
```python
import vtk
from vtk.util import numpy_support

def save_snapshot_vtk(wavefield, filename):
    # Convert NumPy array to VTK
    vtk_data = numpy_support.numpy_to_vtk(wavefield.ravel(), deep=True)

    # Create image data
    image = vtk.vtkImageData()
    image.SetDimensions(wavefield.shape)
    image.GetPointData().SetScalars(vtk_data)

    # Write to file
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(image)
    writer.Write()
```

## Troubleshooting

### Module Import Error

```
AttributeError: module 'pyfuntides.fd_solver' has no attribute 'get_wavefield_snapshot'
```

**Solution:** Rebuild and reinstall the Python module:
```bash
cmake --build buildPy --target fd_solver_py -j8
cp buildPy/src/solver/fd/pywrap/fd_solver.cpython-310-x86_64-linux-gnu.so \
   buildPy/install/python/pyfuntides/fd_solver.cpython-310-x86_64-linux-gnu.so
```

### Out of Memory

If you run out of memory when saving snapshots:
- Increase `snapshot_interval` to save fewer snapshots
- Use compression (HDF5 with gzip)
- Save to a different disk with more space
- Reduce grid resolution

### Incorrect Buffer Index

If wavefields look incorrect, ensure you're using the correct buffer index:
- After swapping buffers with `i1, i2 = i2, i1`, the current field is in `i2`
- Pass `i2` as the `buffer_idx` parameter

## Related Files

- C++ bindings: `src/solver/fd/pywrap/src/bindings.cpp`
- Python example: `tests/benchmarks/python/test_fd_solver_example.py`
- Environment setup: `env.sh`

## See Also

- [FD Python Wrapper Documentation](FD_PYTHON_WRAPPER_COMPLETE.md)
- [Python Kokkos Lifecycle](PYTHON_KOKKOS_LIFECYCLE.md)
- [PML Binding Status](PML_BINDING_STATUS.md)
