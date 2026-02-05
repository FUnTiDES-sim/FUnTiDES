#!/usr/bin/env python3
"""
ADIOS2 3D Pressure Field Visualization Script
==============================================

PURPOSE:
--------
Visualizes 3D pressure field data from ADIOS2 BP files across multiple timesteps.
Creates interactive 3D visualizations showing the pressure distribution in a cubic
domain for SEM geophysical acoustic simulations.

USAGE:
------
    python visualize_pressure.py <nx> <ny> <nz> [--file snapshots.bp] [--var PressureField]

ARGUMENTS:
----------
    nx, ny, nz : Grid dimensions in x, y, z directions
    --file     : Path to ADIOS2 BP file (default: snapshots.bp)
    --var      : Variable name to visualize (default: PressureField)
    --output   : Output directory for plots (default: plots/)
    --cmap     : Colormap (default: seismic)
    --slice    : Plot 2D slices instead of 3D volume (faster)
    --animate  : Create animation of all timesteps

EXAMPLES:
---------
    # Basic usage - visualize 201x201x201 grid
    python visualize_pressure.py 201 201 201

    # Specify custom file
    python visualize_pressure.py 201 201 201 --file my_data.bp

    # Create 2D slices (faster for large data)
    python visualize_pressure.py 201 201 201 --slice

    # Create animation
    python visualize_pressure.py 201 201 201 --animate

OUTPUT:
-------
Creates PNG files for each timestep showing:
- 3D volume rendering with isosurfaces
- Or 2D slice views (XY, XZ, YZ planes)
- Colorbar showing pressure range
- Timestep information

REQUIREMENTS:
-------------
    pip install adios2 numpy matplotlib

Optional for better 3D plots:
    pip install plotly pyvista

AUTHOR: SEM Geophysical Acoustics Simulation Team
VERSION: 2.0 (Fixed)
"""

import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def read_adios2_data(filename, variable_name, nx, ny, nz):
    """
    Read 3D pressure field data from ADIOS2 BP file using h5py or direct binary reading.

    Parameters:
    -----------
    filename : str
        Path to ADIOS2 BP file (can be BP5 directory or h5-compatible)
    variable_name : str
        Name of the variable to read
    nx, ny, nz : int
        Grid dimensions

    Returns:
    --------
    data : list of numpy arrays
        List of 3D arrays, one per timestep
    timesteps : list of int
        Timestep values
    """
    print(f"Reading {filename}...")

    if not os.path.exists(filename):
        print(f"Error: File not found: {filename}")
        return None, None

    data_list = []
    timestep_list = []

    # Try h5py if available (for HDF5 BP files)
    try:
        import h5py
        print("  Attempting to read with h5py...")
        with h5py.File(filename, 'r') as f:
            # List what's in the file
            print(f"  Available datasets: {list(f.keys())}")

            if variable_name not in f:
                print(f"  Error: Variable '{variable_name}' not found")
                print(f"  Available: {list(f.keys())}")
                return None, None

            dataset = f[variable_name]
            print(f"  Shape: {dataset.shape}, dtype: {dataset.dtype}")

            # Read all data or by steps
            if len(dataset.shape) == 4:
                # Shape: (steps, nx, ny, nz)
                for step in range(dataset.shape[0]):
                    data = dataset[step, :, :, :]
                    data_list.append(data.astype(np.float32))
                    timestep_list.append(step)
            elif len(dataset.shape) == 3:
                # Single timestep: (nx, ny, nz)
                data_list.append(dataset[:, :, :].astype(np.float32))
                timestep_list.append(0)
            else:
                print(f"  Error: Unexpected shape {dataset.shape}")
                return None, None

            print(f"Successfully read {len(data_list)} timesteps")
            return data_list, timestep_list

    except ImportError:
        print("  h5py not available, trying BP5 directory...")
    except Exception as e:
        print(f"  Error reading with h5py: {e}")
        print("  Trying BP5 directory...")

    # Fallback: Try reading BP file structure directly
    return read_bp5_directory(filename, variable_name, nx, ny, nz)


def read_bp5_directory(directory, variable_name, nx, ny, nz):
    """
    Read BP5 format data from directory structure.
    BP5 files are stored as directories with data.0, md.0, md.idx, mmd.0 files.

    Parameters:
    -----------
    directory : str
        Path to BP5 directory
    variable_name : str
        Name of the variable to read
    nx, ny, nz : int
        Grid dimensions

    Returns:
    --------
    data : list of numpy arrays
        List of 3D arrays, one per timestep
    timesteps : list of int
        Timestep values
    """
    import struct

    data_list = []
    timestep_list = []

    if not os.path.isdir(directory):
        print(f"  Error: {directory} is not a directory")
        return None, None

    try:
        # Check for BP5 structure files
        data_file = os.path.join(directory, "data.0")
        md_file = os.path.join(directory, "md.0")

        if not os.path.exists(data_file):
            print(f"  Error: {data_file} not found")
            return None, None

        print(f"  Found BP5 files: data.0, md.0, etc.")

        # For now, assume data is stored as consecutive float32 values
        # Read the entire data file
        with open(data_file, 'rb') as f:
            raw_data = f.read()

        # Convert to float32 array
        num_floats = len(raw_data) // 4
        data_array = np.frombuffer(raw_data, dtype=np.float32, count=num_floats)

        print(f"  Read {num_floats} float32 values from data.0")
        print(f"  Expected per timestep: {nx * ny * nz}")

        # Try to determine number of timesteps
        values_per_step = nx * ny * nz
        if num_floats % values_per_step == 0:
            num_steps = num_floats // values_per_step
            print(f"  Detected {num_steps} timesteps")

            for step in range(num_steps):
                start_idx = step * values_per_step
                end_idx = start_idx + values_per_step
                step_data = data_array[start_idx:end_idx].reshape((nx, ny, nz))
                data_list.append(step_data)
                timestep_list.append(step)
        else:
            # Can't evenly divide - just take what we can
            print(f"  Warning: Data size {num_floats} doesn't divide evenly by {values_per_step}")
            print(f"  Attempting to read as single timestep...")

            if num_floats >= values_per_step:
                step_data = data_array[:values_per_step].reshape((nx, ny, nz))
                data_list.append(step_data)
                timestep_list.append(0)

        if len(data_list) > 0:
            print(f"Successfully read {len(data_list)} timestep(s) from BP5 directory")
            return data_list, timestep_list
        else:
            print(f"  No data could be extracted")
            return None, None

    except Exception as e:
        print(f"  Error reading BP5 directory: {e}")
        import traceback
        traceback.print_exc()

    return None, None


def plot_3d_volume(data, timestep, output_file, cmap='seismic', vmin=None, vmax=None):
    """
    Create 3D volume visualization with scatter points.

    Parameters:
    -----------
    data : numpy array (3D)
        Pressure field data
    timestep : int
        Timestep number
    output_file : str
        Output filename
    cmap : str
        Colormap name
    vmin, vmax : float
        Min/max values for colormap
    """
    nx, ny, nz = data.shape

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Create coordinate grids
    x, y, z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing='ij')

    # Plot using scatter (subsample for performance)
    stride = max(1, min(nx, ny, nz) // 50)  # Subsample to ~50 points per dimension

    x_sub = x[::stride, ::stride, ::stride]
    y_sub = y[::stride, ::stride, ::stride]
    z_sub = z[::stride, ::stride, ::stride]
    data_sub = data[::stride, ::stride, ::stride]

    # Flatten for scatter plot
    scatter = ax.scatter(x_sub.ravel(),
                        y_sub.ravel(),
                        z_sub.ravel(),
                        c=data_sub.ravel(),
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                        alpha=0.6,
                        s=20)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Pressure Field - Timestep {timestep}')

    plt.colorbar(scatter, ax=ax, label='Pressure', shrink=0.8)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Saved: {output_file}")


def plot_2d_slices(data, timestep, output_file, cmap='seismic', vmin=None, vmax=None):
    """
    Create 2D slice views through the center of the domain.

    Parameters:
    -----------
    data : numpy array (3D)
        Pressure field data
    timestep : int
        Timestep number
    output_file : str
        Output filename
    cmap : str
        Colormap name
    vmin, vmax : float
        Min/max values for colormap
    """
    nx, ny, nz = data.shape

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    # Get middle slices
    mid_x = nx // 2
    mid_y = ny // 2
    mid_z = nz // 2

    slice_yz = data[mid_x, :, :]  # YZ plane at middle X
    slice_xz = data[:, mid_y, :]  # XZ plane at middle Y
    slice_xy = data[:, :, mid_z]  # XY plane at middle Z

    # Create figure with 3 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # XY slice
    im1 = axes[0, 0].imshow(slice_xy.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 0].set_title(f'XY Slice (Z={mid_z})')
    axes[0, 0].set_xlabel('X')
    axes[0, 0].set_ylabel('Y')
    plt.colorbar(im1, ax=axes[0, 0])

    # XZ slice
    im2 = axes[0, 1].imshow(slice_xz.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 1].set_title(f'XZ Slice (Y={mid_y})')
    axes[0, 1].set_xlabel('X')
    axes[0, 1].set_ylabel('Z')
    plt.colorbar(im2, ax=axes[0, 1])

    # YZ slice
    im3 = axes[1, 0].imshow(slice_yz.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 0].set_title(f'YZ Slice (X={mid_x})')
    axes[1, 0].set_xlabel('Y')
    axes[1, 0].set_ylabel('Z')
    plt.colorbar(im3, ax=axes[1, 0])

    # Statistics
    axes[1, 1].axis('off')
    stats_text = f"""
    Timestep: {timestep}
    Grid: {nx} x {ny} x {nz}

    Statistics:
    Min:  {data.min():.6f}
    Max:  {data.max():.6f}
    Mean: {data.mean():.6f}
    Std:  {data.std():.6f}
    """
    axes[1, 1].text(0.1, 0.5, stats_text, fontsize=12, family='monospace',
                    verticalalignment='center')

    fig.suptitle(f'Pressure Field 2D Slices - Timestep {timestep}', fontsize=14, y=0.995)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Saved: {output_file}")


def create_animation(data_list, timestep_list, output_file, mode='slice', cmap='seismic',
                    vmin=None, vmax=None):
    """
    Create animated GIF of all timesteps.

    Parameters:
    -----------
    data_list : list of numpy arrays
        Pressure field data for all timesteps
    timestep_list : list of int
        Timestep values
    output_file : str
        Output GIF filename
    mode : str
        'slice' or 'volume'
    cmap : str
        Colormap name
    vmin, vmax : float
        Min/max values for colormap (required for animation)
    """
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
    except ImportError:
        print("Error: Animation requires pillow. Install with: pip install pillow")
        return

    print(f"Creating animation with {len(data_list)} frames...")
    print(f"Animation range: [{vmin:.4f}, {vmax:.4f}]")

    if mode == 'slice':
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        def update(frame):
            data = data_list[frame]
            timestep = timestep_list[frame]
            nx, ny, nz = data.shape

            mid_x = nx // 2
            mid_y = ny // 2
            mid_z = nz // 2

            slice_xy = data[:, :, mid_z]
            slice_xz = data[:, mid_y, :]
            slice_yz = data[mid_x, :, :]

            for ax in axes:
                ax.clear()

            im1 = axes[0].imshow(slice_xy.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
            axes[0].set_title(f'XY Slice (Z={mid_z})')
            axes[0].set_xlabel('X')
            axes[0].set_ylabel('Y')

            im2 = axes[1].imshow(slice_xz.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
            axes[1].set_title(f'XZ Slice (Y={mid_y})')
            axes[1].set_xlabel('X')
            axes[1].set_ylabel('Z')

            im3 = axes[2].imshow(slice_yz.T, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
            axes[2].set_title(f'YZ Slice (X={mid_x})')
            axes[2].set_xlabel('Y')
            axes[2].set_ylabel('Z')

            fig.suptitle(f'Timestep {timestep}', fontsize=14)

            return [im1, im2, im3]

        anim = FuncAnimation(fig, update, frames=len(data_list), interval=200, blit=False)

        writer = PillowWriter(fps=5)
        anim.save(output_file, writer=writer)

        plt.close()
        print(f"Animation saved: {output_file}")

    else:
        print("3D volume animation not implemented. Use --slice mode.")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize 3D pressure field from ADIOS2 BP files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('nx', type=int, help='Number of nodes in X direction')
    parser.add_argument('ny', type=int, help='Number of nodes in Y direction')
    parser.add_argument('nz', type=int, help='Number of nodes in Z direction')
    parser.add_argument('--file', type=str, default='snapshots.bp',
                        help='ADIOS2 BP file path (default: snapshots.bp)')
    parser.add_argument('--var', type=str, default='PressureField',
                        help='Variable name to visualize (default: PressureField)')
    parser.add_argument('--output', type=str, default='plots',
                        help='Output directory (default: plots/)')
    parser.add_argument('--cmap', type=str, default='seismic',
                        help='Colormap (default: seismic)')
    parser.add_argument('--slice', action='store_true',
                        help='Plot 2D slices instead of 3D volume')
    parser.add_argument('--animate', action='store_true',
                        help='Create animation of all timesteps')
    parser.add_argument('--global-scale', action='store_true',
                        help='Use global min/max for all timesteps (default: per-timestep scale)')

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Read data
    data_list, timestep_list = read_adios2_data(args.file, args.var, args.nx, args.ny, args.nz)

    if data_list is None or len(data_list) == 0:
        print("Error: No data read from file")
        return 1

    # Determine colorscale range
    if args.global_scale:
        # Global min/max for consistent colorscale across all timesteps
        vmin = min(d.min() for d in data_list)
        vmax = max(d.max() for d in data_list)
        print(f"\nUsing global data range: [{vmin:.6f}, {vmax:.6f}]")
    else:
        # Per-timestep scaling
        vmin = None
        vmax = None
        print(f"\nUsing per-timestep scaling (each plot has its own min/max)")

    # Create animation if requested (always uses global scale for consistency)
    if args.animate:
        anim_file = os.path.join(args.output, 'pressure_animation.gif')
        mode = 'slice' if args.slice else 'volume'
        # Animation needs global scale for consistency
        anim_vmin = min(d.min() for d in data_list)
        anim_vmax = max(d.max() for d in data_list)
        create_animation(data_list, timestep_list, anim_file, mode=mode, cmap=args.cmap,
                        vmin=anim_vmin, vmax=anim_vmax)

    # Create individual plots
    print(f"\nCreating plots in {args.output}/...")

    for i, (data, timestep) in enumerate(zip(data_list, timestep_list)):
        output_file = os.path.join(args.output, f'pressure_step_{timestep:05d}.png')

        if args.slice:
            plot_2d_slices(data, timestep, output_file, cmap=args.cmap, vmin=vmin, vmax=vmax)
        else:
            plot_3d_volume(data, timestep, output_file, cmap=args.cmap, vmin=vmin, vmax=vmax)

    print(f"\nDone! Created {len(data_list)} plots in {args.output}/")

    return 0


if __name__ == '__main__':
    sys.exit(main())
