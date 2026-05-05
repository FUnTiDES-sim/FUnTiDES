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
    Basic usage for visualizing a 201x201x201 grid:
    python visualize_pressure.py 201 201 201

    Specify a custom BP file:
    python visualize_pressure.py 201 201 201 --file my_data.bp

    Create 2D slices which is significantly faster for large datasets:
    python visualize_pressure.py 201 201 201 --slice

    Create an animation covering all timesteps:
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
    Reads time-series data from an ADIOS2 BP file.

    This function opens an ADIOS2 stream, iterates over all available steps,
    and extracts the specified variable. It ensures that the array dimensions
    match the expected grid size and handles matrix transpositions if the C++
    solver exported the data in a different memory layout.

    Parameters:
    -----------
    filename : str
        Path to the ADIOS2 BP file or directory.
    variable_name : str
        Name of the variable to extract from the stream.
    nx, ny, nz : int
        Expected grid dimensions in the X, Y, and Z directions.

    Returns:
    --------
    data_list : list of numpy arrays
        List containing the 3D data array for each extracted timestep.
    timestep_list : list of int
        List containing the corresponding timestep indices.
    """
    print(f"Reading {filename}...")

    if not os.path.exists(filename):
        print(f"Error: File not found: {filename}")
        return None, None

    data_list = []
    timestep_list = []

    try:
        import adios2

        print("  Reading via the official adios2 library...")

        # Utilizing the updated Python API for ADIOS2 streams
        with adios2.Stream(filename, "r") as s:
            for _ in s.steps():
                step_idx = s.current_step()

                try:
                    # Extract the array for the current simulation step
                    data = s.read(variable_name)
                except Exception as e:
                    print(
                        f"  [Warning] Unable to read '{variable_name}' at step {step_idx}: {e}"
                    )
                    continue

                # ADIOS2 preserves the native C++ memory order.
                # If the shape indicates a reversed axis order, transpose it
                # to match the expected (nx, ny, nz) layout in Python.
                if data.shape == (nz, ny, nx):
                    data = data.transpose(2, 1, 0)
                elif data.shape == (ny, nz, nx):
                    data = data.transpose(2, 0, 1)

                data_list.append(data.astype(np.float32))
                timestep_list.append(step_idx)

        print(f"  Successfully read {len(data_list)} timesteps")
        return data_list, timestep_list

    except ImportError:
        print("\n  [ERROR] The Python module 'adios2' is not installed.")
        print("  Please install it by running: pip install adios2\n")
        return None, None
    except Exception as e:
        print(f"  Unexpected error during reading with adios2: {e}")
        import traceback

        traceback.print_exc()
        return None, None


def read_bp5_directory(directory, variable_name, nx, ny, nz):
    """
    Reads raw data from a BP5 format directory structure as a fallback.

    BP5 files are typically stored as directories containing data.0, md.0,
    md.idx, and mmd.0 files. This function reads the raw binary payload from
    data.0 and reshapes it based on the expected grid dimensions.

    Parameters:
    -----------
    directory : str
        Path to the BP5 directory.
    variable_name : str
        Name of the variable to read.
    nx, ny, nz : int
        Expected grid dimensions.

    Returns:
    --------
    data_list : list of numpy arrays
        List of 3D arrays, one per timestep.
    timestep_list : list of int
        Timestep values corresponding to the data arrays.
    """
    import struct

    data_list = []
    timestep_list = []

    if not os.path.isdir(directory):
        print(f"  Error: {directory} is not a directory")
        return None, None

    try:
        # Verify the existence of the core BP5 structure files
        data_file = os.path.join(directory, "data.0")
        md_file = os.path.join(directory, "md.0")

        if not os.path.exists(data_file):
            print(f"  Error: {data_file} not found")
            return None, None

        print(f"  Found BP5 files: data.0, md.0, etc.")

        # Read the raw binary data file entirely into memory
        with open(data_file, "rb") as f:
            raw_data = f.read()

        # Convert the raw bytes into a flat float32 array
        num_floats = len(raw_data) // 4
        data_array = np.frombuffer(raw_data, dtype=np.float32, count=num_floats)

        print(f"  Read {num_floats} float32 values from data.0")
        print(f"  Expected per timestep: {nx * ny * nz}")

        # Calculate the total number of complete timesteps available in the buffer
        values_per_step = nx * ny * nz
        num_steps = num_floats // values_per_step
        remainder = num_floats % values_per_step

        print(
            f"  Detected {num_steps} complete timesteps + {remainder} remainder values"
        )

        if num_steps > 0:
            for step in range(num_steps):
                start_idx = step * values_per_step
                end_idx = start_idx + values_per_step

                # Reshape and adjust axes to match the Python indexing convention
                step_data = (
                    data_array[start_idx:end_idx]
                    .reshape((nz, ny, nx))
                    .transpose(2, 1, 0)
                )
                data_list.append(step_data)
                timestep_list.append(step)
            print(f"  Successfully read {num_steps} complete timesteps")
        else:
            print(f"  Error: Not enough data for even one timestep")
            return None, None

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


def plot_3d_volume(data, timestep, output_file, cmap="seismic", vmin=None, vmax=None):
    """
    Creates a 3D volume visualization using a scattered point cloud representation.

    To maintain performance, the function heavily subsamples the data volume
    before plotting, mapping the scalar field values to a specified colormap.

    Parameters:
    -----------
    data : numpy array (3D)
        The 3D array containing the scalar field values.
    timestep : int
        The current timestep index for titling.
    output_file : str
        Path where the resulting image will be saved.
    cmap : str
        The Matplotlib colormap to apply.
    vmin, vmax : float
        Explicit minimum and maximum values for the colormap normalization.
    """
    nx, ny, nz = data.shape

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Generate 3D grid coordinates corresponding to the data array
    x, y, z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")

    # Determine an appropriate stride to avoid overwhelming the rendering engine
    stride = max(1, min(nx, ny, nz) // 50)

    x_sub = x[::stride, ::stride, ::stride]
    y_sub = y[::stride, ::stride, ::stride]
    z_sub = z[::stride, ::stride, ::stride]
    data_sub = data[::stride, ::stride, ::stride]

    # Map the scalar data to the scatter plot points
    scatter = ax.scatter(
        x_sub.ravel(),
        y_sub.ravel(),
        z_sub.ravel(),
        c=data_sub.ravel(),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        alpha=0.6,
        s=20,
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(f"Pressure Field - Timestep {timestep}")

    plt.colorbar(scatter, ax=ax, label="Pressure", shrink=0.8)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_file}")


def _compute_vmax(data, percentile=100):
    """
    Computes a symmetric maximum value for colormap scaling based on percentiles.
    This helps in dampening extreme outlier values for better visual contrast.
    """
    vmax = float(np.percentile(np.abs(data), percentile))
    return vmax if vmax > 0 else 1.0


def plot_2d_slices(
    data, timestep, output_file, cmap="seismic", vmin=None, vmax=None, percentile=100
):
    """
    Generates orthogonal 2D slice views through the central axes of the domain.

    This function extracts the XY, XZ, and YZ planes intersecting at the exact
    center of the dataset and plots them side-by-side with general statistics.

    Parameters:
    -----------
    data : numpy array (3D)
        The 3D array representing the pressure field.
    timestep : int
        The temporal step index of the data.
    output_file : str
        Target file path for saving the figure.
    cmap : str
        The colormap to map scalar values.
    vmin, vmax : float
        Boundaries for the color scale. Overridden by percentile if left empty.
    percentile : int
        Clips the color scale to ignore extreme outliers.
    """
    nx, ny, nz = data.shape

    if vmin is None or vmax is None:
        _vmax = _compute_vmax(data, percentile)
        vmin = -_vmax
        vmax = _vmax

    # Locate the central indices of the volumetric domain
    mid_x = nx // 2
    mid_y = ny // 2
    mid_z = nz // 2

    # Extract planar cross-sections
    slice_yz = data[mid_x, :, :]
    slice_xz = data[:, mid_y, :]
    slice_xy = data[:, :, mid_z]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Render the XY cross-section
    im1 = axes[0, 0].imshow(slice_xy.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 0].set_title(f"XY Slice (Z={mid_z})")
    axes[0, 0].set_xlabel("X")
    axes[0, 0].set_ylabel("Y")
    plt.colorbar(im1, ax=axes[0, 0])

    # Render the XZ cross-section
    im2 = axes[0, 1].imshow(slice_xz.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    axes[0, 1].set_title(f"XZ Slice (Y={mid_y})")
    axes[0, 1].set_xlabel("X")
    axes[0, 1].set_ylabel("Z")
    plt.colorbar(im2, ax=axes[0, 1])

    # Render the YZ cross-section
    im3 = axes[1, 0].imshow(slice_yz.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    axes[1, 0].set_title(f"YZ Slice (X={mid_x})")
    axes[1, 0].set_xlabel("Y")
    axes[1, 0].set_ylabel("Z")
    plt.colorbar(im3, ax=axes[1, 0])

    # Display dataset statistics in the remaining quadrant
    axes[1, 1].axis("off")
    stats_text = f"""
    Timestep: {timestep}
    Grid: {nx} x {ny} x {nz}

    Statistics:
    Min:  {data.min():.6f}
    Max:  {data.max():.6f}
    Mean: {data.mean():.6f}
    Std:  {data.std():.6f}
    """
    axes[1, 1].text(
        0.1,
        0.5,
        stats_text,
        fontsize=12,
        family="monospace",
        verticalalignment="center",
    )

    fig.suptitle(
        f"Pressure Field 2D Slices - Timestep {timestep}", fontsize=14, y=0.995
    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_file}")


def plot_2d_slices_coupled(
    pn_data, uz_data, timestep, output_file, cmap="seismic", percentile=100
):
    """
    Generates side-by-side 2D slices comparing acousto-elastic coupled fields.

    The visualization uses a grid format. The top row visualizes the acoustic
    pressure field while the bottom row visualizes vertical elastic displacement.
    Independent color ranges are calculated to ensure both fields are visible.

    Parameters:
    -----------
    pn_data : numpy array (3D)
        Acoustic pressure field data.
    uz_data : numpy array (3D)
        Vertical elastic displacement field data.
    timestep : int
        Simulation timestep index.
    output_file : str
        Target file path for the generated image.
    cmap : str
        Colormap to apply to both fields.
    percentile : int
        Used to clip visual extremes in the color mapping process.
    """
    nx, ny, nz = pn_data.shape
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    def _add_row(row, data, label, unit):
        vmax = _compute_vmax(data, percentile)
        vmin = -vmax

        slices = [
            (data[:, :, mid_z], f"{label} XY Slice (Z={mid_z})", "X", "Y"),
            (data[:, mid_y, :], f"{label} XZ Slice (Y={mid_y})", "X", "Z"),
            (data[mid_x, :, :], f"{label} YZ Slice (X={mid_x})", "Y", "Z"),
        ]
        for col, (sl, title, xlabel, ylabel) in enumerate(slices):
            im = axes[row, col].imshow(
                sl.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax
            )
            axes[row, col].set_title(title)
            axes[row, col].set_xlabel(xlabel)
            axes[row, col].set_ylabel(ylabel)
            plt.colorbar(im, ax=axes[row, col], label=unit)

    _add_row(0, pn_data, "Pressure p", "Pa")
    _add_row(1, uz_data, "Displacement uz", "m")

    fig.suptitle(
        f"Acousto-Elastic Fields 2D Slices - Timestep {timestep}", fontsize=14, y=0.995
    )
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_file}")


def create_animation_coupled(
    pn_list, uz_list, timestep_list, output_file, cmap="seismic", percentile=100
):
    """
    Creates an animated GIF highlighting acousto-elastic coupled fields over time.

    Frames are composed of pressure slices (top row) and vertical displacement
    slices (bottom row). The color boundaries remain fixed across the entire
    animation to prevent flickering and maintain consistent physical scaling.

    Parameters:
    -----------
    pn_list : list of numpy arrays
        Time-series list of the acoustic pressure field.
    uz_list : list of numpy arrays
        Time-series list of the elastic displacement field.
    timestep_list : list of int
        Corresponding timeline indices.
    output_file : str
        Target file path for the animated GIF.
    cmap : str
        Colormap name.
    percentile : int
        Threshold to exclude noise spikes from the color scaling calculation.
    """
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
    except ImportError:
        print("Error: Animation requires pillow. Install with: pip install pillow")
        return

    n_frames = len(pn_list)
    print(f"Creating coupled animation with {n_frames} frames...")

    # Calculate global color boundaries across all timesteps to prevent visual jitter
    all_pn = np.concatenate([d.ravel() for d in pn_list])
    all_uz = np.concatenate([d.ravel() for d in uz_list])
    pn_vmax = _compute_vmax(all_pn, percentile)
    uz_vmax = _compute_vmax(all_uz, percentile)

    print(f"  Pressure range:     [{-pn_vmax:.4e}, {pn_vmax:.4e}] Pa")
    print(f"  Displacement range: [{-uz_vmax:.4e}, {uz_vmax:.4e}] m")

    nx, ny, nz = pn_list[0].shape
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    def _slices(data, mid_x, mid_y, mid_z):
        return [data[:, :, mid_z].T, data[:, mid_y, :].T, data[mid_x, :, :].T]

    # Pre-configure image layers to utilize efficient blitting during animation
    ims = []
    for col, (sl, xlabel, ylabel, title_tpl) in enumerate(
        [
            (pn_list[0][:, :, mid_z].T, "X", "Y", "{label} XY (Z={mid_z})"),
            (pn_list[0][:, mid_y, :].T, "X", "Z", "{label} XZ (Y={mid_y})"),
            (pn_list[0][mid_x, :, :].T, "Y", "Z", "{label} YZ (X={mid_x})"),
        ]
    ):
        for row, (vmax, label) in enumerate([(pn_vmax, "p"), (uz_vmax, "uz")]):
            im = axes[row, col].imshow(
                sl, origin="lower", cmap=cmap, vmin=-vmax, vmax=vmax, animated=True
            )
            axes[row, col].set_xlabel(xlabel)
            axes[row, col].set_ylabel(ylabel)
            plt.colorbar(im, ax=axes[row, col], label="Pa" if row == 0 else "m")
            ims.append(im)

    title_obj = fig.suptitle("", fontsize=14)

    def update(frame):
        pn = pn_list[frame]
        uz = uz_list[frame]
        timestep = timestep_list[frame]

        pn_slices = _slices(pn, mid_x, mid_y, mid_z)
        uz_slices = _slices(uz, mid_x, mid_y, mid_z)

        updated = []
        idx = 0
        for col in range(3):
            ims[idx].set_data(pn_slices[col])
            ims[idx + 1].set_data(uz_slices[col])
            updated += [ims[idx], ims[idx + 1]]
            idx += 2

        title_obj.set_text(f"Acousto-Elastic Fields - Timestep {timestep}")
        updated.append(title_obj)
        return updated

    anim = FuncAnimation(fig, update, frames=n_frames, interval=200, blit=False)
    writer = PillowWriter(fps=5)
    anim.save(output_file, writer=writer)
    plt.close()
    print(f"Animation saved: {output_file}")


def create_animation(
    data_list,
    timestep_list,
    output_file,
    mode="slice",
    cmap="seismic",
    vmin=None,
    vmax=None,
):
    """
    Assembles a time-series animation of the standard pressure field.

    Parameters:
    -----------
    data_list : list of numpy arrays
        Sequential list of scalar datasets.
    timestep_list : list of int
        Corresponding timestep values.
    output_file : str
        Target file path for saving the animation.
    mode : str
        Determines the visual format ('slice' or 'volume').
    cmap : str
        Colormap name for spatial distribution.
    vmin, vmax : float
        Absolute limits defining the boundaries of the colormap.
    """
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
    except ImportError:
        print("Error: Animation requires pillow. Install with: pip install pillow")
        return

    print(f"Creating animation with {len(data_list)} frames...")
    print(f"Animation range: [{vmin:.4f}, {vmax:.4f}]")

    if mode == "slice":
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

            im1 = axes[0].imshow(
                slice_xy.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax
            )
            axes[0].set_title(f"XY Slice (Z={mid_z})")
            axes[0].set_xlabel("X")
            axes[0].set_ylabel("Y")

            im2 = axes[1].imshow(
                slice_xz.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax
            )
            axes[1].set_title(f"XZ Slice (Y={mid_y})")
            axes[1].set_xlabel("X")
            axes[1].set_ylabel("Z")

            im3 = axes[2].imshow(
                slice_yz.T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax
            )
            axes[2].set_title(f"YZ Slice (X={mid_x})")
            axes[2].set_xlabel("Y")
            axes[2].set_ylabel("Z")

            fig.suptitle(f"Timestep {timestep}", fontsize=14)

            return [im1, im2, im3]

        anim = FuncAnimation(
            fig, update, frames=len(data_list), interval=200, blit=False
        )

        writer = PillowWriter(fps=5)
        anim.save(output_file, writer=writer)

        plt.close()
        print(f"Animation saved: {output_file}")

    else:
        print("3D volume animation not implemented. Use --slice mode.")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize 3D pressure field from ADIOS2 BP files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("nx", type=int, help="Number of nodes in X direction")
    parser.add_argument("ny", type=int, help="Number of nodes in Y direction")
    parser.add_argument("nz", type=int, help="Number of nodes in Z direction")
    parser.add_argument(
        "--file",
        type=str,
        default="snapshots",
        help="ADIOS2 BP file path (default: snapshots)",
    )
    parser.add_argument(
        "--var",
        type=str,
        default="PressureField",
        help="Variable name to visualize (default: PressureField)",
    )
    parser.add_argument(
        "--output", type=str, default="plots", help="Output directory (default: plots/)"
    )
    parser.add_argument(
        "--cmap", type=str, default="seismic", help="Colormap (default: seismic)"
    )
    parser.add_argument(
        "--slice", action="store_true", help="Plot 2D slices instead of 3D volume"
    )
    parser.add_argument(
        "--animate", action="store_true", help="Create animation of all timesteps"
    )
    parser.add_argument(
        "--global-scale",
        action="store_true",
        help="Use global min/max for all timesteps (default: per-timestep scale)",
    )
    parser.add_argument(
        "--coupled",
        action="store_true",
        help="Acousto-elastic mode: ADIOS2 steps alternate pressure/uz; "
        "display both fields side by side (requires --slice)",
    )
    parser.add_argument(
        "--percentile",
        type=int,
        default=100,
        help="Clip color scale at this percentile of abs values "
        "(default: 100 = true max; try 95 or 98 for more contrast)",
    )

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    data_list, timestep_list = read_adios2_data(
        args.file, args.var, args.nx, args.ny, args.nz
    )

    if data_list is None or len(data_list) == 0:
        print("Error: No data read from file")
        return 1

    # Determines the mapping boundaries depending on the user's scale configuration
    if args.global_scale:
        vmin = min(d.min() for d in data_list)
        vmax = max(d.max() for d in data_list)
        print(f"\nUsing global data range: [{vmin:.6f}, {vmax:.6f}]")
    else:
        vmin = None
        vmax = None
        print(f"\nUsing per-timestep scaling (each plot has its own min/max)")

    # Execute animation sequence rendering
    if args.animate:
        if args.coupled:
            pn_list = data_list[0::2]
            uz_list = data_list[1::2]
            ts_list = timestep_list[0::2]
            anim_file = os.path.join(args.output, "coupled_animation.gif")
            create_animation_coupled(
                pn_list,
                uz_list,
                ts_list,
                anim_file,
                cmap=args.cmap,
                percentile=args.percentile,
            )
        else:
            anim_file = os.path.join(args.output, "pressure_animation.gif")
            mode = "slice" if args.slice else "volume"
            anim_vmin = min(d.min() for d in data_list)
            anim_vmax = max(d.max() for d in data_list)
            create_animation(
                data_list,
                timestep_list,
                anim_file,
                mode=mode,
                cmap=args.cmap,
                vmin=anim_vmin,
                vmax=anim_vmax,
            )

    print(f"\nCreating plots in {args.output}/...")

    # Output standard image processing loop
    if args.coupled and args.slice:
        # Alternating iteration step where even index = pressure, odd index = uz displacement
        for i in range(0, len(data_list) - 1, 2):
            pn_data = data_list[i]
            uz_data = data_list[i + 1]
            timestep = timestep_list[i]
            output_file = os.path.join(args.output, f"coupled_step_{timestep:05d}.png")
            plot_2d_slices_coupled(
                pn_data,
                uz_data,
                timestep,
                output_file,
                cmap=args.cmap,
                percentile=args.percentile,
            )
    else:
        for i, (data, timestep) in enumerate(zip(data_list, timestep_list)):
            output_file = os.path.join(args.output, f"pressure_step_{timestep:05d}.png")

            if args.slice:
                plot_2d_slices(
                    data,
                    timestep,
                    output_file,
                    cmap=args.cmap,
                    vmin=vmin,
                    vmax=vmax,
                    percentile=args.percentile,
                )
            else:
                plot_3d_volume(
                    data, timestep, output_file, cmap=args.cmap, vmin=vmin, vmax=vmax
                )

    print(f"\nDone! Created {len(data_list)} plots in {args.output}/")

    return 0


if __name__ == "__main__":
    sys.exit(main())
