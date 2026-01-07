#!/usr/bin/env python3
"""
Wavefield Snapshot Visualization Tool

This script provides various visualization options for FD solver wavefield snapshots:
- Display individual snapshots
- Extract and display 2D slices (vertical/horizontal planes)
- Create animations from snapshot sequences
- Support for custom snapshot ranges and slice selections

Usage:
    # Display a single snapshot (central vertical plane)
    python3 visualize_snapshots.py --snapshot 100

    # Display a range of snapshots
    python3 visualize_snapshots.py --range 0 200 --step 20

    # Create animation
    python3 visualize_snapshots.py --animate --range 0 890 --step 10 --fps 10

    # Custom slice (XZ plane at y=80)
    python3 visualize_snapshots.py --snapshot 100 --slice xz --slice-index 80

    # Horizontal slice (XY plane at z=50)
    python3 visualize_snapshots.py --snapshot 100 --slice xy --slice-index 50
"""

import argparse
import os
import sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import TwoSlopeNorm


def find_snapshots(snapshot_dir="snapshots", start=None, end=None, step=None):
    """Find snapshot files in the given range.

    Args:
        snapshot_dir: Directory containing snapshot files
        start: Starting time step (None = first available)
        end: Ending time step (None = last available)
        step: Step between snapshots (None = use all in range)

    Returns:
        List of (time_step, filepath) tuples
    """
    # Find all snapshot files
    pattern = os.path.join(snapshot_dir, "wavefield_*.npy")
    files = sorted(glob(pattern))

    if not files:
        print(f"Error: No snapshot files found in {snapshot_dir}")
        sys.exit(1)

    # Extract time steps from filenames
    snapshots = []
    for filepath in files:
        basename = os.path.basename(filepath)
        # Extract number from wavefield_NNNNNN.npy
        time_step = int(basename.split('_')[1].split('.')[0])
        snapshots.append((time_step, filepath))

    # Sort by time step
    snapshots.sort(key=lambda x: x[0])

    # Apply range filtering
    if start is not None:
        snapshots = [(t, f) for t, f in snapshots if t >= start]
    if end is not None:
        snapshots = [(t, f) for t, f in snapshots if t <= end]

    # Apply step filtering
    if step is not None and step > 1:
        filtered = []
        for i, (t, f) in enumerate(snapshots):
            if i % step == 0:
                filtered.append((t, f))
        snapshots = filtered

    return snapshots


def extract_slice(wavefield, slice_type='xz', slice_index=None):
    """Extract a 2D slice from the 3D wavefield.

    Args:
        wavefield: 3D numpy array (nx, ny, nz)
        slice_type: Type of slice ('xy', 'xz', 'yz')
        slice_index: Index along the perpendicular axis (None = central slice)

    Returns:
        2D numpy array containing the slice
        Tuple of (xlabel, ylabel) for the plot
    """
    nx, ny, nz = wavefield.shape

    if slice_type == 'xy':
        # Horizontal slice at constant z
        if slice_index is None:
            slice_index = nz // 2
        if slice_index < 0 or slice_index >= nz:
            slice_index = nz // 2
            print(f"Warning: Invalid z-index, using central slice at z={slice_index}")
        slice_2d = wavefield[:, :, slice_index]
        labels = ('X index', 'Y index')
        title_suffix = f'XY plane at z={slice_index}'

    elif slice_type == 'xz':
        # Vertical slice at constant y (central vertical plane)
        if slice_index is None:
            slice_index = ny // 2
        if slice_index < 0 or slice_index >= ny:
            slice_index = ny // 2
            print(f"Warning: Invalid y-index, using central slice at y={slice_index}")
        slice_2d = wavefield[:, slice_index, :]
        labels = ('X index', 'Z index')
        title_suffix = f'XZ plane at y={slice_index}'

    elif slice_type == 'yz':
        # Vertical slice at constant x
        if slice_index is None:
            slice_index = nx // 2
        if slice_index < 0 or slice_index >= nx:
            slice_index = nx // 2
            print(f"Warning: Invalid x-index, using central slice at x={slice_index}")
        slice_2d = wavefield[slice_index, :, :]
        labels = ('Y index', 'Z index')
        title_suffix = f'YZ plane at x={slice_index}'

    else:
        raise ValueError(f"Unknown slice type: {slice_type}. Use 'xy', 'xz', or 'yz'")

    return slice_2d, labels, title_suffix


def display_snapshot(snapshot_file, slice_type='xz', slice_index=None,
                     vmin=None, vmax=None, cmap='seismic', output=None):
    """Display a single wavefield snapshot.

    Args:
        snapshot_file: Path to the .npy file
        slice_type: Type of slice to display ('xy', 'xz', 'yz')
        slice_index: Index for the slice (None = central)
        vmin, vmax: Color scale limits (None = auto)
        cmap: Colormap name
        output: Output filename (None = display on screen)
    """
    # Load wavefield
    wavefield = np.load(snapshot_file)
    time_step = int(os.path.basename(snapshot_file).split('_')[1].split('.')[0])

    # Extract slice
    slice_2d, labels, title_suffix = extract_slice(wavefield, slice_type, slice_index)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Determine color scale
    if vmin is None or vmax is None:
        data_max = np.abs(slice_2d).max()
        if vmin is None:
            vmin = -data_max
        if vmax is None:
            vmax = data_max

    # Use symmetric colormap centered at zero
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    # Plot
    im = ax.imshow(slice_2d.T, cmap=cmap, norm=norm,
                   aspect='auto', origin='lower', interpolation='bilinear')

    # Labels and title
    ax.set_xlabel(labels[0], fontsize=12)
    ax.set_ylabel(labels[1], fontsize=12)
    ax.set_title(f'Wavefield at time step {time_step}\n{title_suffix}', fontsize=14)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Pressure', fontsize=12)

    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')

    # Statistics
    stats_text = (f'Min: {slice_2d.min():.2e}\n'
                 f'Max: {slice_2d.max():.2e}\n'
                 f'Mean: {slice_2d.mean():.2e}')
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round',
            facecolor='white', alpha=0.8), fontsize=10)

    plt.tight_layout()

    if output:
        plt.savefig(output, dpi=150, bbox_inches='tight')
        print(f"Saved: {output}")
    else:
        plt.show()

    plt.close()


def display_snapshot_range(snapshots, slice_type='xz', slice_index=None,
                           vmin=None, vmax=None, cmap='seismic',
                           output_dir='snapshot_images'):
    """Display a range of snapshots as individual images.

    Args:
        snapshots: List of (time_step, filepath) tuples
        slice_type: Type of slice to display
        slice_index: Index for the slice
        vmin, vmax: Color scale limits
        cmap: Colormap name
        output_dir: Directory to save images
    """
    os.makedirs(output_dir, exist_ok=True)

    # Determine global color scale if not provided
    if vmin is None or vmax is None:
        print("Computing global color scale...")
        data_max = 0
        for time_step, filepath in snapshots:
            wavefield = np.load(filepath)
            slice_2d, _, _ = extract_slice(wavefield, slice_type, slice_index)
            data_max = max(data_max, np.abs(slice_2d).max())

        if vmin is None:
            vmin = -data_max
        if vmax is None:
            vmax = data_max
        print(f"Color scale: [{vmin:.2e}, {vmax:.2e}]")

    # Process each snapshot
    print(f"Processing {len(snapshots)} snapshots...")
    for i, (time_step, filepath) in enumerate(snapshots):
        output = os.path.join(output_dir, f'snapshot_{time_step:06d}.png')
        display_snapshot(filepath, slice_type, slice_index, vmin, vmax, cmap, output)

        if (i + 1) % 10 == 0 or i == len(snapshots) - 1:
            print(f"  Processed {i + 1}/{len(snapshots)}")

    print(f"All snapshots saved to: {output_dir}")


def create_animation(snapshots, slice_type='xz', slice_index=None,
                    vmin=None, vmax=None, cmap='seismic',
                    fps=10, output='wavefield_animation.mp4', dpi=150):
    """Create an animation from a sequence of snapshots.

    Args:
        snapshots: List of (time_step, filepath) tuples
        slice_type: Type of slice to display
        slice_index: Index for the slice
        vmin, vmax: Color scale limits
        cmap: Colormap name
        fps: Frames per second
        output: Output filename
        dpi: Resolution (dots per inch)
    """
    print(f"Creating animation with {len(snapshots)} frames...")

    # Determine slice parameters from first snapshot
    first_wavefield = np.load(snapshots[0][1])
    slice_2d, labels, title_suffix_base = extract_slice(first_wavefield, slice_type, slice_index)

    # Determine global color scale if not provided
    if vmin is None or vmax is None:
        print("Computing global color scale...")
        data_max = 0
        for time_step, filepath in snapshots:
            wavefield = np.load(filepath)
            slice_2d_temp, _, _ = extract_slice(wavefield, slice_type, slice_index)
            data_max = max(data_max, np.abs(slice_2d_temp).max())

        if vmin is None:
            vmin = -data_max
        if vmax is None:
            vmax = data_max
        print(f"Color scale: [{vmin:.2e}, {vmax:.2e}]")

    # Use symmetric colormap centered at zero
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    # Setup figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Initial plot
    im = ax.imshow(slice_2d.T, cmap=cmap, norm=norm,
                   aspect='auto', origin='lower', interpolation='bilinear')

    ax.set_xlabel(labels[0], fontsize=12)
    ax.set_ylabel(labels[1], fontsize=12)
    title = ax.set_title('', fontsize=14)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Pressure', fontsize=12)

    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')

    # Statistics text box
    stats_text = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round',
                        facecolor='white', alpha=0.8), fontsize=10)

    plt.tight_layout()

    # Animation update function
    def update(frame):
        time_step, filepath = snapshots[frame]
        wavefield = np.load(filepath)
        slice_2d, _, title_suffix = extract_slice(wavefield, slice_type, slice_index)

        # Update image data
        im.set_data(slice_2d.T)

        # Update title
        title.set_text(f'Wavefield at time step {time_step}\n{title_suffix}')

        # Update statistics
        stats = (f'Min: {slice_2d.min():.2e}\n'
                f'Max: {slice_2d.max():.2e}\n'
                f'Mean: {slice_2d.mean():.2e}')
        stats_text.set_text(stats)

        # Progress
        if (frame + 1) % 10 == 0 or frame == len(snapshots) - 1:
            print(f"  Frame {frame + 1}/{len(snapshots)}")

        return [im, title, stats_text]

    # Create animation
    anim = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=1000/fps, blit=True)

    # Save animation
    print(f"Saving animation to: {output}")
    try:
        writer = animation.FFMpegWriter(fps=fps, bitrate=2000)
        anim.save(output, writer=writer, dpi=dpi)
        print(f"Animation saved successfully!")
        print(f"  Frames: {len(snapshots)}")
        print(f"  FPS: {fps}")
        print(f"  Duration: {len(snapshots)/fps:.1f} seconds")
    except Exception as e:
        print(f"Error saving animation: {e}")
        print("Note: FFmpeg must be installed for animation export")
        print("  Install with: sudo apt-get install ffmpeg")
    finally:
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize FD solver wavefield snapshots',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input options
    parser.add_argument('--dir', default='snapshots',
                       help='Directory containing snapshot files (default: snapshots)')
    parser.add_argument('--snapshot', type=int,
                       help='Display a single snapshot at this time step')
    parser.add_argument('--range', nargs=2, type=int, metavar=('START', 'END'),
                       help='Range of time steps to process (inclusive)')
    parser.add_argument('--step', type=int, default=1,
                       help='Step between snapshots in range (default: 1)')

    # Slice options
    parser.add_argument('--slice', choices=['xy', 'xz', 'yz'], default='xz',
                       help='Type of 2D slice (default: xz - central vertical plane)')
    parser.add_argument('--slice-index', type=int,
                       help='Index for slice along perpendicular axis (default: central)')

    # Visualization options
    parser.add_argument('--cmap', default='seismic',
                       help='Colormap name (default: seismic)')
    parser.add_argument('--vmin', type=float,
                       help='Minimum value for color scale (default: auto)')
    parser.add_argument('--vmax', type=float,
                       help='Maximum value for color scale (default: auto)')

    # Output options
    parser.add_argument('--animate', action='store_true',
                       help='Create animation instead of individual images')
    parser.add_argument('--fps', type=int, default=10,
                       help='Frames per second for animation (default: 10)')
    parser.add_argument('--output',
                       help='Output filename (for single snapshot or animation)')
    parser.add_argument('--output-dir', default='snapshot_images',
                       help='Output directory for image range (default: snapshot_images)')
    parser.add_argument('--dpi', type=int, default=150,
                       help='Resolution in DPI (default: 150)')

    args = parser.parse_args()

    # Validate options
    if args.snapshot is None and args.range is None:
        parser.error("Must specify either --snapshot or --range")

    if args.snapshot is not None and args.range is not None:
        parser.error("Cannot specify both --snapshot and --range")

    # Process single snapshot
    if args.snapshot is not None:
        snapshots = find_snapshots(args.dir)
        # Find the snapshot file for the requested time step
        snapshot_file = None
        for time_step, filepath in snapshots:
            if time_step == args.snapshot:
                snapshot_file = filepath
                break

        if snapshot_file is None:
            available = [t for t, f in snapshots]
            print(f"Error: Snapshot at time step {args.snapshot} not found")
            print(f"Available time steps: {available[:10]}{'...' if len(available) > 10 else ''}")
            sys.exit(1)

        output = args.output
        if output is None and not args.animate:
            output = f'snapshot_{args.snapshot:06d}.png'

        display_snapshot(snapshot_file, args.slice, args.slice_index,
                        args.vmin, args.vmax, args.cmap, output)

    # Process range of snapshots
    else:
        start, end = args.range
        snapshots = find_snapshots(args.dir, start, end, args.step)

        if not snapshots:
            print(f"Error: No snapshots found in range [{start}, {end}] with step {args.step}")
            sys.exit(1)

        print(f"Found {len(snapshots)} snapshots in range [{start}, {end}]")

        if args.animate:
            output = args.output if args.output else 'wavefield_animation.mp4'
            create_animation(snapshots, args.slice, args.slice_index,
                           args.vmin, args.vmax, args.cmap,
                           args.fps, output, args.dpi)
        else:
            display_snapshot_range(snapshots, args.slice, args.slice_index,
                                  args.vmin, args.vmax, args.cmap,
                                  args.output_dir)


if __name__ == '__main__':
    main()
