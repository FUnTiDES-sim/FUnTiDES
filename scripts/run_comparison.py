#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_comparison.py
=================
Launch SEM, DG, and DG-SEM acoustic simulations on a 2000×2000×2000 m³ cube
and generate comparison plots.

What the script produces (inside --outdir, default: "comparison_results"):
  sem/    dg/    dgsem/   — raw simulation outputs (receiver_trace.txt, slice*.dat, …)
  plots/
    receiver_sem.png          — SEM receiver pressure vs. time
    receiver_dg.png           — DG  receiver pressure vs. time
    receiver_dgsem.png        — DG-SEM receiver pressure vs. time
    compare_dg_sem.png        — DG + SEM overlay
    compare_dgsem_sem.png     — DG-SEM + SEM overlay
    snapshots_sem/            — SEM cross-section (XY @ z=1000 m) at every snapshot step
    snapshots_dg/             — DG  cross-section (XY @ z≈1000 m) at every snapshot step
    snapshots_dgsem/          — DG-SEM combined XZ cross-section at every snapshot step

Usage:
  python run_comparison.py [--outdir DIR] [--skip-sim] [--no-snapshots] [--np N]

Options:
  --outdir DIR        Output root directory (default: comparison_results)
  --skip-sim          Skip simulations (only regenerate plots from existing output)
  --no-snapshots      Disable snapshot saving (faster runs, no wavefield cross-section plots)
  --np N              Number of MPI processes (default: 1)
  --snap-interval N   Iterations between snapshots (default: 100; 50 snapshots over 5000 steps)
"""
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — safe on HPC
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

ORDER   = 3          # polynomial approximation order
EX = EY = EZ = 100  # elements per direction  (h = 20 m → ~10 pts/wavelength at 7 Hz)
LX = LY = LZ = 2000.0   # domain size in metres

DT      = 0.0003    # time step (s)
TIMEMAX = 1.5       # total simulation time (s)
F0      = 7.0       # dominant Ricker frequency (Hz)
TPEAK   = 1.0 / F0 # Ricker peak time (s)
RICKER_ORDER = 2    # second-order Ricker wavelet

SRCX, SRCY, SRCZ = 1000.0, 1000.0,  500.0
RCVX, RCVY, RCVZ = 1000.0, 1000.0, 1500.0

# DG-SEM coupling interface (source in DG zone, receiver in SEM zone)
DGSEM_BOUNDARY_Z = 1000.0

# Absorbing boundary — size=0 uses only Clayton-Engquist ABC (same as DG), no sponge
BOUNDARIES_SIZE = 0.0
TAPER_DELTA     = 100.0

# Snapshot saving
DEFAULT_SNAP_INTERVAL = 125  # save every 125 steps → 40 snapshots over 5000 steps

# Derived mesh sizes
N_STEPS = int(round(TIMEMAX / DT))   # 5000
N1D     = ORDER + 1                  # GLL points per element per direction
NX = NY = NZ = EX * ORDER + 1        # number of nodes per direction (151)

EZ_DG         = EZ // 2             # DG elements in z for DG-SEM  (25)
NZ_DG_NODES   = EZ_DG * ORDER       # DG nodal z rows  (75)

# ---------------------------------------------------------------------------
# Binary auto-detection
# ---------------------------------------------------------------------------

BINARY_NAME = "funtides-sem"
PROJECT_ROOT = Path(__file__).resolve().parent.parent  # scripts/../


def find_executable(root: Path) -> Path:
    candidates = sorted(
        root.glob(f"build*/bin/{BINARY_NAME}"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No '{BINARY_NAME}' found in any build*/bin/ under {root}.\n"
            "Build the project first (cmake + make)."
        )
    return candidates[0]


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------

def build_cli_args(method: str, snap_interval: int, use_snapshots: bool) -> list[str]:
    """Build funtides-sem CLI argument list for the requested method."""
    args = [
        # Domain
        "--order", str(ORDER),
        "--ex", str(EX), "--ey", str(EY), "--ez", str(EZ),
        "--lx", str(LX), "--ly", str(LY), "--lz", str(LZ),
        "--mesh", "cartesian",
        "--implem", "makutu",
        "--method", method,
        # Time
        "--dt", str(DT),
        "--timemax", str(TIMEMAX),
        # Source
        "--f0", str(F0),
        "--tpeak", str(TPEAK),
        "--ricker-order", str(RICKER_ORDER),
        "--srcx", str(SRCX), "--srcy", str(SRCY), "--srcz", str(SRCZ),
        # Receiver
        "--rcvx", str(RCVX), "--rcvy", str(RCVY), "--rcvz", str(RCVZ),
        # Absorbing boundaries
        "--boundaries-size", str(BOUNDARIES_SIZE),
        "--taper-delta", str(TAPER_DELTA),
    ]

    if use_snapshots:
        args += ["--snapshots", "--snap-interval", str(snap_interval)]

    if method == "dg-sem":
        args += ["--dg-sem-boundary-z", str(DGSEM_BOUNDARY_Z)]

    return args


def run_simulation(binary: Path, method: str, outdir: Path,
                   snap_interval: int, use_snapshots: bool, np_procs: int) -> int:
    """Run funtides-sem for *method* with outputs going to *outdir*."""
    outdir.mkdir(parents=True, exist_ok=True)

    cli = build_cli_args(method, snap_interval, use_snapshots)

    if np_procs > 1:
        cmd = ["srun", "-n", str(np_procs), str(binary)] + cli
    else:
        cmd = [str(binary)] + cli

    print(f"\n{'='*60}")
    print(f"  Running {method.upper()} simulation")
    print(f"  Output dir : {outdir}")
    print(f"  Command    : {' '.join(cmd)}")
    print(f"{'='*60}\n", flush=True)

    result = subprocess.run(cmd, cwd=outdir)
    return result.returncode


# ---------------------------------------------------------------------------
# Output readers
# ---------------------------------------------------------------------------

def read_receiver_trace(dirpath: Path) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """Parse receiver_trace.txt → (time_array, pressure_array)."""
    fpath = dirpath / "receiver_trace.txt"
    if not fpath.exists():
        print(f"  [WARN] {fpath} not found — skipping")
        return None, None
    data = []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                data.append((float(parts[0]), float(parts[1])))
    if not data:
        return None, None
    arr = np.array(data, dtype=np.float32)
    return arr[:, 0], arr[:, 1]


def read_dat_slice(filepath: Path) -> np.ndarray | None:
    """
    Load a whitespace-delimited 2-D slice saved by funtides-sem.
    Returns a (nrows, ncols) float32 array, or None on error.
    """
    if not filepath.exists():
        return None
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append([float(v) for v in line.split()])
    if not rows:
        return None
    # Rows may have different lengths if the last line is shorter — pad with 0
    maxlen = max(len(r) for r in rows)
    grid = np.zeros((len(rows), maxlen), dtype=np.float32)
    for i, r in enumerate(rows):
        grid[i, :len(r)] = r
    return grid


def list_snapshot_files(dirpath: Path, pattern: str) -> list[Path]:
    """Return sorted list of snapshot files matching *pattern* in *dirpath*."""
    return sorted(dirpath.glob(pattern))


def read_adios_snapshots(bp_path: Path, nx: int, ny: int, nz: int):
    """
    Read SEM snapshots from an ADIOS2 BP5 file.
    Returns (data_list, timestep_list) where each data is shape (nx, ny, nz),
    or (None, None) if adios2 is unavailable.
    """
    if not bp_path.exists():
        return None, None
    try:
        import adios2
        data_list, ts_list = [], []
        with adios2.Stream(str(bp_path), "r") as s:
            for _ in s.steps():
                step = s.current_step()
                try:
                    data = s.read("PressureField")
                except Exception:
                    continue
                # Reorder axes to (nx, ny, nz) if needed
                if data.shape == (nz, ny, nx):
                    data = data.transpose(2, 1, 0)
                elif data.shape == (ny, nz, nx):
                    data = data.transpose(2, 0, 1)
                data_list.append(data.astype(np.float32))
                ts_list.append(step)
        return (data_list, ts_list) if data_list else (None, None)
    except ImportError:
        print("  [INFO] adios2 Python module not available — SEM snapshots via .bp skipped")
        return None, None
    except Exception as e:
        print(f"  [WARN] Could not read {bp_path}: {e}")
        return None, None


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _clim(grid: np.ndarray, pct: float = 98.0) -> float:
    """Symmetric colour limit at the given percentile of absolute values."""
    v = float(np.percentile(np.abs(grid[np.isfinite(grid)]), pct))
    return v


def _global_vmax(grids: list[np.ndarray], pct: float = 98.0) -> float:
    """Global colour limit across all frames — avoids early zero-frames inflating vmax."""
    v = max(_clim(g, pct) for g in grids)
    return v if v > 0 else 1.0


def plot_xz_slice(grid: np.ndarray, title: str, outfile: Path,
                  lx: float = LX, lz_total: float = LZ,
                  src_z: float | None = None, rcv_z: float | None = None,
                  src_x: float | None = None, rcv_x: float | None = None,
                  vmax: float | None = None):
    """
    Plot a 2-D (nrows × ncols) XZ cross-section and save to *outfile*.

    The grid has rows = z direction (row 0 = z=0, surface) and
    cols = x direction (col 0 = x=0).  Z increases downward in the figure.
    """
    if vmax is None:
        vmax = _clim(grid)

    fig, ax = plt.subplots(figsize=(8, 6))
    # extent=[left, right, bottom, top] with origin="upper"
    # → row 0 shown at top (y=top=0), last row at bottom (y=bottom=lz_total)
    ext = [0, lx, lz_total, 0]
    im = ax.imshow(grid, cmap="seismic", vmin=-vmax, vmax=vmax,
                   aspect="auto", extent=ext, origin="upper")
    if src_z is not None and src_x is not None:
        ax.plot(src_x, src_z, "k^", ms=7, zorder=5, label=f"Source z={src_z:.0f} m")
    if rcv_z is not None and rcv_x is not None:
        ax.plot(rcv_x, rcv_z, "kv", ms=7, zorder=5, label=f"Receiver z={rcv_z:.0f} m")
    if src_z is not None or rcv_z is not None:
        ax.legend(fontsize=8)
    plt.colorbar(im, ax=ax, label="Pressure (a.u.)")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Z (m)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()


def plot_xy_slice(grid: np.ndarray, title: str, outfile: Path,
                  lx: float = LX, ly: float = LY,
                  vmax: float | None = None):
    """Plot a horizontal (ny × nx) XY slice and save to *outfile*."""
    if vmax is None:
        vmax = _clim(grid)
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(grid, cmap="seismic", vmin=-vmax, vmax=vmax,
                   aspect="auto", extent=[0, lx, ly, 0], origin="upper")
    plt.colorbar(im, ax=ax, label="Pressure (a.u.)")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()


# ---------------------------------------------------------------------------
# Per-method snapshot plotting
# ---------------------------------------------------------------------------

def plot_snapshots_sem(sim_dir: Path, snap_dir: Path, use_adios: bool):
    """
    SEM: prefer XZ slice from snapshots.bp (ADIOS2); fall back to XY .dat files.
    The XZ slice (at y=mid) shows the vertical wave propagation from source to receiver.
    """
    snap_dir.mkdir(parents=True, exist_ok=True)
    used_adios = False

    if use_adios:
        data_list, ts_list = read_adios_snapshots(sim_dir / "snapshots.bp", NX, NY, NZ)
        if data_list:
            used_adios = True
            xz_list = [data[:, NY // 2, :].T for data in data_list]  # (NZ, NX) each
            global_vmax = _global_vmax(xz_list)
            print(f"  SEM: {len(data_list)} ADIOS2 snapshots → XZ cross-sections (vmax={global_vmax:.4g})")
            for xz, ts in zip(xz_list, ts_list):
                t_phys = ts * DT
                title = (f"SEM — XZ cross-section @ y={LY/2:.0f} m  "
                         f"(t={t_phys:.4f} s, step {ts})")
                outfile = snap_dir / f"sem_xz_{ts:05d}.png"
                plot_xz_slice(xz, title, outfile,
                              lx=LX, lz_total=LZ,
                              src_z=SRCZ, src_x=SRCX,
                              rcv_z=RCVZ, rcv_x=RCVX,
                              vmax=global_vmax)
            return

    # Fallback: horizontal XY .dat slices at z ≈ LZ/2
    dat_files = list_snapshot_files(sim_dir, "slice_?????.dat")
    if not dat_files:
        print("  SEM: no snapshot files found")
        return
    grids_sem = [(f, read_dat_slice(f)) for f in dat_files]
    grids_sem = [(f, g) for f, g in grids_sem if g is not None]
    global_vmax = _global_vmax([g for _, g in grids_sem])
    print(f"  SEM: {len(grids_sem)} .dat snapshots → XY cross-sections @ z≈{LZ/2:.0f} m (vmax={global_vmax:.4g})")
    for fpath, grid in grids_sem:
        ts_str = fpath.stem.replace("slice_", "")
        ts = int(ts_str) if ts_str.isdigit() else 0
        t_phys = ts * DT
        title = (f"SEM — XY cross-section @ z≈{LZ/2:.0f} m  "
                 f"(t={t_phys:.4f} s, step {ts})")
        outfile = snap_dir / f"sem_xy_{ts:05d}.png"
        plot_xy_slice(grid, title, outfile, vmax=global_vmax)


def plot_snapshots_dg(sim_dir: Path, snap_dir: Path):
    """
    DG: horizontal XY slices saved as slice_dg_XXXXX.dat.
    Each file shape: (EY*N1D, EX*N1D) = (200, 200).
    The slice is taken at the mid-depth of the domain (z≈1000 m).
    """
    snap_dir.mkdir(parents=True, exist_ok=True)
    dat_files = list_snapshot_files(sim_dir, "slice_dg_?????.dat")
    if not dat_files:
        print("  DG: no snapshot files found")
        return
    grids_dg = [(f, read_dat_slice(f)) for f in dat_files]
    grids_dg = [(f, g) for f, g in grids_dg if g is not None]
    global_vmax = _global_vmax([g for _, g in grids_dg])
    print(f"  DG: {len(grids_dg)} snapshots → XY cross-sections @ z≈{LZ/2:.0f} m (vmax={global_vmax:.4g})")
    for fpath, grid in grids_dg:
        ts_str = fpath.stem.replace("slice_dg_", "")
        ts = int(ts_str) if ts_str.isdigit() else 0
        t_phys = ts * DT
        title = (f"DG — XY slice @ z≈{LZ/2:.0f} m  "
                 f"(t={t_phys:.4f} s, step {ts})")
        outfile = snap_dir / f"dg_xy_{ts:05d}.png"
        plot_xy_slice(grid, title, outfile, vmax=global_vmax)


def plot_snapshots_dgsem(sim_dir: Path, snap_dir: Path):
    """
    DG-SEM: combined XZ slice from slice_dgsem_xz_XXXXX.dat.
    The file has DG rows (z=0..1000 m) followed by SEM rows (z=1000..2000 m),
    with 151 columns (x nodes).  This gives a clear view of the wave crossing
    the DG-SEM interface.
    """
    snap_dir.mkdir(parents=True, exist_ok=True)
    dat_files = list_snapshot_files(sim_dir, "slice_dgsem_xz_?????.dat")
    if not dat_files:
        print("  DG-SEM: no snapshot files found")
        return

    # Total physical z for a combined grid
    # DG part: EZ_DG*N1D rows → LZ/2 = 1000 m
    # SEM part: NZ - NZ_DG_NODES rows → LZ/2 = 1000 m
    n_dg_rows = EZ_DG * N1D         # 100
    n_sem_rows = NZ - NZ_DG_NODES   # 76
    n_total_rows = n_dg_rows + n_sem_rows  # 176

    # Build a uniform physical z axis (combined, 176 rows → 2000 m)
    # DG rows: z=0..1000 with step = LZ/2 / n_dg_rows
    # SEM rows: z=1000..2000 with step = LZ/2 / n_sem_rows
    z_dg  = np.linspace(0,       LZ / 2, n_dg_rows,  endpoint=False)
    z_sem = np.linspace(LZ / 2,  LZ,     n_sem_rows + 1)[1:]
    z_combined = np.concatenate([z_dg, z_sem])  # (176,)

    grids_dgsem = [(f, read_dat_slice(f)) for f in dat_files]
    grids_dgsem = [(f, g) for f, g in grids_dgsem if g is not None]
    global_vmax = _global_vmax([
        g[:min(g.shape[0], n_total_rows), :min(g.shape[1], NX)]
        for _, g in grids_dgsem
    ])
    print(f"  DG-SEM: {len(grids_dgsem)} XZ snapshots → combined cross-sections (vmax={global_vmax:.4g})")
    for fpath, grid in grids_dgsem:
        ts_str = fpath.stem.replace("slice_dgsem_xz_", "")
        ts = int(ts_str) if ts_str.isdigit() else 0
        t_phys = ts * DT

        rows = min(grid.shape[0], n_total_rows)
        cols = min(grid.shape[1], NX)
        sub = grid[:rows, :cols]

        x_axis = np.linspace(0, LX, cols)

        fig, ax = plt.subplots(figsize=(8, 7))
        ax.pcolormesh(x_axis, z_combined[:rows], sub,
                      cmap="seismic", vmin=-global_vmax, vmax=global_vmax, shading="auto")
        ax.axhline(DGSEM_BOUNDARY_Z, color="k", lw=1.4, ls="--",
                   label=f"DG/SEM interface z={DGSEM_BOUNDARY_Z:.0f} m")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Z (m)")
        ax.set_ylim(LZ, 0)  # z increases downward
        ax.set_title(f"DG-SEM — XZ combined slice @ y={LY/2:.0f} m\n"
                     f"(t={t_phys:.4f} s, step {ts})")
        sm = plt.cm.ScalarMappable(cmap="seismic",
                                   norm=plt.Normalize(vmin=-global_vmax, vmax=global_vmax))
        plt.colorbar(sm, ax=ax, label="Pressure (a.u.)")
        plt.tight_layout()
        outfile = snap_dir / f"dgsem_xz_{ts:05d}.png"
        plt.savefig(outfile, dpi=150)
        plt.close()


# ---------------------------------------------------------------------------
# Receiver pressure plots
# ---------------------------------------------------------------------------

COLORS = {"sem": "#1f77b4", "dg": "#ff7f0e", "dgsem": "#2ca02c"}
LABELS = {"sem": "SEM", "dg": "DG", "dgsem": "DG-SEM"}


def plot_single_receiver(t: np.ndarray, p: np.ndarray, method: str, outfile: Path):
    """One receiver trace on its own figure."""
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t, p, color=COLORS[method], lw=1.0)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Pressure at receiver")
    ax.set_title(f"{LABELS[method]} — receiver pressure vs. time\n"
                 f"Source ({SRCX:.0f},{SRCY:.0f},{SRCZ:.0f}) m  "
                 f"→  Receiver ({RCVX:.0f},{RCVY:.0f},{RCVZ:.0f}) m")
    ax.grid(True, lw=0.4, alpha=0.5)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"  Saved: {outfile}")


def plot_overlay(traces: dict[str, tuple], methods: list[str], title: str, outfile: Path):
    """Overlay two or more receiver traces on one figure."""
    fig, ax = plt.subplots(figsize=(10, 4))
    for m in methods:
        if m not in traces or traces[m][0] is None:
            print(f"  [WARN] No data for {m} — skipping in overlay {outfile.name}")
            continue
        t, p = traces[m]
        ax.plot(t, p, color=COLORS[m], lw=1.0, label=LABELS[m])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Pressure at receiver (a.u.)")
    ax.set_title(f"{title}\n"
                 f"Source ({SRCX:.0f},{SRCY:.0f},{SRCZ:.0f}) m  "
                 f"→  Receiver ({RCVX:.0f},{RCVY:.0f},{RCVZ:.0f}) m")
    ax.legend()
    ax.grid(True, lw=0.4, alpha=0.5)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"  Saved: {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="comparison_results",
                   help="Root output directory (default: comparison_results)")
    p.add_argument("--skip-sim", action="store_true",
                   help="Skip simulations, only regenerate plots")
    p.add_argument("--no-snapshots", action="store_true",
                   help="Disable snapshot saving (faster, no cross-section plots)")
    p.add_argument("--np", type=int, default=1, dest="np_procs",
                   help="Number of MPI processes (default: 1)")
    p.add_argument("--snap-interval", type=int, default=DEFAULT_SNAP_INTERVAL, dest="snap_interval",
                   help=f"Steps between snapshots (default: {DEFAULT_SNAP_INTERVAL})")
    return p.parse_args()


def main():
    args = parse_args()
    use_snapshots = not args.no_snapshots

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    plots_dir = outdir / "plots"
    plots_dir.mkdir(exist_ok=True)

    # ------------------------------------------------------------------ binary
    env_bin = os.environ.get("FUNTIDES_BIN")
    if env_bin:
        binary = Path(env_bin)
    else:
        try:
            binary = find_executable(PROJECT_ROOT)
        except FileNotFoundError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)
    print(f"Binary: {binary}")

    methods = [
        ("sem",   outdir / "sem"),
        ("dg",    outdir / "dg"),
        ("dg-sem", outdir / "dgsem"),
    ]

    # ------------------------------------------------------------------ simulations
    if not args.skip_sim:
        print("\n" + "="*60)
        print(" SIMULATION PARAMETERS")
        print("="*60)
        print(f"  Domain        : {LX:.0f} × {LY:.0f} × {LZ:.0f} m³")
        print(f"  Elements      : {EX} × {EY} × {EZ}  (h = {LX/EX:.1f} m)")
        print(f"  Order         : {ORDER}  (N1D = {N1D})")
        print(f"  dt / timemax  : {DT} s / {TIMEMAX} s  ({N_STEPS} steps)")
        print(f"  f0 / tpeak    : {F0} Hz / {TPEAK:.4f} s")
        print(f"  Source        : ({SRCX},{SRCY},{SRCZ}) m")
        print(f"  Receiver      : ({RCVX},{RCVY},{RCVZ}) m")
        print(f"  Sponge size   : {BOUNDARIES_SIZE} m  (taper delta = {TAPER_DELTA})")
        print(f"  DG-SEM iface  : z = {DGSEM_BOUNDARY_Z} m")
        if use_snapshots:
            n_snaps = N_STEPS // args.snap_interval
            print(f"  Snapshots     : every {args.snap_interval} steps ({n_snaps} total)")
        else:
            print(f"  Snapshots     : disabled")
        print("="*60)

        for method, sim_dir in methods:
            rc = run_simulation(binary, method, sim_dir,
                                args.snap_interval, use_snapshots, args.np_procs)
            if rc != 0:
                print(f"\n[ERROR] {method.upper()} simulation failed (return code {rc})",
                      file=sys.stderr)

    # ------------------------------------------------------------------ receiver plots
    print("\n" + "="*60)
    print(" GENERATING RECEIVER TRACE PLOTS")
    print("="*60)

    traces: dict[str, tuple] = {}
    key_map = {"sem": "sem", "dg": "dg", "dg-sem": "dgsem"}

    for method, sim_dir in methods:
        key = key_map[method]
        t, p = read_receiver_trace(sim_dir)
        traces[key] = (t, p)

    for method, sim_dir in methods:
        key = key_map[method]
        t, p = traces[key]
        if t is None:
            print(f"  [WARN] No receiver data for {method} — skipping")
            continue
        plot_single_receiver(t, p, key, plots_dir / f"receiver_{key}.png")

    # DG vs SEM overlay
    plot_overlay(traces, ["dg", "sem"],
                 "Receiver pressure — DG vs. SEM",
                 plots_dir / "compare_dg_sem.png")

    # DG-SEM vs SEM overlay
    plot_overlay(traces, ["dgsem", "sem"],
                 "Receiver pressure — DG-SEM vs. SEM",
                 plots_dir / "compare_dgsem_sem.png")

    # ------------------------------------------------------------------ snapshot plots
    if use_snapshots:
        print("\n" + "="*60)
        print(" GENERATING SNAPSHOT CROSS-SECTION PLOTS")
        print("="*60)

        print("\n-- SEM snapshots")
        plot_snapshots_sem(outdir / "sem",
                           plots_dir / "snapshots_sem",
                           use_adios=True)

        print("\n-- DG snapshots")
        plot_snapshots_dg(outdir / "dg",
                          plots_dir / "snapshots_dg")

        print("\n-- DG-SEM snapshots")
        plot_snapshots_dgsem(outdir / "dgsem",
                             plots_dir / "snapshots_dgsem")

    # ------------------------------------------------------------------ summary
    print("\n" + "="*60)
    print(" DONE")
    print("="*60)
    print(f"All plots saved in: {plots_dir}")
    print("\nGenerated plots:")
    for f in sorted(plots_dir.rglob("*.png")):
        print(f"  {f.relative_to(plots_dir)}")


if __name__ == "__main__":
    main()
