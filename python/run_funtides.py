#!/usr/bin/env python3
"""
run_funtides.py — Wrapper script to launch funtides-sem from a YAML config file.

Usage:
    python run_funtides.py config.yaml [extra CLI args passed to the binary]
    python run_funtides.py --help

YAML keys map 1-to-1 to the CLI flags of funtides-sem (dashes and underscores
are treated identically). CLI arguments provided after the config file override
the corresponding YAML values.

The reserved key ``mpi-np`` (integer) is not passed to the binary — it controls
whether the launch is prefixed with ``mpirun -np <mpi-np>``.

"""

import os
import subprocess
import sys
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# Executable auto-detection
# ---------------------------------------------------------------------------

BINARY_NAME = "funtides-sem"

# The script lives in <project_root>/python/, so go one level up.
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def find_executable(root: Path) -> Path:
    """Search for funtides-sem in build*/bin/ directories under *root*.

    If several copies are found, returns the most recently modified one.

    Args:
        root: Project root directory.

    Returns:
        Path to the funtides-sem binary.

    Raises:
        FileNotFoundError: If no binary is found.
    """
    candidates = sorted(
        root.glob(f"build*/bin/{BINARY_NAME}"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No '{BINARY_NAME}' found in any build*/bin/ directory under {root}.\n"
            "Build the project first (cmake + make)."
        )
    if len(candidates) > 1:
        print(
            f"[run_funtides] Multiple binaries found, using the most recent:\n"
            f"  {candidates[0]}\n"
            f"  (others: {', '.join(str(c) for c in candidates[1:])})\n"
        )
    return candidates[0]


# ---------------------------------------------------------------------------
# MPI launch prefix
# ---------------------------------------------------------------------------

# Keys consumed by this script — never forwarded to the binary.
_SCRIPT_KEYS = {"mpi-np", "mpi-launcher", "launcher"}

# Flag used by each supported launcher to set the number of processes.
_NP_FLAG: dict[str, str] = {
    "mpirun":  "-np",
    "mpiexec": "-n",
    "srun":    "-n",
    "jsrun":   "-n",
}


def _detect_scheduler() -> str | None:
    """Return the scheduler name if running inside a known job allocation."""
    if "SLURM_JOB_ID" in os.environ:
        return "slurm"
    if "LSB_JOBID" in os.environ:
        return "lsf"
    return None


def build_mpi_prefix(config: dict) -> list[str]:
    """Return the MPI launcher prefix to prepend to the command.

    Resolution order for the launcher:
    1. ``FUNTIDES_MPI_LAUNCHER`` environment variable
    2. ``mpi-launcher`` YAML key
    3. Auto-detected from the scheduler (SLURM → ``srun``, LSF → ``mpirun``)
    4. ``mpirun`` (default for direct launches)

    Resolution order for the number of processes:
    1. ``FUNTIDES_MPI_NP`` environment variable
    2. ``mpi-np`` YAML key (default: 1)

    When running under SLURM with ``srun`` and no explicit ``mpi-np``, the
    number of tasks is omitted so that srun inherits the allocation
    (``#SBATCH --ntasks``).

    Args:
        config: Parsed YAML configuration dictionary.

    Returns:
        Launcher prefix as a list of strings, or ``[]`` for a direct launch.
    """
    # mpi-np / mpi-launcher may live under a "launcher:" section or at root
    launcher_section = config.get("launcher") or {}

    # Resolve number of processes
    env_np = os.environ.get("FUNTIDES_MPI_NP")
    np = int(env_np) if env_np is not None else int(
        launcher_section.get("mpi-np", config.get("mpi-np", 1))
    )

    # Resolve launcher
    launcher = (
        os.environ.get("FUNTIDES_MPI_LAUNCHER")
        or launcher_section.get("mpi-launcher")
        or config.get("mpi-launcher")
    )
    if launcher is None:
        scheduler = _detect_scheduler()
        if scheduler == "slurm":
            launcher = "srun"
        elif scheduler == "lsf":
            launcher = "mpirun"
        elif np > 1:
            launcher = "mpirun"
        else:
            return []  # single-process direct launch

    launcher = str(launcher)
    np_flag = _NP_FLAG.get(launcher, "-np")

    # Under srun inside a SLURM allocation with no explicit np,
    # omit -n and let srun inherit --ntasks from the batch header.
    if launcher == "srun" and np <= 1 and _detect_scheduler() == "slurm":
        return ["srun"]

    # Warn if mpi-np conflicts with the scheduler allocation size.
    if np > 1:
        scheduler = _detect_scheduler()
        allocated_str = None
        scheduler_label = ""
        if scheduler == "slurm":
            allocated_str = os.environ.get("SLURM_NTASKS")
            scheduler_label = "SLURM (SLURM_NTASKS)"
        elif scheduler == "lsf":
            allocated_str = os.environ.get("LSB_DJOB_NUMPROC")
            scheduler_label = "LSF (LSB_DJOB_NUMPROC)"

        if allocated_str is not None and int(allocated_str) != np:
            allocated = int(allocated_str)
            print(
                f"[run_funtides] WARNING: mpi-np={np} differs from the "
                f"{scheduler_label} allocation ({allocated}).",
                file=sys.stderr,
            )
            if np > allocated:
                print(
                    f"[run_funtides] ERROR: requesting more tasks ({np}) than "
                    f"allocated ({allocated}). The scheduler will reject this.",
                    file=sys.stderr,
                )
                sys.exit(1)
            else:
                print(
                    f"[run_funtides] WARNING: only {np}/{allocated} allocated "
                    "tasks will be used.",
                    file=sys.stderr,
                )

    if np > 1:
        return [launcher, np_flag, str(np)]
    return [launcher]


def _normalize_key(key: str) -> str:
    """Convert underscores to dashes so both styles are accepted in YAML."""
    return key.replace("_", "-")


def yaml_to_cli_args(config: dict) -> list[str]:
    """Convert a YAML config dict to a list of CLI arguments.

    The config may be flat or organized into arbitrary sections (nested dicts).
    Section names are ignored — only leaf values are converted.

    Mapping rules for leaf values:
    - bool True  → ``--key``
    - bool False → omitted (leaves binary default unchanged)
    - list       → ``--key val1,val2,val3``
    - empty list → omitted
    - dict       → recurse (section grouping, name ignored)
    - other      → ``--key value``

    Args:
        config: Parsed YAML configuration dictionary.

    Returns:
        List of CLI argument strings.
    """
    args: list[str] = []
    for raw_key, value in config.items():
        key = _normalize_key(raw_key)
        if key in _SCRIPT_KEYS:
            continue
        if value is None:
            continue  # commented-out or empty YAML section
        if isinstance(value, dict):
            args.extend(yaml_to_cli_args(value))  # section — recurse, name ignored
        elif isinstance(value, bool):
            if value:
                args.append(f"--{key}")
        elif isinstance(value, list):
            if value:
                args.extend([f"--{key}", ",".join(str(v) for v in value)])
        else:
            args.extend([f"--{key}", str(value)])
    return args


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------

HELP_TEXT = f"""\
Usage:
  python run_funtides.py <config.yaml> [extra args ...]
  python run_funtides.py --help

Arguments:
  config.yaml   Path to a YAML configuration file. Keys match the CLI flags
                of {BINARY_NAME} (dashes and underscores interchangeable).
  extra args    Any additional arguments are appended as-is to the command
                line, overriding the YAML values where applicable.

Binary selection:
  If the environment variable FUNTIDES_BIN is set, it is used directly.
  Otherwise the script searches for '{BINARY_NAME}' in build*/bin/ under the
  project root, and picks the most recently modified binary if several exist.

MPI:
  mpi-np: N          Number of MPI processes (default: 1 = direct launch).
  mpi-launcher: X    Launcher to use: mpirun | mpiexec | srun | jsrun.
                     Auto-detected when inside a SLURM or LSF allocation
                     (SLURM → srun, LSF → mpirun).  Under srun, omitting
                     mpi-np lets the job inherit --ntasks from the batch header.
  Environment overrides: FUNTIDES_MPI_NP, FUNTIDES_MPI_LAUNCHER.
"""


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    args = sys.argv[1:]

    if not args or args[0] in ("-h", "--help"):
        print(HELP_TEXT)
        return 0

    config_path = Path(args[0])
    extra_args = args[1:]

    if not config_path.exists():
        print(f"[run_funtides] Error: config file not found: {config_path}", file=sys.stderr)
        return 1

    try:
        with config_path.open() as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as exc:
        print(f"[run_funtides] Error parsing {config_path}: {exc}", file=sys.stderr)
        return 1

    # Locate executable — honour FUNTIDES_BIN if set, otherwise auto-detect
    env_bin = os.environ.get("FUNTIDES_BIN")
    if env_bin:
        binary = Path(env_bin)
        if not binary.exists():
            print(f"[run_funtides] FUNTIDES_BIN points to a non-existent path: {binary}", file=sys.stderr)
            return 1
        print(f"[run_funtides] Binary (FUNTIDES_BIN): {binary}\n")
    else:
        try:
            binary = find_executable(PROJECT_ROOT)
        except FileNotFoundError as exc:
            print(f"[run_funtides] {exc}", file=sys.stderr)
            return 1

    mpi_prefix = build_mpi_prefix(config)
    cli_args = yaml_to_cli_args(config)
    cmd = mpi_prefix + [str(binary)] + cli_args + extra_args

    print(f"[run_funtides] Binary : {binary}")
    print(f"[run_funtides] Command: {' '.join(cmd)}\n")

    result = subprocess.run(cmd)
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
