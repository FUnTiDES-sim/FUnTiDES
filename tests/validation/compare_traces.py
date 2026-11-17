#!/usr/bin/env python3
"""
Compare a FUnTiDES receiver trace against a reference trace (typically written
by analytical_acoustic_homogeneous.py or analytical_elastic_homogeneous.py).

Both files carry a '#'-prefixed header line naming the columns, e.g.
    # time pressure        or        # time ux uy uz

The reference is linearly resampled onto the numerical trace's time grid, over
the overlapping time range only, so a different dt / sample count between the
two is handled. No per-trace normalization is applied: the analytical scripts
match the PDE and source convention of FUnTiDES, so amplitudes must agree too.

Reports the relative L2 error of each component and the total one

    sqrt( sum_c int (num_c - ref_c)^2 dt ) / sqrt( sum_c int ref_c^2 dt )

Usage:
    python3 compare_traces.py receiver_trace.txt analytical_trace.txt
    python3 compare_traces.py num.txt ref.txt --components ux --plot
"""
import argparse

import numpy as np

trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)


def read_header_names(path):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                return line.lstrip("#").split() if line.startswith("#") else None
    return None


def load_trace(path, components=None):
    """Return (time, dict{component: values}). Columns resolved by header name, all of them by default."""
    names = read_header_names(path)
    if names is None:
        raise SystemExit(f"{path}: no '#'-prefixed header line to resolve column names")
    if "time" not in names:
        raise SystemExit(f"{path}: header has no 'time' column: {names}")
    if components is None:
        components = [n for n in names if n != "time"]
    data = np.loadtxt(path, ndmin=2)
    out = {}
    for c in components:
        if c not in names:
            raise SystemExit(f"{path}: column '{c}' not found in header: {names}")
        out[c] = data[:, names.index(c)]
    return data[:, names.index("time")], out


def relative_l2_errors(t_num, num, t_ref, ref, lag=0):
    """
    Relative L2 error of each component of num against ref, and the total one.

    ref is resampled onto t_num, shifted by `lag` samples. Returns
    (t, ref_resampled, {component: error}, total_error), t being the part of
    t_num inside the time range of ref.
    """
    mask = (t_num >= max(t_num.min(), t_ref.min())) & (t_num <= min(t_num.max(), t_ref.max()))
    t = t_num[mask]
    shift = lag * float(np.median(np.diff(t))) if lag else 0.0
    ref_r = {c: np.interp(t + shift, t_ref, ref[c]) for c in num}
    errors = {}
    diff_tot = ref_tot = 0.0
    for c in num:
        diff_sq = trapz((num[c][mask] - ref_r[c]) ** 2, t)
        ref_sq = trapz(ref_r[c] ** 2, t)
        errors[c] = np.sqrt(diff_sq / ref_sq) if ref_sq > 0 else float("nan")
        diff_tot += diff_sq
        ref_tot += ref_sq
    total = np.sqrt(diff_tot / ref_tot) if ref_tot > 0 else float("nan")
    return t, ref_r, errors, total


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("numerical", help="FUnTiDES trace (receiver_trace.txt)")
    p.add_argument("reference", help="analytical trace (analytical_trace.txt)")
    p.add_argument("--components", default=None, help="comma-separated component names (default: every column but time)")
    p.add_argument("--plot", action="store_true", help="show an overlay + difference plot")
    p.add_argument("--kinematics-check", action="store_true", help="also test num against d/dt(ref) and int(ref) -- catches a velocity/acceleration receiver")
    p.add_argument("--align", action="store_true", help="find & apply the integer-sample lag minimizing L2 (detects an off-by-one time convention)")
    p.add_argument("--align-window", type=int, default=50, help="max |lag| in samples searched by --align (default 50)")
    args = p.parse_args()

    comps = [c.strip() for c in args.components.split(",") if c.strip()] if args.components else None
    tn, num = load_trace(args.numerical, comps)
    tr, ref = load_trace(args.reference, list(num))

    lag = 0
    if args.align:
        w = abs(args.align_window)
        lag = min(range(-w, w + 1), key=lambda s: relative_l2_errors(tn, num, tr, ref, s)[3])
        dt = float(np.median(np.diff(tn)))
        print(f"best integer-sample lag: {lag}  ({lag * dt:+.4g} s)  -- ref shifted by this to match num")
        if abs(lag) == w:
            print(f"  WARNING: lag hit the search window edge (+/-{w}); widen with --align-window or check the wave velocities")

    t, ref_r, errors, total = relative_l2_errors(tn, num, tr, ref, lag)
    mask = (tn >= t.min()) & (tn <= t.max())
    num_r = {c: num[c][mask] for c in num}

    print(f"samples compared : {len(t)}  (time range [{t.min():.6g}, {t.max():.6g}])")
    for c, e in errors.items():
        print(f"{c:<10} L2 relative : {e:.6e}  ({100 * e:.4f} %)")
    print(f"{'-' * 60}")
    print(f"total L2 relative     : {total:.6e}  ({100 * total:.4f} %)")

    if args.kinematics_check:
        # Is FUnTiDES really outputting displacement, or velocity / acceleration?
        # Compare the numeric trace against the reference, its time-derivative and
        # its time-integral (each with both signs). Lowest number wins.
        def total_rel(transform):
            n2 = r2 = 0.0
            for c in num_r:
                g = transform(ref_r[c])
                n2 += trapz((num_r[c] - g) ** 2, t)
                r2 += trapz(g**2, t)
            return np.sqrt(n2) / np.sqrt(r2) if r2 > 0 else float("nan")

        cumtrapz = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])
        variants = {
            "+ref": lambda y: y,
            "-ref": lambda y: -y,
            "+d/dt ref": lambda y: np.gradient(y, t),
            "-d/dt ref": lambda y: -np.gradient(y, t),
            "+int ref": cumtrapz,
            "-int ref": lambda y: -cumtrapz(y),
        }
        print("kinematics check (total L2 rel vs each transform of ref):")
        for name, fn in variants.items():
            print(f"  num vs {name:<10} : {100 * total_rel(fn):8.4f} %")

    if args.plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(len(num_r), 2, figsize=(12, 2.6 * len(num_r)), squeeze=False)
        for row, c in enumerate(num_r):
            ax_o, ax_d = axes[row]
            ax_o.plot(t, num_r[c], lw=1, label="FUnTiDES")
            ax_o.plot(t, ref_r[c], lw=1, ls="--", label="analytic")
            ax_o.set_ylabel(c)
            ax_o.legend(loc="upper right", fontsize=8)
            ax_d.plot(t, num_r[c] - ref_r[c], color="crimson", lw=1)
            ax_d.set_ylabel(f"{c} diff")
        axes[-1][0].set_xlabel("time")
        axes[-1][1].set_xlabel("time")
        fig.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
