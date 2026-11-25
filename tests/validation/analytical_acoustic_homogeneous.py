#!/usr/bin/env python3
"""
Closed-form analytical solution of the 2nd-order acoustic wave equation for a
homogeneous infinite medium with a point source, matching EXACTLY the PDE
solved by FUnTiDES's acoustic SEM solver (same source convention, same Ricker
formula/order, no unit conversion).

FUnTiDES semi-discrete system (see src/solver/fe/impl/common/include/sem_solver_impl.h):
    M p'' + K p = F(t),
    M_I = int N_I^2 / (rho * Vp^2) dV   (computeGlobalMassMatrix, model_factor)
    K_I = int (1/rho) grad(N_I).grad(N_J) dV
    F_I(t) = f(t) * N_I(xs)             (applyRHSTerm, raw Ricker, no scaling)

For a homogeneous medium this is the strong form:
    (1/(rho c^2)) d2p/dt2 - div((1/rho) grad p) = f(t) delta(x - xs)
<=> d2p/dt2 - c^2 lap(p) = rho c^2 f(t) delta(x - xs)

Using the retarded Green's function of (d2/dt2 - c^2 lap) [G = delta(t-r/c)/(4 pi c^2 r)],
the rho*c^2 cancels the Green's function's 1/c^2, leaving the closed form:

    p(x, t) = rho * f(t - r/c) / (4 pi r)

with r = |x - xs| and f() the exact same Ricker wavelet FUnTiDES injects.
"""
import argparse

import numpy as np


def evaluate_ricker(time_n, f0, order, tpeak):
    """Bit-for-bit port of SolverUtils::evaluateRicker (src/acquisition/include/source_time_function.h)."""
    time_n = np.asarray(time_n, dtype=np.float64)
    pulse = np.zeros_like(time_n)
    mask = (time_n > -0.9 * tpeak) & (time_n < 2.9 * tpeak)

    lam = (f0 * np.pi) ** 2
    t1 = time_n[mask] - tpeak

    if order == 4:
        val = (4.0 * lam * lam * (3.0 - 12.0 * lam * t1**2 + 4.0 * lam * lam * t1**4) * np.exp(-lam * t1**2))
    elif order == 3:
        val = 4.0 * lam * lam * t1 * (3.0 - 2.0 * lam * t1**2) * np.exp(-lam * t1**2)
    elif order == 2:
        val = 2.0 * lam * (2.0 * lam * t1**2 - 1.0) * np.exp(-lam * t1**2)
    elif order == 1:
        val = -2.0 * lam * t1 * np.exp(-lam * t1**2)
    elif order == 0:
        val = -t1 * np.exp(-2.0 * lam * t1**2)
    else:
        raise ValueError("ricker order must be 0..4")

    pulse[mask] = val
    return pulse


def analytical_pressure(t, src, rcv, f0, order, tpeak, vp, rho):
    """p(x, t) = rho * f(t - r/c) / (4 pi r) for a homogeneous infinite medium."""
    r = float(np.linalg.norm(np.asarray(rcv, dtype=np.float64) - np.asarray(src, dtype=np.float64)))
    if r == 0.0:
        raise ValueError("source and receiver coincide: r=0 is a singularity of the Green's function")
    retarded_t = t - r / vp
    return rho * evaluate_ricker(retarded_t, f0, order, tpeak) / (4.0 * np.pi * r)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--f0", type=float, required=True, help="central frequency of the Ricker source (Hz)")
    p.add_argument("--ricker-order", type=int, required=True, choices=[0, 1, 2, 3, 4], dest="ricker_order", help="Ricker wavelet derivative order (matches --ricker-order passed to funtides)")
    p.add_argument("--tpeak", type=float, required=True, help="peak time of the Ricker wavelet (matches --tpeak)")
    p.add_argument("--vp", type=float, required=True, help="P-wave velocity of the homogeneous medium")
    p.add_argument("--rho", type=float, required=True, help="density of the homogeneous medium")
    p.add_argument("--src", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"), help="source coordinates")
    p.add_argument("--rcv", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"), help="receiver coordinates")
    p.add_argument("--dt", type=float, required=True, help="time step (matches --dt)")
    p.add_argument("--nsamples", type=int, required=True, help="number of time samples (matches solver's num_samples_)")
    p.add_argument("--t0", type=float, default=0.0, help="start time of the first sample (default 0)")
    p.add_argument("-o", "--output", default="analytical_trace.txt", help="output file (same format as receiver_trace.txt)")
    args = p.parse_args()

    t = args.t0 + args.dt * np.arange(args.nsamples)
    p_t = analytical_pressure(t, args.src, args.rcv, args.f0, args.ricker_order, args.tpeak, args.vp, args.rho)

    with open(args.output, "w") as fout:
        fout.write("# time pressure_at_receiver\n")
        for ti, pi in zip(t, p_t):
            fout.write(f"{ti} {pi}\n")

    print(f"wrote {args.nsamples} samples to {args.output}")
    r = float(np.linalg.norm(np.asarray(args.rcv) - np.asarray(args.src)))
    print(f"source-receiver distance r={r:.6g}, arrival time r/c={r / args.vp:.6g}s")


if __name__ == "__main__":
    main()
