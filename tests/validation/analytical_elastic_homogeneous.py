#!/usr/bin/env python3
"""
Closed-form analytical solution of the 2nd-order elastodynamic wave equation for
a homogeneous infinite isotropic medium with a point force, matching EXACTLY the
PDE and source convention solved by FUnTiDES's elastic SEM solver (same Ricker
formula/order, same force direction, no unit conversion, no derivative shift).

FUnTiDES semi-discrete system (src/solver/fe/impl/common/include/sem_solver_impl.h):
    M u'' + K u = F(t),
    M_I     = int rho * N_I^2 dV                     (computeGlobalMassMatrix, model_factor = rho)
    K       = int grad(N) : C : grad(N) dV           (elastic stiffness, iso)
    F_I,c(t) = - f(t) * N_I(xs)   for c = x, y, z    (applyRHSTerm uses '-=';
                                                      InitSource sets rhs_term_{x,y,z} = +f(t))

For a homogeneous medium this is the strong form (i = x, y, z):
    rho d2u_i/dt2 - div(sigma)_i = f(t) delta(x - xs)

i.e. a point force with time history X(t) = f(t) * (1, 1, 1) and f() the exact
same Ricker wavelet FUnTiDES injects, so no
derivative-order shift or amplitude factor is needed to compare against a
FUnTiDES trace. (The '-=' in applyRHSTerm is cancelled by the residual sign in
the time update -- same as the validated acoustic analytic, which also carries
no minus on f.)

Using the full-space elastodynamic Green's function (Aki & Richards, 2002,
eq. 4.29) for a point force X_p(t) at xs, the displacement at x is:

    gamma = (x - xs) / r,  r = |x - xs|,  alpha = Vp,  beta = Vs
    S = gamma_x + gamma_y + gamma_z          (force along (1, 1, 1))
    X_p(t) = f(t)   for every p

    u_n(x, t) =   (1 / (4 pi rho r^3)) (3 gamma_n S - 1) * J(t)              near field
                + (1 / (4 pi rho alpha^2 r)) gamma_n S * f(t - r/alpha)     far field P
                - (1 / (4 pi rho beta^2 r)) (gamma_n S - 1) * f(t - r/beta) far field S

    J(t) = int_{r/alpha}^{r/beta} tau * f(t - tau) dtau   (numeric quadrature)
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


def analytical_displacement(t, src, rcv, f0, order, tpeak, vp, rho, vs, terms="all"):
    """Full-space elastic displacement u(x, t) for a (1, 1, 1) point force X(t) = f(t).

    Returns an array of shape (len(t), 3) holding (ux, uy, uz).
    `terms` selects which contributions to keep: "all", "near", "far", "farp", "fars"
    (diagnostic -- to isolate a near-field vs far-field disagreement).
    """
    src = np.asarray(src, dtype=np.float64)
    rcv = np.asarray(rcv, dtype=np.float64)
    d = rcv - src
    r = float(np.linalg.norm(d))
    if r == 0.0:
        raise ValueError("source and receiver coincide: r=0 is a singularity of the Green's function")

    gamma = d / r  # direction cosines
    s = float(gamma.sum())  # gamma . (1, 1, 1)
    t = np.asarray(t, dtype=np.float64)

    # --- far-field P and S ---
    # The FUnTiDES point-force time history is +f(t): the '-=' in applyRHSTerm is
    # cancelled by the residual sign in the update, exactly as for the (validated)
    # acoustic analytic which also carries no minus on f.
    fp = evaluate_ricker(t - r / vp, f0, order, tpeak)
    fs = evaluate_ricker(t - r / vs, f0, order, tpeak)

    far_p = (1.0 / (4.0 * np.pi * rho * vp * vp * r)) * (gamma * s)[None, :] * fp[:, None]
    far_s = -(1.0 / (4.0 * np.pi * rho * vs * vs * r)) * (gamma * s - 1.0)[None, :] * fs[:, None]

    # --- near field: J(t) = int_{r/vp}^{r/vs} tau * f(t - tau) dtau ---
    tau_lo, tau_hi = r / vp, r / vs
    dt_samp = float(t[1] - t[0]) if t.size > 1 else (tau_hi - tau_lo)
    nsub = int(max(400, 20.0 * (tau_hi - tau_lo) / max(dt_samp, 1e-12)))
    tau = np.linspace(tau_lo, tau_hi, nsub)
    # integrand[i, j] = tau_j * f(t_i - tau_j)
    integrand = tau[None, :] * evaluate_ricker(t[:, None] - tau[None, :], f0, order, tpeak)
    trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
    j_t = trapz(integrand, tau, axis=1)  # shape (len(t),)

    near = (1.0 / (4.0 * np.pi * rho * r**3)) * (3.0 * gamma * s - 1.0)[None, :] * j_t[:, None]

    keep_near = terms in ("all", "near")
    keep_farp = terms in ("all", "far", "farp")
    keep_fars = terms in ("all", "far", "fars")
    return (near if keep_near else 0.0) + (far_p if keep_farp else 0.0) + (far_s if keep_fars else 0.0)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--f0", type=float, required=True, help="central frequency of the Ricker source (Hz)")
    p.add_argument("--ricker-order", type=int, required=True, choices=[0, 1, 2, 3, 4], dest="ricker_order", help="Ricker wavelet derivative order (matches --ricker-order passed to funtides)")
    p.add_argument("--tpeak", type=float, required=True, help="peak time of the Ricker wavelet (matches --tpeak)")
    p.add_argument("--vp", type=float, required=True, help="P-wave velocity of the homogeneous medium")
    p.add_argument("--vs", type=float, required=True, help="S-wave velocity of the homogeneous medium")
    p.add_argument("--rho", type=float, required=True, help="density of the homogeneous medium")
    p.add_argument("--src", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"), help="source coordinates")
    p.add_argument("--rcv", type=float, nargs=3, required=True, metavar=("X", "Y", "Z"), help="receiver coordinates")
    p.add_argument("--dt", type=float, required=True, help="time step (matches --dt)")
    p.add_argument("--nsamples", type=int, required=True, help="number of time samples (matches solver's num_samples_)")
    p.add_argument("--t0", type=float, default=0.0, help="start time of the first sample (default 0)")
    p.add_argument("--terms", choices=["all", "near", "far", "farp", "fars"], default="all", help="keep only these Green's-function terms (diagnostic)")
    p.add_argument("-o", "--output", default="analytical_elastic_trace.txt", help="output file (# time ux uy uz)")
    args = p.parse_args()

    t = args.t0 + args.dt * np.arange(args.nsamples)
    u = analytical_displacement(t, args.src, args.rcv, args.f0, args.ricker_order, args.tpeak, args.vp, args.rho, args.vs, args.terms)

    with open(args.output, "w") as fout:
        fout.write("# time ux uy uz\n")
        for ti, ui in zip(t, u):
            fout.write(f"{ti} {ui[0]} {ui[1]} {ui[2]}\n")

    r = float(np.linalg.norm(np.asarray(args.rcv) - np.asarray(args.src)))
    print(f"wrote {args.nsamples} samples to {args.output}")
    print(f"source-receiver distance r={r:.6g}")
    print(f"P arrival r/Vp={r / args.vp:.6g}s, S arrival r/Vs={r / args.vs:.6g}s")


if __name__ == "__main__":
    main()
