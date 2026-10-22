#!/usr/bin/env python3
"""
Validate the SEM, DG, DG-SEM and p-adaptive DG solvers against closed-form
homogeneous solutions.

For each selected case (see cases.py), three steps:
  1. run the solver through pyfuntides, write <output-dir>/<case>/receiver_trace.txt
  2. evaluate the analytical solution on the same time grid, write analytical_trace.txt
  3. relative L2 error between both traces (compare_traces.py), checked against
     the case tolerance

Both trace files can be compared again by hand with compare_traces.py
(--plot, --align, --kinematics-check).

Usage:
    python3 run_validation.py                          # every case
    python3 run_validation.py --method dgsem dgpadaptive --order 3   # p-adaptive: --order is pMax
    python3 run_validation.py --physics acoustic --mesh cartesian --order 2
    python3 run_validation.py --on-nodes yes --list
"""
import argparse
import gc
import os
import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import cases
from analytical_acoustic_homogeneous import analytical_pressure
from analytical_elastic_homogeneous import analytical_displacement
from compare_traces import relative_l2_errors


@dataclass
class CaseResult:
    case: cases.ValidationCase
    errors: dict
    total_error: float
    elapsed: float

    @property
    def passed(self):
        return bool(self.total_error < self.case.tolerance)


def analytical_trace(case, t, vp, vs, rho):
    """Analytical solution at the receiver, one column per entry of case.components."""
    common = (t, case.source, case.receiver, cases.F0, cases.RICKER_ORDER, cases.TPEAK, vp, rho)
    if case.is_elastic:
        return analytical_displacement(*common, vs)
    return analytical_pressure(*common)[:, None]


def write_trace(path, t, trace, components):
    np.savetxt(path, np.column_stack([t, trace]), header="time " + " ".join(components))


def validate_case(case, output_dir):
    """Run, evaluate and compare one case; the traces are written under output_dir/<case name>."""
    import simulation  # imports kokkos and pyfuntides, deferred so that --list works without them

    case_dir = Path(output_dir) / case.name
    case_dir.mkdir(parents=True, exist_ok=True)

    start = time.perf_counter()
    t, trace, (vp, vs, rho) = simulation.run(case)
    elapsed = time.perf_counter() - start
    write_trace(case_dir / "receiver_trace.txt", t, trace, case.components)

    reference = analytical_trace(case, t, vp, vs, rho)
    write_trace(case_dir / "analytical_trace.txt", t, reference, case.components)

    num = {c: trace[:, i] for i, c in enumerate(case.components)}
    ref = {c: reference[:, i] for i, c in enumerate(case.components)}
    _, _, errors, total = relative_l2_errors(t, num, t, ref)
    return CaseResult(case, errors, total, elapsed)


def select_cases(args):
    on_nodes = {"yes": True, "no": False}
    return [
        c
        for c in cases.all_cases()
        if c.method in args.method
        and c.physics in args.physics
        and c.mesh in args.mesh
        and c.order in args.order
        and c.on_nodes in {on_nodes[v] for v in args.on_nodes}
    ]


def print_summary(results, failures):
    print(f"\n{'=' * 78}")
    print(f"{'case':<42} {'L2 rel':>10} {'tol':>8} {'time (s)':>9}  status")
    print(f"{'-' * 78}")
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        print(f"{r.case.name:<42} {r.total_error:>10.3e} {r.case.tolerance:>8.2g} {r.elapsed:>9.1f}  {status}")
    for case in failures:
        print(f"{case.name:<42} {'':>10} {case.tolerance:>8.2g} {'':>9}  ERROR")
    print(f"{'=' * 78}")
    n_ok = sum(r.passed for r in results)
    print(f"{n_ok}/{len(results) + len(failures)} cases passed")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--method", nargs="+", choices=cases.METHODS, default=list(cases.METHODS))
    p.add_argument("--physics", nargs="+", choices=cases.PHYSICS, default=list(cases.PHYSICS))
    p.add_argument("--mesh", nargs="+", choices=cases.MESHES, default=list(cases.MESHES))
    p.add_argument("--order", nargs="+", type=int, choices=cases.ORDERS, default=list(cases.ORDERS))
    p.add_argument("--on-nodes", nargs="+", choices=["yes", "no"], default=["no", "yes"], dest="on_nodes",
                   help="model given on nodes (yes), on elements (no), or both (default)")
    p.add_argument("--output-dir", default="validation_output", help="where the traces are written (default: ./validation_output)")
    p.add_argument("--kokkos-num-threads", type=int, default=6, help="host threads given to Kokkos (default: 6)")
    p.add_argument("--dg-cfl", type=float, default=None, dest="dg_cfl",
                   help="DG time step as a fraction of h / (vp * order^2), for every order "
                        f"(default: {cases.DG_CFL})")
    p.add_argument("--list", action="store_true", help="print the selected cases and exit")
    args = p.parse_args()

    if args.dg_cfl is not None:
        cases.DG_CFL = {o: args.dg_cfl for o in cases.ORDERS}
    selected = select_cases(args)
    if args.list or not selected:
        for c in selected:
            print(c.name)
        return 0

    for var in ("KOKKOS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT"):
        os.environ[var] = str(args.kokkos_num_threads)
    import kokkos

    results, failures = [], []
    kokkos.initialize()
    try:
        for i, case in enumerate(selected, 1):
            print(f"\n[{i}/{len(selected)}] {case.name}: {case.elements_per_direction}^3 elements, "
                  f"{case.n_samples} steps of {case.dt:g} s", flush=True)
            try:
                r = validate_case(case, args.output_dir)
            except Exception:
                traceback.print_exc()
                failures.append(case)
                continue
            results.append(r)
            for c, e in r.errors.items():
                print(f"  {c:<22} L2 rel {e:.3e}")
            print(f"  total L2 rel {r.total_error:.3e} (tolerance {case.tolerance:g}) "
                  f"-> {'PASS' if r.passed else 'FAIL'}  [{r.elapsed:.1f} s]", flush=True)
    finally:
        # Every Kokkos view must be released before finalize.
        gc.collect()
        kokkos.finalize()

    print_summary(results, failures)
    return 0 if not failures and all(r.passed for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
