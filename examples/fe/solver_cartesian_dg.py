#!/usr/bin/env python3
"""
solver_cartesian_dg.py — DG acoustic solver on a Cartesian mesh.

Fields are stored per-element: shape (n_elem, n_dof_per_elem).
For help run with the --help option.
"""

import argparse
import os
import time
from datetime import datetime
from enum import Enum

import matplotlib.pyplot as plt
import numpy as np
import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver

CartesianStructBuilderFI1 = Model.CartesianStructBuilder_f32_i32_O1
CartesianStructBuilderFI2 = Model.CartesianStructBuilder_f32_i32_O2
CartesianStructBuilderFI3 = Model.CartesianStructBuilder_f32_i32_O3
CartesianUnstructBuilder = Model.CartesianUnstructBuilder_f32_i32
CartesianParams = Model.CartesianParams_f32_i32

os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OMP_THREAD_LIMIT", "6")
os.environ.setdefault("KOKKOS_NUM_THREADS", "6")


class MemSpace(Enum):
    CPU = "HostSpace"
    GPU = "CudaUVMSpace"


class ModelType(Enum):
    STRUCTURED = "Structured"
    UNSTRUCTURED = "Unstructured"


def detect_default_memspace():
    try:
        doc = Solver.DGWavefieldAcoustic.__init__.__doc__
        if doc and ("CudaUVMSpace" in doc or "CudaSpace" in doc):
            return MemSpace.GPU.name
        if doc and "HostSpace" in doc:
            return MemSpace.CPU.name
    except Exception:
        pass
    return MemSpace.CPU.name


def parse_args():
    default_mem = detect_default_memspace()
    parser = argparse.ArgumentParser(
        description="Run DG acoustic solver on a Cartesian unstructured mesh."
    )
    parser.add_argument("--mem", choices=[e.name for e in MemSpace], default=default_mem)
    parser.add_argument(
        "--model", choices=[e.name for e in ModelType], default=ModelType.UNSTRUCTURED.name
    )
    parser.add_argument("--order", type=int, default=2, choices=range(1, 4))
    parser.add_argument("--domain_size", type=float, default=1500.0)
    parser.add_argument("--ex", type=int, default=10)
    parser.add_argument("--ey", type=int, default=10)
    parser.add_argument("--ez", type=int, default=10)
    parser.add_argument("--f0", type=float, default=5.0)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--n_time_steps", type=int, default=500)
    parser.add_argument("--n_rhs", type=int, default=1)
    parser.add_argument("--on_nodes", action="store_true", default=False)
    return parser.parse_args()


def select_kokkos_memspace(memspace_arg):
    enum_value = MemSpace[memspace_arg]
    if enum_value == MemSpace.CPU:
        return kokkos.HostSpace, kokkos.LayoutRight
    return kokkos.CudaUVMSpace, kokkos.LayoutLeft


def source_term(time_n, f0):
    o_tpeak = 1.0 / f0
    if time_n <= -0.9 * o_tpeak or time_n >= 2.9 * o_tpeak:
        return 0.0
    pi = 3.14157
    lam = (f0 * pi) ** 2
    return (
        2.0 * lam * (2.0 * lam * (time_n - o_tpeak) ** 2 - 1.0)
        * np.exp(-lam * (time_n - o_tpeak) ** 2)
    )


def create_model(model_type, ex, ey, ez, lx, ly, lz, order, on_nodes):
    enum_value = ModelType[model_type]
    if enum_value == ModelType.STRUCTURED:
        builders = {
            1: CartesianStructBuilderFI1,
            2: CartesianStructBuilderFI2,
            3: CartesianStructBuilderFI3,
        }
        if order not in builders:
            raise ValueError(f"Order {order} not supported for structured DG")
        return builders[order](ex, lx, ey, ly, ez, lz, on_nodes, False).get_model()
    else:
        params = CartesianParams()
        params.ex, params.ey, params.ez = ex, ey, ez
        params.lx, params.ly, params.lz = lx, ly, lz
        params.order = order
        params.is_model_on_nodes = on_nodes
        params.is_elastic = False
        return CartesianUnstructBuilder(params).get_model()


def allocate_pressure_dg(n_elem, n_dof_per_elem, memspace, layout):
    """Allocate 2D pressure arrays (n_elem, n_dof_per_elem) for DG."""
    kk_prev = kokkos.array(
        [n_elem, n_dof_per_elem], dtype=kokkos.float32, space=memspace, layout=layout
    )
    kk_curr = kokkos.array(
        [n_elem, n_dof_per_elem], dtype=kokkos.float32, space=memspace, layout=layout
    )
    np.array(kk_prev, copy=False)[:] = 0.0
    np.array(kk_curr, copy=False)[:] = 0.0
    return kk_prev, kk_curr


def allocate_rhs_term(n_rhs, n_time_steps, dt, f0, memspace, layout):
    kk_term = kokkos.array(
        [n_rhs, n_time_steps], dtype=kokkos.float32, space=memspace, layout=layout
    )
    term = np.array(kk_term, copy=False)
    for i in range(n_time_steps):
        for r in range(n_rhs):
            term[r, i] = source_term(i * dt, f0)
    return kk_term


def allocate_rhs_weights(n_rhs, n_dof_per_elem, memspace, layout):
    kk_w = kokkos.array(
        [n_rhs, n_dof_per_elem], dtype=kokkos.float32, space=memspace, layout=layout
    )
    w = np.array(kk_w, copy=False)
    w[:] = 1.0 / n_dof_per_elem
    return kk_w


def allocate_rhs_element(n_rhs, ex, ey, ez, memspace, layout):
    kk_e = kokkos.array([n_rhs], dtype=kokkos.int32, space=memspace, layout=layout)
    e = np.array(kk_e, copy=False)
    # Place source near center
    e[0] = int(ex / 2) + int(ey / 2) * ex + int(ez / 2) * ey * ex
    for r in range(1, n_rhs):
        e[r] = int(ex / 3) + int(ey / 2) * ex + int(ez / 2) * ey * ex
    return kk_e


def get_snapshot_dg(ex, ez, ey, field_np, dof_idx=0):
    """Extract 2D slice (xz plane at y=ey//2) from DG field (n_elem, n_dof)."""
    grid = np.zeros((ex, ez))
    y_mid = ey // 2
    for ix in range(ex):
        for iz in range(ez):
            e = ix + y_mid * ex + iz * ey * ex
            grid[ix, iz] = field_np[e, dof_idx]
    return grid


def setup_plot(ex, ez):
    grid = np.zeros((ex, ez))
    fig, ax = plt.subplots()
    im = ax.imshow(grid, cmap="viridis", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Pressure")
    plt.title("DG acoustic — xz slice (y=mid), dof 0")
    plt.xlabel("X"); plt.ylabel("Z")
    plt.ioff()
    return fig, ax, im


def main():
    args = parse_args()

    order = args.order
    ex, ey, ez = args.ex, args.ey, args.ez
    lx = ly = lz = args.domain_size
    n_elem = ex * ey * ez
    n_dof_per_elem = (order + 1) ** 3
    dt = args.dt
    n_time_steps = args.n_time_steps
    n_rhs = args.n_rhs
    on_nodes = args.on_nodes

    print("==========SIMULATION PARAMETERS==========")
    print(f"method            : DG")
    print(f"model             : {args.model}")
    print(f"order             : {order}")
    print(f"elements          : {ex}x{ey}x{ez} = {n_elem}")
    print(f"dof/elem          : {n_dof_per_elem}")
    print(f"total dof         : {n_elem * n_dof_per_elem}")
    print(f"dt                : {dt}")
    print(f"n_time_steps      : {n_time_steps}")
    print(f"n_rhs             : {n_rhs}")
    print(f"mem               : {args.mem}")
    print("=========================================")

    _, _, im = setup_plot(ex, ez)

    kokkos.initialize()
    memspace, layout = select_kokkos_memspace(args.mem)

    print("Creating model...")
    model = create_model(args.model, ex, ey, ez, lx, ly, lz, order, on_nodes)
    print("Model created")

    mesh_type = (
        Solver.MeshType.STRUCT
        if ModelType[args.model] == ModelType.STRUCTURED
        else Solver.MeshType.UNSTRUCT
    )
    print("Creating DG solver...")
    solver = Solver.create_solver(
        Solver.MethodType.DG,
        Solver.ImplemType.MAKUTU,
        mesh_type,
        Solver.ModelLocationType.ONELEMENTS if not on_nodes else Solver.ModelLocationType.ONNODES,
        Solver.PhysicType.ACOUSTIC,
        order,
    )
    print("Solver created")

    print("Initializing solver (compute_fe_init)...")
    solver.compute_fe_init(model)
    print("Solver initialized")

    print("Allocating arrays...")
    kk_prev, kk_curr = allocate_pressure_dg(n_elem, n_dof_per_elem, memspace, layout)
    kk_term   = allocate_rhs_term(n_rhs, n_time_steps, dt, args.f0, memspace, layout)
    kk_weight = allocate_rhs_weights(n_rhs, n_dof_per_elem, memspace, layout)
    kk_elem   = allocate_rhs_element(n_rhs, ex, ey, ez, memspace, layout)
    print("Arrays allocated")

    wavefield = Solver.DGWavefieldAcoustic(kk_prev, kk_curr)
    rhs       = Solver.RhsAcoustic(kk_term, kk_elem, kk_weight)
    data      = Solver.DGsolverDataAcoustic(wavefield, rhs)
    data.print()

    pn_np = np.array(kk_curr, copy=False)

    start = time.time()
    iter_times = []

    for t in range(n_time_steps):
        t0 = time.time()
        solver.compute_forces(dt, t, data)
        solver.update_solution(dt, data)
        iter_times.append(time.time() - t0)

        data.swap_wavefields()

        if t % 100 == 0:
            print(f"  step {t}/{n_time_steps}")
        if t % 50 == 0:
            grid = get_snapshot_dg(ex, ez, ey, pn_np)
            im.set_array(grid)
            plt.draw()
            plt.ioff()
            plt.savefig(f"dg_snap{t:05d}.png")

    total = time.time() - start
    print("==========SIMULATION STATISTICS==========")
    print(f"Total time       : {total:.2f} s")
    print(f"Avg iter time    : {np.mean(iter_times):.4f} s")
    print(f"Min iter time    : {np.min(iter_times):.4f} s")
    print(f"Max iter time    : {np.max(iter_times):.4f} s")
    print("=========================================")

    del pn_np
    del data
    del wavefield, rhs
    del solver, model
    del kk_prev, kk_curr, kk_term, kk_weight, kk_elem
    kokkos.finalize()
    print("Done")


if __name__ == "__main__":
    main()
