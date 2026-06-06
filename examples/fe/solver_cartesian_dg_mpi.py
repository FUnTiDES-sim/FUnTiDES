#!/usr/bin/env python3
"""
solver_cartesian_dg_mpi.py — DG acoustic solver, MPI 1D X-decomposition.

Ghost element approach: each rank owns [n_local + 2*n_ghost, n_dof] wavefield rows.
Ghost rows are filled via MPI Sendrecv before each compute_forces call.
Verlet update is restricted to rows 0..n_local-1 via set_n_local_elem.

Note: CartesianUnstructBuilder always starts at origin (0,0,0).
      For a homogeneous medium this is fine (all elements identical geometry).

Usage:
  mpirun -n 4 python solver_cartesian_dg_mpi.py --ex 40 --ey 10 --ez 10
"""

import argparse
import gc
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

_total_threads = int(os.environ.get("TOTAL_OMP_THREADS", "6"))
_ranks = int(os.environ.get("OMPI_COMM_WORLD_SIZE", os.environ.get("PMI_SIZE", "1")))
_threads_per_rank = max(1, _total_threads // _ranks)
os.environ.setdefault("OMP_NUM_THREADS", str(_threads_per_rank))
os.environ.setdefault("OMP_PROC_BIND", "false")
os.environ.setdefault("KOKKOS_NUM_THREADS", str(_threads_per_rank))

import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

CartesianUnstructBuilder = Model.CartesianUnstructBuilder_f32_i32
CartesianParams = Model.CartesianParams_f32_i32


def parse_args():
    parser = argparse.ArgumentParser(
        description="DG acoustic MPI solver, 1D X-decomposition."
    )
    parser.add_argument("--order", type=int, default=2, choices=range(1, 4))
    parser.add_argument("--ex", type=int, default=40, help="Total elements in X")
    parser.add_argument("--ey", type=int, default=10)
    parser.add_argument("--ez", type=int, default=10)
    parser.add_argument("--domain_size", type=float, default=1500.0)
    parser.add_argument("--f0", type=float, default=5.0)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--n_time_steps", type=int, default=500)
    parser.add_argument("--snap_interval", type=int, default=50)
    parser.add_argument("--threads", type=int, default=0,
                        help="Kokkos threads per rank (0=auto: TOTAL_OMP_THREADS/nranks)")
    return parser.parse_args()


def get_partition_1d(global_ex, global_lx, rank, size):
    base = global_ex // size
    rem = global_ex % size
    local_ex = base + (1 if rank < rem else 0)
    offset_x = rank * base + min(rank, rem)
    dx = global_lx / global_ex
    return local_ex, local_ex * dx, offset_x * dx, offset_x


def source_term(time_n, f0):
    o_tpeak = 1.0 / f0
    if time_n <= -0.9 * o_tpeak or time_n >= 2.9 * o_tpeak:
        return 0.0
    pi = 3.14157
    lam = (f0 * pi) ** 2
    return (2.0 * lam * (2.0 * lam * (time_n - o_tpeak) ** 2 - 1.0)
            * np.exp(-lam * (time_n - o_tpeak) ** 2))


class DGBoundaryComm:
    """Exchange DG element DOF slices for 1D X-decomposition (COPY semantics)."""

    def __init__(self, ex_local, ey, ez, n_dof_per_elem, n_local):
        self.n_ghost = ey * ez
        self.n_local = n_local
        self.n_dof = n_dof_per_elem
        self.buf_size = self.n_ghost * n_dof_per_elem
        # Boundary element indices, ordered k-outer j-inner (must match ghost_elems in set_partition_faces_from_elems).
        self.left_elems = [j * ex_local + k * ex_local * ey
                           for k in range(ez) for j in range(ey)]
        self.right_elems = [(ex_local - 1) + j * ex_local + k * ex_local * ey
                            for k in range(ez) for j in range(ey)]

    def exchange(self, field_np):
        """
        Fill ghost rows from MPI neighbours (two paired Sendrecv calls).

        field_np: numpy view of shape [n_local + 2*n_ghost, n_dof_per_elem].
          rows n_local .. n_local+n_ghost-1         : left ghost  (rank-1's right slice)
          rows n_local+n_ghost .. n_local+2*n_ghost-1: right ghost (rank+1's left slice)
        """
        left_rank  = rank - 1 if rank > 0       else MPI.PROC_NULL
        right_rank = rank + 1 if rank < size - 1 else MPI.PROC_NULL
        n = self.n_ghost
        ghost_left_start  = self.n_local
        ghost_right_start = self.n_local + n

        # Phase 1: each rank sends RIGHT slice to right_rank, receives from left_rank → left ghost.
        send_right = field_np[self.right_elems, :].copy().ravel()
        recv_left  = np.empty(self.buf_size, dtype=np.float32)
        comm.Sendrecv(send_right, dest=right_rank, sendtag=0,
                      recvbuf=recv_left, source=left_rank, recvtag=0)
        if rank > 0:
            field_np[ghost_left_start:ghost_left_start + n, :] = recv_left.reshape(n, self.n_dof)

        # Phase 2: each rank sends LEFT slice to left_rank, receives from right_rank → right ghost.
        send_left  = field_np[self.left_elems, :].copy().ravel()
        recv_right = np.empty(self.buf_size, dtype=np.float32)
        comm.Sendrecv(send_left, dest=left_rank, sendtag=1,
                      recvbuf=recv_right, source=right_rank, recvtag=1)
        if rank < size - 1:
            field_np[ghost_right_start:ghost_right_start + n, :] = recv_right.reshape(n, self.n_dof)


def build_partition_elem_lists(ex_local, ey, ez, n_local, n_ghost):
    """
    Build element + ghost index lists for set_partition_faces_from_elems.

    Ordering: k-outer, j-inner (matches DGBoundaryComm ghost rows).
      left_elems[m]   = elem at (i=0,         j=m%ey, k=m//ey)
      right_elems[m]  = elem at (i=ex_local-1, j=m%ey, k=m//ey)
      left_ghosts[m]  = n_local + k*ey + j  (left ghost rows)
      right_ghosts[m] = n_local + n_ghost + k*ey + j  (right ghost rows)
    """
    left_elems   = [j * ex_local + k * ex_local * ey              for k in range(ez) for j in range(ey)]
    left_ghosts  = [n_local + k * ey + j                          for k in range(ez) for j in range(ey)]
    right_elems  = [(ex_local - 1) + j * ex_local + k * ex_local * ey for k in range(ez) for j in range(ey)]
    right_ghosts = [n_local + n_ghost + k * ey + j                for k in range(ez) for j in range(ey)]

    le = left_elems  if rank > 0       else []
    lg = left_ghosts if rank > 0       else []
    re = right_elems  if rank < size - 1 else []
    rg = right_ghosts if rank < size - 1 else []
    return le, lg, re, rg


def get_local_xz_slice(field_np, ex_local, ey, ez):
    """XZ slice at mid-Y from local rows only."""
    grid = np.zeros((ex_local, ez))
    y_mid = ey // 2
    for ix in range(ex_local):
        for iz in range(ez):
            e = ix + y_mid * ex_local + iz * ey * ex_local
            grid[ix, iz] = field_np[e, 0]
    return grid


def save_snapshot(field_np, ex_local, ey, ez, t):
    local_slice = get_local_xz_slice(field_np, ex_local, ey, ez)
    gathered = comm.gather(local_slice, root=0)
    if rank == 0:
        global_grid = np.concatenate(gathered, axis=0)
        plt.figure(figsize=(6, 6))
        vm = max(np.max(np.abs(global_grid)), 1e-10)
        plt.imshow(global_grid.T, origin='lower', cmap='seismic', aspect='equal',
                   vmin=-vm, vmax=vm)
        plt.colorbar(label='Pressure')
        plt.title(f"DG MPI — step {t}")
        plt.savefig(f"dg_mpi_snap{t:05d}.png", dpi=100)
        plt.close()


def run_simulation(args):
    """All Kokkos objects live here; function return frees them before finalize."""
    order = args.order
    global_ex, ey, ez = args.ex, args.ey, args.ez
    lx = ly = lz = args.domain_size
    n_dof_per_elem = (order + 1) ** 3
    dt = args.dt

    local_ex, local_lx, origin_x, offset_x = get_partition_1d(global_ex, lx, rank, size)
    n_local = local_ex * ey * ez
    n_ghost = ey * ez
    n_total_elem = n_local + 2 * n_ghost

    comm.Barrier()
    print(f"[Rank {rank}] ex_local={local_ex}, origin_x={origin_x:.1f}, "
          f"n_local={n_local}, n_total_elem={n_total_elem}")
    comm.Barrier()
    if rank == 0:
        print(f"--- FUnTiDES DG MPI ({size} ranks) "
              f"grid={global_ex}x{ey}x{ez} order={order} ---")

    params = CartesianParams()
    params.ex, params.ey, params.ez = local_ex, ey, ez
    params.lx, params.ly, params.lz = local_lx, ly, lz
    params.order = order
    params.is_model_on_nodes = False
    params.is_elastic = False
    model = CartesianUnstructBuilder(params).get_model()

    solver = Solver.create_solver(
        Solver.MethodType.DG, Solver.ImplemType.MAKUTU,
        Solver.MeshType.UNSTRUCT, Solver.ModelLocationType.ONELEMENTS,
        Solver.PhysicType.ACOUSTIC, order,
    )
    solver.compute_fe_init(model)

    le, lg, re, rg = build_partition_elem_lists(local_ex, ey, ez, n_local, n_ghost)
    solver.set_partition_faces_from_elems(le, lg, re, rg)
    solver.set_n_local_elem(n_local)

    kk_prev = kokkos.array([n_total_elem, n_dof_per_elem], dtype=kokkos.float32,
                            space=kokkos.HostSpace, layout=kokkos.LayoutRight)
    kk_curr = kokkos.array([n_total_elem, n_dof_per_elem], dtype=kokkos.float32,
                            space=kokkos.HostSpace, layout=kokkos.LayoutRight)
    np.array(kk_prev, copy=False)[:] = 0.0
    np.array(kk_curr, copy=False)[:] = 0.0

    kk_term   = kokkos.array([1, args.n_time_steps], dtype=kokkos.float32,
                              space=kokkos.HostSpace, layout=kokkos.LayoutRight)
    kk_elem   = kokkos.array([1], dtype=kokkos.int32,
                              space=kokkos.HostSpace, layout=kokkos.LayoutRight)
    kk_weight = kokkos.array([1, n_dof_per_elem], dtype=kokkos.float32,
                              space=kokkos.HostSpace, layout=kokkos.LayoutRight)

    term_np = np.array(kk_term, copy=False)
    for i in range(args.n_time_steps):
        term_np[0, i] = source_term(i * dt, args.f0)

    np_elem   = np.array(kk_elem, copy=False)
    np_weight = np.array(kk_weight, copy=False)
    src_global_x = global_ex // 2
    if offset_x <= src_global_x < offset_x + local_ex:
        local_src_x = src_global_x - offset_x
        np_elem[0]   = local_src_x + (ey // 2) * local_ex + (ez // 2) * local_ex * ey
        np_weight[0, :] = 1.0 / n_dof_per_elem
    else:
        np_elem[0]   = 0
        np_weight[0, :] = 0.0  # no source on this rank

    wavefield = Solver.DGWavefieldAcoustic(kk_prev, kk_curr)
    rhs       = Solver.RhsAcoustic(kk_term, kk_elem, kk_weight)
    data      = Solver.DGsolverDataAcoustic(wavefield, rhs)

    # Track which numpy array is the "current" field (swaps each step).
    curr_np = np.array(kk_curr, copy=False)
    prev_np = np.array(kk_prev, copy=False)

    dg_comm = DGBoundaryComm(local_ex, ey, ez, n_dof_per_elem, n_local)

    comm.Barrier()
    if rank == 0:
        print("Starting simulation...")
    start = time.time()

    for t in range(args.n_time_steps):
        dg_comm.exchange(curr_np)            # fill ghost rows from neighbours
        solver.compute_forces(dt, t, data)   # flux reads ghost rows
        solver.update_solution(dt, data)     # Verlet on 0..n_local-1 only
        data.swap_wavefields()
        curr_np, prev_np = prev_np, curr_np  # track active field

        if t % 100 == 0 and rank == 0:
            print(f"  step {t}/{args.n_time_steps}")
        if t % args.snap_interval == 0:
            save_snapshot(curr_np, local_ex, ey, ez, t)

    comm.Barrier()
    if rank == 0:
        print(f"Done. Total: {time.time() - start:.2f}s")
    # All locals freed here on return.


def main():
    args = parse_args()
    if args.threads > 0:
        os.environ["OMP_NUM_THREADS"] = str(args.threads)
        os.environ["KOKKOS_NUM_THREADS"] = str(args.threads)
    kokkos.initialize()
    run_simulation(args)
    gc.collect()
    kokkos.finalize()


if __name__ == "__main__":
    main()
