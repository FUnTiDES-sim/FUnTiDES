#!/usr/bin/env python3
"""
FUnTiDES: MPI + Python + Kokkos Solver (CPU/MPI Version)
========================================================

Usage:
  mpirun -n 2 python3 examples/fe/solver_cartesian_mpi.py --ex 100
"""

import argparse
import os
import sys
import time
from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

# Environment Setup
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OMP_PROC_BIND", "false")
os.environ.setdefault("KOKKOS_NUM_THREADS", "1")

import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

CartesianStructBuilderFI2 = Model.CartesianStructBuilder_f32_i32_O2
CartesianParams = Model.CartesianParams_f32_i32

class MemSpace(Enum):
    CPU = "HostSpace"
    GPU = "CudaUVMSpace"

def select_kokkos_memspace(memspace_arg):
    """
    Selects valid memory space. Falls back to CPU if GPU requested but unavailable.
    """
    # Check if GPU requested AND Kokkos was compiled with CUDA
    if memspace_arg == MemSpace.GPU.name and hasattr(kokkos, 'CudaUVMSpace'):
        return kokkos.CudaUVMSpace, kokkos.LayoutRight

    # Default to CPU
    if memspace_arg == MemSpace.GPU.name and rank == 0:
        print("WARNING: CudaUVMSpace requested but not found in Kokkos. Falling back to HostSpace.")

    return kokkos.HostSpace, kokkos.LayoutRight

def parse_args():
    parser = argparse.ArgumentParser()
    # Default to CPU if auto-detection fails or simple default needed
    parser.add_argument("--mem", choices=[e.name for e in MemSpace], default="CPU")
    parser.add_argument("--ex", type=int, default=100)
    parser.add_argument("--ey", type=int, default=50)
    parser.add_argument("--ez", type=int, default=50)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--n_time_steps", type=int, default=1000)
    return parser.parse_args()

def get_partition_1d(global_ex, global_lx, rank, size):
    base = global_ex // size
    rem = global_ex % size
    local_ex = base + (1 if rank < rem else 0)
    offset = rank * base + min(rank, rem)
    dx = global_lx / global_ex
    local_lx = local_ex * dx
    origin_x = offset * dx
    return local_ex, local_lx, origin_x, offset

class BoundaryComm:
    def __init__(self, nx, ny, nz):
        # 1D X-decomposition boundaries
        self.idx_left = [j*nx + k*nx*ny for k in range(nz) for j in range(ny)]
        self.idx_right = [(nx-1) + i for i in self.idx_left]
        self.idx_left = np.array(self.idx_left, dtype=np.int32)
        self.idx_right = np.array(self.idx_right, dtype=np.int32)
        self.buf_size = len(self.idx_left)

    def sync_accumulate(self, field_view):
        # field_view is now a numpy array (wrapping C++ memory)
        if rank < size - 1:
            recv = np.zeros(self.buf_size, dtype=np.float32)
            comm.Sendrecv(sendbuf=field_view[self.idx_right].copy(), dest=rank+1,
                          recvbuf=recv, source=rank+1)
            field_view[self.idx_right] += recv
        if rank > 0:
            recv = np.zeros(self.buf_size, dtype=np.float32)
            comm.Sendrecv(sendbuf=field_view[self.idx_left].copy(), dest=rank-1,
                          recvbuf=recv, source=rank-1)
            field_view[self.idx_left] += recv

def main():
    args = parse_args()
    kokkos.initialize()

    mem_space, layout = select_kokkos_memspace(args.mem)
    if rank == 0:
        print(f"Using Memory Space: {mem_space}")

    try:
        local_ex, local_lx, origin_x, offset_elem = get_partition_1d(args.ex, 1500.0, rank, size)
        local_nx = local_ex * 2 + 1
        ny = args.ey * 2 + 1
        nz = args.ez * 2 + 1
        n_dof = local_nx * ny * nz

        if rank == 0: print(f"Global: {args.ex}x{args.ey}x{args.ez}")
        print(f"[Rank {rank}] Local: {local_ex} elems, OriginX: {origin_x:.1f}")

        # Model
        builder = CartesianStructBuilderFI2(local_ex, local_lx, args.ey, 1500.0, args.ez, 1500.0,
                                            False, False, origin_x, 0.0, 0.0)
        model = builder.get_model()

        # Solver
        solver = Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU,
                                      Solver.MeshType.STRUCT, Solver.ModelLocationType.ONELEMENTS,
                                      Solver.PhysicType.ACOUSTIC, 2)
        solver.compute_fe_init(model, [0.0]*3, False, 0.0)

        # Mass Matrix Sync
        comm_sync = BoundaryComm(local_nx, ny, nz)
        mass_matrix = solver.get_mass_matrix()
        comm_sync.sync_accumulate(mass_matrix)

        # Allocations (Zero init)
        # Use kokkos.float32 instead of np.float32
        kk_prev = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)
        kk_curr = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)
        np.array(kk_prev, copy=False)[:] = 0.0
        np.array(kk_curr, copy=False)[:] = 0.0

        wavefield = Solver.WavefieldAcoustic(kk_prev, kk_curr)

        # RHS
        n_steps = args.n_time_steps
        # Use kokkos.float32 / kokkos.int32
        kk_rhs_t = kokkos.array([1, n_steps], dtype=kokkos.float32, space=mem_space, layout=layout)
        kk_rhs_e = kokkos.array([1], dtype=kokkos.int32, space=mem_space, layout=layout)
        kk_rhs_w = kokkos.array([1, 27], dtype=kokkos.float32, space=mem_space, layout=layout)

        # Source Signal
        t_vals = np.linspace(0, n_steps*args.dt, n_steps)
        f0 = 5.0
        src_sig = (1 - 2*(np.pi*f0*(t_vals-1.2/f0))**2) * np.exp(-(np.pi*f0*(t_vals-1.2/f0))**2)
        np.array(kk_rhs_t, copy=False)[0, :] = src_sig

        # Place source on Rank 0
        if rank == 0:
            np.array(kk_rhs_e, copy=False)[0] = local_ex//2 + (args.ey//2)*local_ex + (args.ez//2)*local_ex*args.ey
            np.array(kk_rhs_w, copy=False)[:] = 1.0/27.0
        else:
            np.array(kk_rhs_e, copy=False)[0] = 0
            np.array(kk_rhs_w, copy=False)[:] = 0.0

        rhs = Solver.RhsAcoustic(kk_rhs_t, kk_rhs_e, kk_rhs_w)
        data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

        comm.Barrier()
        if rank == 0: print("Starting Time Loop...")

        start_time = time.time()
        for t in range(n_steps):
            solver.compute_forces(args.dt, t, data)

            # Sync forces (C++ returns numpy wrapper via pybind11)
            force_vec = solver.get_force_vector(0)
            comm_sync.sync_accumulate(force_vec)

            solver.update_solution(args.dt, data)
            data.swap_wavefields()

            if rank == 0 and t % 100 == 0:
                print(f"Step {t}")

        if rank == 0:
            print(f"Done. Total time: {time.time() - start_time:.2f}s")

    finally:
        # Cleanup handles to ensure C++ destructors run before Finalize
        # Using simple try-except to avoid UnboundLocalError during cleanup if init failed
        try:
            del data
            del wavefield
            del rhs
            del solver
            del model
            del kk_prev
            del kk_curr
            del kk_rhs_t
            del kk_rhs_e
            del kk_rhs_w
            if 'mass_matrix' in locals(): del mass_matrix
            if 'force_vec' in locals(): del force_vec
        except:
            pass

        kokkos.finalize()

if __name__ == "__main__":
    main()
