#!/usr/bin/env python3
"""
FUnTiDES: MPI + Python + Kokkos Solver
======================================

Runs the Acoustic SEM solver distributed across N MPI ranks.
Features:
- 1D Domain Decomposition (X-axis).
- 2 Ricker sources placed at 25% and 75% of the domain.
- Visualization snapshots gathered to Rank 0.

Usage:
  mpirun -n 4 python3 examples/fe/solver_cartesian_mpi.py --ex 200 --snap_interval 20
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
# Pin 1 thread per rank to avoid OpenMP oversubscription
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
    if memspace_arg == MemSpace.GPU.name and hasattr(kokkos, 'CudaUVMSpace'):
        return kokkos.CudaUVMSpace, kokkos.LayoutRight
    if memspace_arg == MemSpace.GPU.name and rank == 0:
        print("WARNING: CudaUVMSpace requested but not found. Falling back to HostSpace.")
    return kokkos.HostSpace, kokkos.LayoutRight

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mem", choices=[e.name for e in MemSpace], default="CPU")
    parser.add_argument("--ex", type=int, default=100, help="Total elements in X")
    parser.add_argument("--ey", type=int, default=50, help="Elements in Y")
    parser.add_argument("--ez", type=int, default=50, help="Elements in Z")
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--n_time_steps", type=int, default=1000)
    parser.add_argument("--snap_interval", type=int, default=50)
    parser.add_argument("--output_dir", type=str, default="snapshots")
    return parser.parse_args()

def get_partition_1d(global_ex, global_lx, rank, size):
    """
    Splits the global X elements among MPI ranks.
    Returns: local_ex, local_lx, local_origin_x, global_offset_elem
    """
    base = global_ex // size
    rem = global_ex % size

    # Ranks < rem get 1 extra element to balance load
    local_ex = base + (1 if rank < rem else 0)

    # Calculate offset
    offset = rank * base + min(rank, rem)

    dx = global_lx / global_ex
    local_lx = local_ex * dx
    origin_x = offset * dx

    return local_ex, local_lx, origin_x, offset

class BoundaryComm:
    """Synchronizes boundary values (Sum/Assembly) for 1D X-decomposition."""
    def __init__(self, nx, ny, nz):
        # Indices for Left Face (local index 0 for X)
        self.idx_left = [j*nx + k*nx*ny for k in range(nz) for j in range(ny)]
        # Indices for Right Face (local index nx-1 for X)
        self.idx_right = [(nx-1) + i for i in self.idx_left]

        self.idx_left = np.array(self.idx_left, dtype=np.int32)
        self.idx_right = np.array(self.idx_right, dtype=np.int32)
        self.buf_size = len(self.idx_left)

    def sync_accumulate(self, field_view):
        """
        Exchange and sum boundary values.
        1. Rank i sends Right -> Rank i+1 receives Left.
        2. Rank i sends Left -> Rank i-1 receives Right.
        """
        # field_view is a numpy wrapper around C++ memory

        # Phase 1: Exchange with Right Neighbor (Rank i <-> Rank i+1)
        if rank < size - 1:
            recv = np.zeros(self.buf_size, dtype=np.float32)
            # Send my Right face, Receive neighbor's Left face
            comm.Sendrecv(sendbuf=field_view[self.idx_right].copy(), dest=rank+1,
                          recvbuf=recv, source=rank+1)
            field_view[self.idx_right] += recv # Accumulate

        # Phase 2: Exchange with Left Neighbor (Rank i <-> Rank i-1)
        if rank > 0:
            recv = np.zeros(self.buf_size, dtype=np.float32)
            # Send my Left face, Receive neighbor's Right face
            comm.Sendrecv(sendbuf=field_view[self.idx_left].copy(), dest=rank-1,
                          recvbuf=recv, source=rank-1)
            field_view[self.idx_left] += recv # Accumulate

def save_snapshot(field_view, local_nx, ny, nz, t, output_dir):
    """
    Gathers distributed fields to Rank 0, stitches them, and saves an image.
    """
    # 1. Get local data as numpy array
    data_1d = np.array(field_view, copy=False)

    try:
        # Reshape to 3D (assuming C++ LayoutRight / RowMajor mapping)
        local_3d = data_1d.reshape((nz, ny, local_nx))
    except ValueError:
        return

    # 2. Trim Overlap
    # In SEM, the last node of Rank i is the same physical point as the first node of Rank i+1.
    # To stitch images, Ranks > 0 skip their first X-plane.
    if rank == 0:
        send_data = local_3d
    else:
        send_data = local_3d[:, :, 1:]

    # 3. Gather all blocks to Rank 0
    gathered_blocks = comm.gather(send_data, root=0)

    # 4. Plot (Rank 0 only)
    if rank == 0:
        # Concatenate along X (axis 2)
        global_field = np.concatenate(gathered_blocks, axis=2)

        # Slice XZ plane at middle Y
        mid_y = ny // 2
        slice_xz = global_field[:, mid_y, :]

        plt.figure(figsize=(10, 5))
        plt.imshow(slice_xz, origin='lower', cmap='seismic', aspect='auto')

        # Normalize color scale to see waves clearly
        vm = np.max(np.abs(slice_xz))
        if vm > 1e-6:
            plt.clim(-vm, vm)

        plt.colorbar(label='Pressure')
        plt.title(f"Time Step {t} (XZ Slice)")
        plt.xlabel("X (Nodes)")
        plt.ylabel("Z (Nodes)")

        filename = os.path.join(output_dir, f"snap_{t:05d}.png")
        plt.savefig(filename, dpi=100)
        plt.close()

def main():
    args = parse_args()
    kokkos.initialize()

    mem_space, layout = select_kokkos_memspace(args.mem)

    if rank == 0:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"--- FUnTiDES MPI Solver ({size} ranks) ---")
        print(f"Memory Space: {mem_space}")
        print(f"Grid: {args.ex}x{args.ey}x{args.ez}")

    try:
        # 1. Domain Decomposition
        local_ex, local_lx, origin_x, offset_elem = get_partition_1d(args.ex, 1500.0, rank, size)

        # Local dimensions (SEM Order 2 -> 2 nodes per element + 1 shared)
        local_nx = local_ex * 2 + 1
        ny = args.ey * 2 + 1
        nz = args.ez * 2 + 1
        n_dof = local_nx * ny * nz

        # Print Local Partition Info
        print(f"[Rank {rank}] Elements: {local_ex}, Nodes: {local_nx}x{ny}x{nz}, Offset: {offset_elem}")
        comm.Barrier()

        # 2. Build Model & Solver
        builder = CartesianStructBuilderFI2(local_ex, local_lx, args.ey, 1500.0, args.ez, 1500.0,
                                            False, False, origin_x, 0.0, 0.0)
        model = builder.get_model()

        solver = Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU,
                                      Solver.MeshType.STRUCT, Solver.ModelLocationType.ONELEMENTS,
                                      Solver.PhysicType.ACOUSTIC, 2)

        # Init internal structures
        solver.compute_fe_init(model, [0.0]*3, False, 0.0)

        # 3. Physics Assembly (Mass Matrix)
        # Synchronize mass matrix boundaries once at startup
        comm_sync = BoundaryComm(local_nx, ny, nz)
        mass_matrix = solver.get_mass_matrix()
        comm_sync.sync_accumulate(mass_matrix)

        # 4. Allocations
        kk_prev = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)
        kk_curr = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)
        np.array(kk_prev, copy=False)[:] = 0.0
        np.array(kk_curr, copy=False)[:] = 0.0

        wavefield = Solver.WavefieldAcoustic(kk_prev, kk_curr)

        # 5. Source Setup (2 Sources Distributed)
        n_sources = 2
        n_steps = args.n_time_steps
        kk_rhs_t = kokkos.array([n_sources, n_steps], dtype=kokkos.float32, space=mem_space, layout=layout)
        kk_rhs_e = kokkos.array([n_sources], dtype=kokkos.int32, space=mem_space, layout=layout)
        kk_rhs_w = kokkos.array([n_sources, 27], dtype=kokkos.float32, space=mem_space, layout=layout)

        # Generate Source Signal (Ricker)
        t_vals = np.linspace(0, n_steps*args.dt, n_steps)
        f0 = 5.0
        src_sig = (1 - 2*(np.pi*f0*(t_vals-1.2/f0))**2) * np.exp(-(np.pi*f0*(t_vals-1.2/f0))**2)
        np.array(kk_rhs_t, copy=False)[:] = src_sig

        # Define Source Locations (Global Element Indices)
        # Source 1: 25% of X
        # Source 2: 75% of X
        src_global_locations = [args.ex // 4, 3 * args.ex // 4]

        np_rhs_e = np.array(kk_rhs_e, copy=False)
        np_rhs_w = np.array(kk_rhs_w, copy=False)

        # Determine if sources are on this rank
        for i, src_glob in enumerate(src_global_locations):
            if offset_elem <= src_glob < (offset_elem + local_ex):
                # Calculate local element index
                local_elem_x = src_glob - offset_elem
                # Flatten 3D index: x + y*ex + z*ex*ey
                elem_idx = local_elem_x + (args.ey//2)*local_ex + (args.ez//2)*local_ex*args.ey

                np_rhs_e[i] = elem_idx
                np_rhs_w[i, :] = 1.0/27.0 # Distribute weight
            else:
                # Not on this rank
                np_rhs_e[i] = 0
                np_rhs_w[i, :] = 0.0

        rhs = Solver.RhsAcoustic(kk_rhs_t, kk_rhs_e, kk_rhs_w)
        data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

        # 6. Time Loop
        comm.Barrier()
        if rank == 0: print("Starting Simulation...")

        start_time = time.time()
        for t in range(n_steps):
            # A. Compute Local Forces
            solver.compute_forces(args.dt, t, data)

            # B. Synchronize Forces (Assembly)
            force_vec = solver.get_force_vector(0)
            comm_sync.sync_accumulate(force_vec)

            # C. Update Solution
            solver.update_solution_forward(args.dt, data)

            # D. Snapshot
            if t % args.snap_interval == 0:
                # Determine which buffer holds the valid data (leapfrog swap logic)
                # t=0 computes into kk_prev (passed as 'prev' to C++, which writes to 'prev' buffer for next step)
                # Then swap happens.
                # Just saving one consistent view (kk_prev or kk_curr) is usually enough for visualization
                # For strict correctness, we track the active buffer.
                active = kk_prev if (t % 2 == 0) else kk_curr
                save_snapshot(active, local_nx, ny, nz, t, args.output_dir)
                if rank == 0: print(f"Step {t} saved")

            # E. Swap
            data.swap_wavefields()

        if rank == 0:
            print(f"Done. Total time: {time.time() - start_time:.2f}s")

    finally:
        # Cleanup
        try:
            del data, wavefield, rhs, solver, model
            del kk_prev, kk_curr, kk_rhs_t, kk_rhs_e, kk_rhs_w
            if 'mass_matrix' in locals(): del mass_matrix
            if 'force_vec' in locals(): del force_vec
        except:
            pass
        kokkos.finalize()

if __name__ == "__main__":
    main()
