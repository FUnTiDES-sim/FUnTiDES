#!/usr/bin/env python3
"""
FUnTiDES: Multi-GPU Acoustic Wave Simulation
============================================

This example demonstrates a distributed multi-GPU simulation.
It uses:
 1. MPI4Py for inter-process communication.
 2. Kokkos for on-device computation.
 3. C++ Partitioner to handle domain decomposition automatically.
 4. Matplotlib for visualization.

Usage:
  mpirun -n <N_GPUS> python examples/fe/solver_multigpu.py
"""

import argparse
import os
import sys
import time
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. MPI & GPU Setup (Before importing Kokkos)
# ------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Simple Round-Robin GPU assignment
# This works for both CUDA and HIP without vendor-specific headers
local_rank = int(os.environ.get("OMP_COMM_WORLD_LOCAL_RANK", rank % 4))
os.environ["KOKKOS_DEVICE_ID"] = str(local_rank)
os.environ["CUDA_VISIBLE_DEVICES"] = str(local_rank)
os.environ["HIP_VISIBLE_DEVICES"] = str(local_rank)

# Import bindings after env setup
import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver

# Type aliases
Partitioner = Model.CartesianXPartitioner_f32_i32
Params = Model.CartesianParams_f32_i32
Builder = Model.CartesianStructBuilder_f32_i32_O2 # Order 2

# ------------------------------------------------------------------------------
# 2. Helper Classes
# ------------------------------------------------------------------------------

class GhostExchanger:
    """
    Handles synchronization of boundary values between MPI ranks.
    Assumes 1D decomposition along X.
    """
    def __init__(self, nx, ny, nz):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        # Pre-calculate indices for Left (X=0) and Right (X=nx-1) faces
        self.face_indices_left = []
        self.face_indices_right = []

        nx_val = int(nx)
        ny_val = int(ny)
        nz_val = int(nz)

        # Loop order matches standard flattening (but indices are what matters)
        for k in range(nz_val):
            for j in range(ny_val):
                # Flat index = i + j*nx + k*nx*ny (X-Fastest / F-Order)

                # Left face: i = 0
                idx_l = 0 + j * nx_val + k * nx_val * ny_val

                # Right face: i = nx - 1
                idx_r = (nx_val - 1) + j * nx_val + k * nx_val * ny_val

                self.face_indices_left.append(idx_l)
                self.face_indices_right.append(idx_r)

        self.face_indices_left = np.array(self.face_indices_left, dtype=np.int32)
        self.face_indices_right = np.array(self.face_indices_right, dtype=np.int32)
        self.buf_size = len(self.face_indices_left)

    def sync_accumulate(self, field_array):
        """
        Exchange ghost values and ADD them to local boundaries.
        Used for Mass Matrix assembly and Force assembly.
        """
        # 1. Right -> Next Left
        if self.rank < self.size - 1:
            send_buf = field_array[self.face_indices_right].copy()
            recv_buf = np.zeros(self.buf_size, dtype=np.float32)

            self.comm.Sendrecv(sendbuf=send_buf, dest=self.rank+1,
                               recvbuf=recv_buf, source=self.rank+1)

            field_array[self.face_indices_right] += recv_buf

        # 2. Left -> Prev Right
        if self.rank > 0:
            send_buf = field_array[self.face_indices_left].copy()
            recv_buf = np.zeros(self.buf_size, dtype=np.float32)

            self.comm.Sendrecv(sendbuf=send_buf, dest=self.rank-1,
                               recvbuf=recv_buf, source=self.rank-1)

            field_array[self.face_indices_left] += recv_buf

def save_snapshot_distributed(field, nx, ny, nz, t, output_dir):
    """
    Gathers data from all GPUs to Rank 0, stitches it, and saves a PNG plot.
    """
    # 1. Reshape 1D array to 3D with F-order (X-fastest)
    # Shape becomes (nx, ny, nz)
    try:
        local_3d = field.reshape((nx, ny, nz), order='F')
    except Exception as e:
        if rank == 0: print(f"Reshape failed: {e}")
        return

    # 2. Trim ghost/overlap layers
    # Rank 0 keeps everything. Ranks > 0 skip the first X plane (i=0)
    # because it represents the same physical nodes as Rank i-1's last plane.
    if rank == 0:
        data_to_send = local_3d
    else:
        data_to_send = local_3d[1:, :, :] # Skip X=0

    # 3. Gather blocks to Rank 0
    # We use simple gather() which pickles numpy arrays.
    # While slightly slower than Gatherv for huge data, it simplifies 3D stitching significantly.
    blocks = comm.gather(data_to_send, root=0)

    # 4. Stitch and Plot (Rank 0)
    if rank == 0:
        # Stitch along X axis (axis 0)
        full_field = np.concatenate(blocks, axis=0) # Result shape: (GlobalX, Y, Z)

        # Extract middle Z slice for XY plotting
        mid_z = full_field.shape[2] // 2
        slice_xy = full_field[:, :, mid_z]

        # Transpose for plotting:
        # imshow expects (Rows/Y, Cols/X). Our array is (X, Y).
        slice_plot = slice_xy.T

        plt.figure(figsize=(10, 4))

        # Dynamic color scaling based on wave amplitude
        vm = np.max(np.abs(slice_plot))
        if vm < 1e-6: vm = 1e-6

        plt.imshow(slice_plot, cmap='seismic', vmin=-vm, vmax=vm, origin='lower', aspect='auto')
        plt.colorbar(label="Pressure")
        plt.title(f"Wavefield at Step {t}")
        plt.xlabel("X (Global)")
        plt.ylabel("Y")

        filename = os.path.join(output_dir, f"snap_{t:05d}.png")
        plt.savefig(filename, bbox_inches='tight', dpi=100)
        plt.close()
        print(f"  -> Saved {filename}")

# ------------------------------------------------------------------------------
# 3. Main Simulation Logic
# ------------------------------------------------------------------------------

def run_simulation(args):
    """
    Runs the simulation loop.
    All Kokkos objects created here will go out of scope when this function returns.
    """

    # --- Detect Memory Space & Layout ---
    layout = kokkos.LayoutRight
    mem_space = kokkos.HostSpace

    if hasattr(kokkos, 'CudaUVMSpace'):
        mem_space = kokkos.CudaUVMSpace
        layout = kokkos.LayoutLeft
    elif hasattr(kokkos, 'HipSpace'):
        mem_space = kokkos.HipSpace
        layout = kokkos.LayoutLeft

    if rank == 0:
        print(f"--- Multi-GPU Simulation ({size} ranks) ---")
        print(f"Global Grid: {args.ex} x {args.ey} x {args.ez}")
        print(f"Memory: {mem_space}, Layout: {layout}")
        os.makedirs("snapshots", exist_ok=True)

    # --- A. Domain Decomposition ---
    global_params = Params()
    global_params.order = 2
    global_params.ex = args.ex
    global_params.ey = args.ey
    global_params.ez = args.ez
    global_params.lx = 2000.0
    global_params.ly = 500.0
    global_params.lz = 500.0
    global_params.is_model_on_nodes = False
    global_params.is_elastic = False

    partitioner = Partitioner()
    local_params = partitioner.partition(global_params, rank, size)

    # --- B. Build Mesh & Solver ---
    builder = Builder(
        local_params.ex, local_params.lx,
        local_params.ey, local_params.ly,
        local_params.ez, local_params.lz,
        False, False,
        local_params.origin_x,
        local_params.origin_y,
        local_params.origin_z
    )
    mesh = builder.get_model()

    solver = Solver.create_solver(
        Solver.MethodType.SEM,
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.STRUCT,
        Solver.ModelLocationType.ONELEMENTS,
        Solver.PhysicType.ACOUSTIC,
        2
    )

    # --- C. Initialization & Assembly ---
    solver.compute_fe_init(mesh)
    mass_matrix = solver.get_mass_matrix()

    nx = int(local_params.ex * 2 + 1)
    ny = int(local_params.ey * 2 + 1)
    nz = int(local_params.ez * 2 + 1)
    n_dof = int(mesh.get_number_of_nodes())
    n_elem = int(mesh.get_number_of_elements())

    exchanger = GhostExchanger(nx, ny, nz)
    exchanger.sync_accumulate(mass_matrix)

    # --- D. Allocations ---
    # Create Kokkos Views with explicit layout to match C++ bindings
    kk_prev = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)
    kk_curr = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)

    np.array(kk_prev, copy=False)[:] = 0.0
    np.array(kk_curr, copy=False)[:] = 0.0

    wavefield = Solver.WavefieldAcoustic(kk_prev, kk_curr)

    # Source Terms
    n_rhs = 1
    kk_rhs_t = kokkos.array([n_rhs, args.steps], dtype=kokkos.float32, space=mem_space, layout=layout)
    kk_rhs_e = kokkos.array([n_rhs], dtype=kokkos.int32, space=mem_space, layout=layout)
    kk_rhs_w = kokkos.array([n_rhs, 27], dtype=kokkos.float32, space=mem_space, layout=layout)

    # Setup source on Rank 0
    if rank == 0:
        t_axis = np.linspace(0, 1.0, args.steps)
        f0 = 5.0
        src_val = (1 - 2*(np.pi*f0*(t_axis-0.2))**2) * np.exp(-(np.pi*f0*(t_axis-0.2))**2)

        np_rhs_t = np.array(kk_rhs_t, copy=False)
        np_rhs_t[0, :] = src_val.astype(np.float32)

        np_rhs_e = np.array(kk_rhs_e, copy=False)
        # Place source in the middle of Rank 0
        np_rhs_e[0] = n_elem // 2

        np_rhs_w = np.array(kk_rhs_w, copy=False)
        np_rhs_w[0, :] = 1.0 / 27.0
    else:
        np.array(kk_rhs_w, copy=False)[:] = 0.0

    rhs = Solver.RhsAcoustic(kk_rhs_t, kk_rhs_e, kk_rhs_w)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    # --- E. Time Loop ---
    dt = 0.001
    comm.Barrier()

    if rank == 0:
        print(f"Starting simulation for {args.steps} steps...")
        start_t = time.time()

    for t in range(args.steps):
        # 1. Compute Local
        solver.compute_forces(dt, t, data)

        # 2. Sync Boundaries
        force_vec = solver.get_force_vector(0)
        exchanger.sync_accumulate(force_vec)

        # 3. Update
        solver.update_solution_forward(dt, data)

        # 4. Snapshot
        if t % 100 == 0:
            if rank == 0: print(f"Step {t} / {args.steps}")
            active_field = np.array(kk_curr, copy=False)
            save_snapshot_distributed(active_field, nx, ny, nz, t, "snapshots")

        # 5. Swap
        data.swap_wavefields()

    comm.Barrier()
    if rank == 0:
        dur = time.time() - start_t
        print(f"Done in {dur:.2f}s")

# ------------------------------------------------------------------------------
# 4. Init / Finalize Wrapper
# ------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ex", type=int, default=200, help="Global Elements X")
    parser.add_argument("--ey", type=int, default=50,  help="Global Elements Y")
    parser.add_argument("--ez", type=int, default=50,  help="Global Elements Z")
    # Updated default to 1500
    parser.add_argument("--steps", type=int, default=1500, help="Time steps")
    args = parser.parse_args()

    # Initialize Kokkos
    kokkos.initialize()
    try:
        # Run simulation in separate function to ensure Python variables
        # holding Kokkos Views are destroyed *before* finalize is called.
        run_simulation(args)
    finally:
        # Finalize Kokkos
        kokkos.finalize()

if __name__ == "__main__":
    main()
