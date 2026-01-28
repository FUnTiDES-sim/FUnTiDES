#!/usr/bin/env python3
"""
FUnTiDES: Multi-GPU & Multi-Node Acoustic Wave Simulation
=========================================================

Run on 2 Nodes, 2 GPUs each (4 GPUs total):
  sbatch run_multinode.sbatch

"""

import argparse
import os
import sys
import time
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. MPI & GPU Setup
# ------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def get_local_rank():
    # Case 1: SLURM
    if "SLURM_LOCALID" in os.environ:
        return int(os.environ["SLURM_LOCALID"])
    # Case 2: OpenMPI
    if "OMPI_COMM_WORLD_LOCAL_RANK" in os.environ:
        return int(os.environ["OMPI_COMM_WORLD_LOCAL_RANK"])
    # Case 3: MVAPICH
    if "MV2_COMM_WORLD_LOCAL_RANK" in os.environ:
        return int(os.environ["MV2_COMM_WORLD_LOCAL_RANK"])
    # Fallback
    return rank % 2 

local_rank = get_local_rank()

os.environ["KOKKOS_DEVICE_ID"] = str(local_rank)
os.environ["CUDA_VISIBLE_DEVICES"] = str(local_rank)
os.environ["HIP_VISIBLE_DEVICES"] = str(local_rank)

import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver

# Type aliases
Partitioner = Model.CartesianXPartitioner_f32_i32
Params = Model.CartesianParams_f32_i32
Builder = Model.CartesianStructBuilder_f32_i32_O2

# ------------------------------------------------------------------------------
# 2. Communication Helper
# ------------------------------------------------------------------------------
class GhostExchanger:
    def __init__(self, nx, ny, nz):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        
        self.face_indices_left = []
        self.face_indices_right = []
        
        nx_val = int(nx)
        ny_val = int(ny)
        nz_val = int(nz)

        # X-Fastest indexing
        for k in range(nz_val):
            for j in range(ny_val):
                idx_l = 0 + j * nx_val + k * nx_val * ny_val
                idx_r = (nx_val - 1) + j * nx_val + k * nx_val * ny_val
                self.face_indices_left.append(idx_l)
                self.face_indices_right.append(idx_r)
                
        self.face_indices_left = np.array(self.face_indices_left, dtype=np.int32)
        self.face_indices_right = np.array(self.face_indices_right, dtype=np.int32)

    def sync_accumulate(self, field_array):
        # 1. Right -> Next Left
        if self.rank < self.size - 1:
            send_buf = field_array[self.face_indices_right].copy()
            recv_buf = np.zeros(len(self.face_indices_right), dtype=np.float32)
            self.comm.Sendrecv(sendbuf=send_buf, dest=self.rank+1,
                               recvbuf=recv_buf, source=self.rank+1)
            field_array[self.face_indices_right] += recv_buf

        # 2. Left -> Prev Right
        if self.rank > 0:
            send_buf = field_array[self.face_indices_left].copy()
            recv_buf = np.zeros(len(self.face_indices_left), dtype=np.float32)
            self.comm.Sendrecv(sendbuf=send_buf, dest=self.rank-1,
                               recvbuf=recv_buf, source=self.rank-1)
            field_array[self.face_indices_left] += recv_buf

# ------------------------------------------------------------------------------
# 3. I/O Helper
# ------------------------------------------------------------------------------
def save_snapshot_distributed(field, nx, ny, nz, t, output_dir):
    t_io_start = time.time()
    
    # 1. Reshape 1D -> 3D (View)
    try:
        local_3d = field.reshape((nx, ny, nz), order='F')
    except Exception as e:
        if rank == 0: print(f"Reshape failed: {e}", flush=True)
        return 0.0

    # 2. Select valid region (Trim ghost/overlap)
    if rank == 0:
        valid_view = local_3d
    else:
        valid_view = local_3d[1:, :, :] # Skip X=0

    # 3. Extract 2D Slice Locally & COPY TO HOST
    mid_z = nz // 2
    local_slice_xy = np.array(valid_view[:, :, mid_z], copy=True)

    # 4. Gather 2D slices
    blocks = comm.gather(local_slice_xy, root=0)

    # 5. Stitch and Plot (Rank 0)
    if rank == 0:
        full_slice = np.concatenate(blocks, axis=0)
        slice_plot = full_slice.T 

        plt.figure(figsize=(12, 4))
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
        print(f"  -> Saved {filename}", flush=True)

    return time.time() - t_io_start

# ------------------------------------------------------------------------------
# 4. Main Simulation
# ------------------------------------------------------------------------------
def run_simulation(args):
    
    layout = kokkos.LayoutRight 
    mem_space = kokkos.HostSpace
    if hasattr(kokkos, 'CudaUVMSpace'):
        mem_space = kokkos.CudaUVMSpace
        layout = kokkos.LayoutLeft
    elif hasattr(kokkos, 'HipSpace'):
        mem_space = kokkos.HipSpace
        layout = kokkos.LayoutLeft 

    if rank == 0:
        print(f"--- FUnTiDES Multi-Node ({size} GPUs total) ---", flush=True)
        print(f"Global Elements: {args.ex} x {args.ey} x {args.ez}", flush=True)
        print(f"Memory Space:    {mem_space}", flush=True)
        os.makedirs("snapshots", exist_ok=True)

    # --- Domain Decomposition ---
    global_params = Params()
    global_params.order = 2
    global_params.ex = args.ex
    global_params.ey = args.ey
    global_params.ez = args.ez
    global_params.lx = 8000.0 
    global_params.ly = 2000.0
    global_params.lz = 2000.0
    global_params.is_model_on_nodes = False
    global_params.is_elastic = False

    partitioner = Partitioner()
    local_params = partitioner.partition(global_params, rank, size)
    
    if rank == 0:
        print("Partitioning & Topology setup...", flush=True)

    # --- Mesh & Solver ---
    builder = Builder(
        local_params.ex, local_params.lx,
        local_params.ey, local_params.ly,
        local_params.ez, local_params.lz,
        False, False, 
        local_params.origin_x, local_params.origin_y, local_params.origin_z
    )
    mesh = builder.get_model()

    solver = Solver.create_solver(
        Solver.MethodType.SEM, Solver.ImplemType.MAKUTU,
        Solver.MeshType.STRUCT, Solver.ModelLocationType.ONELEMENTS,
        Solver.PhysicType.ACOUSTIC, 2
    )

    # --- Init ---
    solver.compute_fe_init(mesh)
    mass_matrix = solver.get_mass_matrix()
    
    nx = int(local_params.ex * 2 + 1)
    ny = int(local_params.ey * 2 + 1)
    nz = int(local_params.ez * 2 + 1)
    n_dof = int(mesh.get_number_of_nodes())
    n_elem = int(mesh.get_number_of_elements())
    
    exchanger = GhostExchanger(nx, ny, nz)
    exchanger.sync_accumulate(mass_matrix)

    # --- Allocations ---
    kk_prev = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)
    kk_curr = kokkos.array([n_dof], dtype=kokkos.float32, space=mem_space, layout=layout)
    np.array(kk_prev, copy=False)[:] = 0.0
    np.array(kk_curr, copy=False)[:] = 0.0

    wavefield = Solver.WavefieldAcoustic(kk_prev, kk_curr)

    # --- Sources ---
    n_rhs = 1
    kk_rhs_t = kokkos.array([n_rhs, args.steps], dtype=kokkos.float32, space=mem_space, layout=layout)
    kk_rhs_e = kokkos.array([n_rhs], dtype=kokkos.int32, space=mem_space, layout=layout)
    kk_rhs_w = kokkos.array([n_rhs, 27], dtype=kokkos.float32, space=mem_space, layout=layout)

    t_axis = np.linspace(0, 1.0, args.steps)
    f0 = 5.0
    src_val = (1 - 2*(np.pi*f0*(t_axis-0.2))**2) * np.exp(-(np.pi*f0*(t_axis-0.2))**2)
    
    # 1. Time Function
    np_rhs_t = np.array(kk_rhs_t, copy=False)
    np_rhs_t[0, :] = src_val.astype(np.float32)
    
    # 2. Source Location (Explicit Geometric Center)
    # Get local element counts
    lex = local_params.ex
    ley = local_params.ey
    lez = local_params.ez
    
    # Calculate center indices (integers)
    cx = lex // 2
    cy = ley // 2
    cz = lez // 2
    
    # Convert to Linear Index (C++ Layout: X-Fastest)
    # index = i + j*nx + k*nx*ny
    center_idx = cx + cy * lex + cz * lex * ley
    
    np_rhs_e = np.array(kk_rhs_e, copy=False)
    np_rhs_e[0] = center_idx
    
    # 3. Weights
    np_rhs_w = np.array(kk_rhs_w, copy=False)
    np_rhs_w[0, :] = 1.0 / 27.0

    rhs = Solver.RhsAcoustic(kk_rhs_t, kk_rhs_e, kk_rhs_w)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    # --- Loop ---
    dt = 0.0005 
    total_compute_time = 0.0
    total_io_time = 0.0
    
    comm.Barrier()
    if rank == 0: 
        print(f"Starting simulation on {size} ranks...", flush=True)
        print(f"Injecting {size} sources (one per GPU) at local center.", flush=True)
        start_t = time.time()

    for t in range(args.steps):
        
        t_start = time.time()
        
        solver.compute_forces(dt, t, data)
        force_vec = solver.get_force_vector(0)
        exchanger.sync_accumulate(force_vec)
        solver.update_solution(dt, data)
        
        t_compute = time.time() - t_start
        total_compute_time += t_compute

        # --- I/O BLOCK ---
        if t % 100 == 0:
            if rank == 0: print(f"Step {t} / {args.steps}", flush=True)
            active_field = np.array(kk_curr, copy=False)
            t_io = save_snapshot_distributed(active_field, nx, ny, nz, t, "snapshots")
            total_io_time += t_io
        elif t % 10 == 0 and rank == 0:
            print(".", end="", flush=True)

        data.swap_wavefields()

    comm.Barrier()
    if rank == 0:
        total_time = total_compute_time + total_io_time
        print("\n" + "="*40, flush=True)
        print(f"BENCHMARK RESULTS ({size} GPUs)", flush=True)
        print(f"Grid: {args.ex}x{args.ey}x{args.ez}", flush=True)
        print("="*40, flush=True)
        print(f"Total GPU Compute Time : {total_compute_time:.4f} s", flush=True)
        print(f"Total I/O (Plot) Time  : {total_io_time:.4f} s", flush=True)
        print(f"Total Wall Time        : {total_time:.4f} s", flush=True)
        print(f"Avg Time per Step      : {total_compute_time/args.steps*1000:.2f} ms", flush=True)
        print("="*40, flush=True)

def main():
    parser = argparse.ArgumentParser()
    # Optimized grid size for P100 (enough work, but fits in memory)
    parser.add_argument("--ex", type=int, default=800, help="Global Elements X")
    parser.add_argument("--ey", type=int, default=300, help="Global Elements Y")
    parser.add_argument("--ez", type=int, default=300, help="Global Elements Z")
    parser.add_argument("--steps", type=int, default=1500, help="Time steps")
    args = parser.parse_args()

    kokkos.initialize()
    try:
        run_simulation(args)
    finally:
        kokkos.finalize()

if __name__ == "__main__":
    main()
