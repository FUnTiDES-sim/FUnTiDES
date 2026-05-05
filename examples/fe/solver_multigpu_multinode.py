#!/usr/bin/env python3
"""
FUnTiDES: Multi-GPU & Multi-Node Acoustic Wave Simulation
=========================================================
"""

import argparse
import os
import sys
import time
import socket
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt

# Try importing CuPy for GPU-Direct support
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def get_local_rank():
    if "SLURM_LOCALID" in os.environ: return int(os.environ["SLURM_LOCALID"])
    if "OMPI_COMM_WORLD_LOCAL_RANK" in os.environ: return int(os.environ["OMPI_COMM_WORLD_LOCAL_RANK"])
    if "MV2_COMM_WORLD_LOCAL_RANK" in os.environ: return int(os.environ["MV2_COMM_WORLD_LOCAL_RANK"])
    return rank % 2

local_rank = get_local_rank()

os.environ["KOKKOS_DEVICE_ID"] = str(local_rank)
os.environ["CUDA_VISIBLE_DEVICES"] = str(local_rank)
os.environ["HIP_VISIBLE_DEVICES"] = str(local_rank)

import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver

Partitioner = Model.CartesianXPartitioner_f32_i32
Params = Model.CartesianParams_f32_i32
Builder = Model.CartesianStructBuilder_f32_i32_O2

class GhostExchanger:
    def __init__(self, nx, ny, nz):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        indices_left = []
        indices_right = []

        nx_val, ny_val, nz_val = int(nx), int(ny), int(nz)

        for k in range(nz_val):
            for j in range(ny_val):
                indices_left.append(0 + j * nx_val + k * nx_val * ny_val)
                indices_right.append((nx_val - 1) + j * nx_val + k * nx_val * ny_val)

        if HAS_CUPY:
            self.idx_l = cp.array(indices_left, dtype=np.int32)
            self.idx_r = cp.array(indices_right, dtype=np.int32)
        else:
            self.idx_l = np.array(indices_left, dtype=np.int32)
            self.idx_r = np.array(indices_right, dtype=np.int32)

        self.buf_size = len(indices_left)

    def sync_accumulate(self, field_array):
        if HAS_CUPY:
            d_field = cp.array(field_array, copy=False)
            if self.rank < self.size - 1:
                d_send_r = d_field[self.idx_r]
                d_recv_r = cp.empty_like(d_send_r)
                self.comm.Sendrecv(sendbuf=d_send_r, dest=self.rank+1, recvbuf=d_recv_r, source=self.rank+1)
                d_field[self.idx_r] += d_recv_r

            if self.rank > 0:
                d_send_l = d_field[self.idx_l]
                d_recv_l = cp.empty_like(d_send_l)
                self.comm.Sendrecv(sendbuf=d_send_l, dest=self.rank-1, recvbuf=d_recv_l, source=self.rank-1)
                d_field[self.idx_l] += d_recv_l
            cp.cuda.Stream.null.synchronize()
        else:
            if self.rank < self.size - 1:
                send_buf = field_array[self.idx_r].copy()
                recv_buf = np.zeros(self.buf_size, dtype=np.float32)
                self.comm.Sendrecv(sendbuf=send_buf, dest=self.rank+1, recvbuf=recv_buf, source=self.rank+1)
                field_array[self.idx_r] += recv_buf

            if self.rank > 0:
                send_buf = field_array[self.idx_l].copy()
                recv_buf = np.zeros(self.buf_size, dtype=np.float32)
                self.comm.Sendrecv(sendbuf=send_buf, dest=self.rank-1, recvbuf=recv_buf, source=self.rank-1)
                field_array[self.idx_l] += recv_buf

def save_snapshot_distributed(field, nx, ny, nz, t, output_dir):
    t_io_start = time.time()
    try:
        local_3d = field.reshape((nx, ny, nz), order='F')
    except Exception as e:
        if rank == 0: print(f"Reshape failed: {e}", flush=True)
        return 0.0

    if rank == 0:
        valid_view = local_3d
    else:
        valid_view = local_3d[1:, :, :]

    mid_z = nz // 2
    local_slice_xy = np.array(valid_view[:, :, mid_z], copy=True)
    blocks = comm.gather(local_slice_xy, root=0)

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

def run_simulation(args):

    # --- START INIT TIMER ---
    t_init_start = time.time()

    layout = kokkos.LayoutRight
    mem_space = kokkos.HostSpace
    if hasattr(kokkos, 'CudaUVMSpace'):
        mem_space = kokkos.CudaUVMSpace
        layout = kokkos.LayoutLeft
    elif hasattr(kokkos, 'HipSpace'):
        mem_space = kokkos.HipSpace
        layout = kokkos.LayoutLeft

    H_SIZE = 5.0

    if rank == 0:
        print(f"--- FUnTiDES Benchmark ({args.type_scaling}) ---", flush=True)
        print(f"Ranks (GPUs):    {size}", flush=True)
        print(f"Target Elements: {args.ex} x {args.ey} x {args.ez}", flush=True)
        if not args.benchmark:
            os.makedirs("snapshots", exist_ok=True)

    global_params = Params()
    global_params.order = 2
    global_params.ex = args.ex
    global_params.ey = args.ey
    global_params.ez = args.ez
    global_params.lx = args.ex * H_SIZE
    global_params.ly = args.ey * H_SIZE
    global_params.lz = args.ez * H_SIZE
    global_params.is_model_on_nodes = False
    global_params.is_elastic = False

    partitioner = Partitioner()
    local_params = partitioner.partition(global_params, rank, size)

    # Verify Partition
    local_ex = local_params.ex
    partition_info = comm.gather(local_ex, root=0)
    if rank == 0:
        total_ex = sum(partition_info)
        if total_ex != args.ex:
            print(f"ERROR: Partition mismatch {total_ex} != {args.ex}", flush=True)
            comm.Abort(1)

    builder = Builder(local_params.ex, local_params.lx, local_params.ey, local_params.ly, local_params.ez, local_params.lz, False, False, local_params.origin_x, local_params.origin_y, local_params.origin_z)
    mesh = builder.get_model()

    solver = Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.STRUCT, Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ACOUSTIC, 2)

    solver.compute_fe_init(mesh)
    mass_matrix = solver.get_mass_matrix()

    nx = int(local_params.ex * 2 + 1)
    ny = int(local_params.ey * 2 + 1)
    nz = int(local_params.ez * 2 + 1)
    n_nodes_local = int(mesh.get_number_of_nodes())
    n_elem = int(mesh.get_number_of_elements())

    exchanger = GhostExchanger(nx, ny, nz)
    exchanger.sync_accumulate(mass_matrix)

    kk_prev = kokkos.array([n_nodes_local], dtype=kokkos.float32, space=mem_space, layout=layout)
    kk_curr = kokkos.array([n_nodes_local], dtype=kokkos.float32, space=mem_space, layout=layout)
    np.array(kk_prev, copy=False)[:] = 0.0
    np.array(kk_curr, copy=False)[:] = 0.0

    wavefield = Solver.WavefieldAcoustic(kk_prev, kk_curr)

    n_rhs = 1
    kk_rhs_t = kokkos.array([n_rhs, args.steps], dtype=kokkos.float32, space=mem_space, layout=layout)
    kk_rhs_e = kokkos.array([n_rhs], dtype=kokkos.int32, space=mem_space, layout=layout)
    kk_rhs_w = kokkos.array([n_rhs, 27], dtype=kokkos.float32, space=mem_space, layout=layout)

    t_axis = np.linspace(0, 1.0, args.steps)
    f0 = 5.0
    t_peak = 0.15 + (rank * 0.05)
    src_val = (1 - 2*(np.pi*f0*(t_axis-t_peak))**2) * np.exp(-(np.pi*f0*(t_axis-t_peak))**2)
    np_rhs_t = np.array(kk_rhs_t, copy=False)
    np_rhs_t[0, :] = src_val.astype(np.float32)
    lex = local_params.ex
    ley = local_params.ey
    lez = local_params.ez
    center_idx = (lex // 2) + (ley // 2) * lex + (lez // 2) * lex * ley
    np_rhs_e = np.array(kk_rhs_e, copy=False)
    np_rhs_e[0] = center_idx
    np_rhs_w = np.array(kk_rhs_w, copy=False)
    np_rhs_w[0, :] = 1.0 / 27.0

    rhs = Solver.RhsAcoustic(kk_rhs_t, kk_rhs_e, kk_rhs_w)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    # Sync before loop
    comm.Barrier()

    # --- END INIT TIMER ---
    t_init_duration = time.time() - t_init_start

    # --- Loop ---
    dt = 0.0005
    total_compute_time = 0.0
    total_io_time = 0.0

    if rank == 0:
        print(f"Initialization done in {t_init_duration:.4f}s. Starting loop...", flush=True)

    for t in range(args.steps):

        # === COMPUTE TIMER START ===
        t_c_start = time.time()

        solver.compute_forces(dt, t, data)
        force_vec = solver.get_force_vector(0)
        exchanger.sync_accumulate(force_vec)
        solver.update_solution(dt, data)

        # MPI sync acts as fence
        total_compute_time += (time.time() - t_c_start)
        # === COMPUTE TIMER END ===

        if not args.benchmark and t % 100 == 0:
            if rank == 0: print(f"Step {t}", flush=True)
            active_field = np.array(kk_curr, copy=False)
            # IO function has its own internal timer
            t_io = save_snapshot_distributed(active_field, nx, ny, nz, t, "snapshots")
            total_io_time += t_io
        elif t % 50 == 0 and rank == 0:
            print(".", end="", flush=True)

        data.swap_wavefields()

    comm.Barrier()

    # --- Reporting & CSV ---
    if rank == 0:
        total_wall_time = t_init_duration + total_compute_time + total_io_time

        print("\n" + "="*40, flush=True)
        print(f"RESULTS ({args.type_scaling} scaling)", flush=True)
        print("="*40, flush=True)
        print(f"Init Time      : {t_init_duration:.4f} s", flush=True)
        print(f"Compute Time   : {total_compute_time:.4f} s", flush=True)
        print(f"I/O Plot Time  : {total_io_time:.4f} s", flush=True)
        print(f"Total Wall     : {total_wall_time:.4f} s", flush=True)
        print("="*40, flush=True)

        # Write to CSV
        csv_file = f"results_{args.type_scaling}.csv"
        file_exists = os.path.isfile(csv_file)

        with open(csv_file, "a") as f:
            if not file_exists:
                f.write("gpus,global_ex,ey,ez,steps,t_init,t_compute,t_io,t_total\n")
            f.write(f"{size},{args.ex},{args.ey},{args.ez},{args.steps},"
                    f"{t_init_duration:.5f},{total_compute_time:.5f},{total_io_time:.5f},{total_wall_time:.5f}\n")

        print(f"Data appended to {csv_file}", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ex", type=int, default=400, help="Global Elements X")
    parser.add_argument("--ey", type=int, default=300, help="Global Elements Y")
    parser.add_argument("--ez", type=int, default=300, help="Global Elements Z")
    parser.add_argument("--steps", type=int, default=1500, help="Time steps")
    parser.add_argument("--benchmark", action="store_true", help="Disable plotting")
    parser.add_argument("--type-scaling", type=str, default="custom", help="Label for CSV (weak/strong)")
    args = parser.parse_args()

    kokkos.initialize()
    try:
        run_simulation(args)
    finally:
        kokkos.finalize()

if __name__ == "__main__":
    main()
