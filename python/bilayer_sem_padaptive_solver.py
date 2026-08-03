#!/usr/bin/env python3
"""
Scenario to solve: (acoustic case)
----------------------------------
    - water layer at the top (Vp=1500, needs high order) 
    - "rock" layer at the bottom (Vp=4000, needs low order)
    - one horizontal interface between the two layers
Target architecture: SEM_high(water) -> DG_high -> DG_low -> SEM_low(rock).

Sources/receivers:
 One source is placed in the SEM_low(rock) layer. On receiver is placed in the SEM_high(water) layer.
 Source and reveiver can be moved as long as they stay in their respective SEM domain

Implementation 
--------------
Solvers:
 Three independent solver instances, each on its OWN minimally-sized local mesh:
    dgsem_top    = DGSEMsolver       (SEM_high bulk + thin DG_high sub-region)
    dg_pad_mid   = DGPAdaptiveSolver (pMax=DG_high sub-region + pMin=DG_low sub-region)
    dgsem_bottom = DGSEMsolver       (thin DG_low sub-region + SEM_low bulk)

Meshes:
 Three overlapping meshes. Each local mesh carries a 1-element-deep GHOST layer 
 at every interface between meshes, mirroring the neighboring solver's boundary 
 elements.

Time step loop:
  - Each solver calls compute_one_step()
  - Each solver's ghost layer are overwritten with the value from its neighboring solver
  - Call swap_wavefields() on each data set

"""

import os
import sys
import time
from pathlib import Path

import numpy as np
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))

os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OMP_THREAD_LIMIT", "6")
os.environ.setdefault("KOKKOS_NUM_THREADS", "6")

import kokkos  # noqa: E402  (from pykokkos-base, built with the TPLs)
import cupy as cp  # noqa: E402  (device-side ghost exchange / trace gather, see as_cupy)

from pyfuntides import model, solver  # noqa: E402

from bilayer_mesh_common import (  # noqa: E402
    build_bilayer_model, compute_rhs_weights, detect_default_memspace, kk_zeros,
    dense_grid_rows, dg_local_block_getter, sem_local_block_getter, as_cupy, PhaseTimer,
)

# =============================================================================
# Parameters
# =============================================================================
EX, EY, EZ = 100, 45, 60                # elements per direction
LX, LY, LZ = 4000.0, 1800.0, 1500.0     # domain size [m], per validation scenario doc
WATER_ROCK_ZBOUNDARY = 1200.0   # water_rock interface in the global mesh (top mesh + middle mesh + bottom mesh)
                                 # thin water (300m/20%) over thick rock (1200m/80%) -- realistic
                                 # marine-seismic proportion, matches bilayer_uniform_solver.py /
                                 # bilayer_dg_padaptive_solver.py
VEL_WATER = 1500.0
VEL_ROCK = 4000.0

DT = 0.0001
N_SAMPLES = 10000
F0 = 10.0
PRINT_INTERVAL = 100     # stdout |p|_max diagnostics every N steps
SNAP_INTERVAL = 100      # x-z slice snapshot every N steps (gnuplot: plot 'file' matrix with image)
SYNC_TIMERS = True       # fence the device before stopping each phase's clock (see PhaseTimer);
                         # set False to check the fences aren't inflating the reported total
SRC_COORD = (2000.0, 900.0, 1450.0)     # water region, matches bilayer_uniform_solver.py
RCV_Y = 900.0
RCV_Z_TOP = 1400.0                              # water region (top mesh)
RCV_Z_BOT = 200.0                               # rock region (bottom mesh)
RCV_X_COORDS = np.linspace(200.0, 3800.0, 7)    # same X line reused for both domains
N_SRC = 1
N_RCV = len(RCV_X_COORDS)


ORDER_MIN = 2   # fast region, bottom(rock)
ORDER_MAX = 6   # slow region, top(water) -- hack-padaptive26.patch only builds the (2,6) pair
N_DOF_MIN = (ORDER_MIN + 1) ** 3
N_DOF_MAX = (ORDER_MAX + 1) ** 3

# CartesianStructBuilder is templated on order; the mesh's order must match its solver's
# order exactly (DGSEMsolver::computeFEInit does a dynamic_cast<MESH_TYPE*>, which fails
# silently-as-a-RuntimeError if the mesh was built at a different order than the solver).
CartesianStructBuilderTop = getattr(model, f"CartesianStructBuilder_f32_i32_O{ORDER_MAX}")
CartesianStructBuilderBot = getattr(model, f"CartesianStructBuilder_f32_i32_O{ORDER_MIN}")
N_ELEMENTS = EX * EY * EZ

ELEM_SIZE_Z = LZ / EZ
N_ELEM_TO_BOUNDARY = WATER_ROCK_ZBOUNDARY / ELEM_SIZE_Z
assert abs(N_ELEM_TO_BOUNDARY - round(N_ELEM_TO_BOUNDARY)) < 1e-6, (
    f"WATER_ROCK_ZBOUNDARY={WATER_ROCK_ZBOUNDARY} is not aligned to an element boundary.")

# Top Mesh (water)
LZ_top = LZ - WATER_ROCK_ZBOUNDARY  # z axe domain size [m]
EZ_top = int(EZ * LZ_top/LZ)        # elements in z direction
WATER_DGSEM_ZBOUNDARY = 2 * LZ/EZ   # DG-SEM boundary: 2 z-directed element in DG (1 ghost + 1 truly solved)
N_ELEMENTS_top = EX * EY * EZ_top

# Middle Mesh (water-rock interface)
EZ_mid = 4                # elements in z direction (2 elements truly solved + 1 ghost on each top and bottom interfaces)
LZ_mid = EZ_mid * LZ/EZ   # z axe domain size [m]
Z_ELEM_INTERFACE = 2      # element index of (p-adaptive) interface: used to build the bilayer model
PADAPTIVE_ZBOUNDARY = Z_ELEM_INTERFACE * LZ/EZ # p-adaptive boundary: middle of the local mesh
N_ELEMENTS_mid = EX * EY * EZ_mid

# Bottom Mesh (rock) - this domain is symmetrized along the plan z=1000 to adapt to the DG-SEM solver
LZ_bot = WATER_ROCK_ZBOUNDARY     # z axe domain size [m]
EZ_bot = int(EZ * LZ_bot/LZ)      # elements in z direction
ROCK_DGSEM_ZBOUNDARY = 2 * LZ/EZ  # DG-SEM boundary: 2 z-directed element in DG (1 ghost + 1 truly solved)
N_ELEMENTS_bot = EX * EY * EZ_bot

assert (SRC_COORD[2] > WATER_ROCK_ZBOUNDARY + WATER_DGSEM_ZBOUNDARY or
        SRC_COORD[2] < WATER_ROCK_ZBOUNDARY - ROCK_DGSEM_ZBOUNDARY), \
    f"SRC_COORD z={SRC_COORD[2]} falls in the DG/ghost cap, not in a SEM domain"

for _name, _z in (("RCV_Z_TOP", RCV_Z_TOP), ("RCV_Z_BOT", RCV_Z_BOT)):
    assert (_z > WATER_ROCK_ZBOUNDARY + WATER_DGSEM_ZBOUNDARY or
            _z < WATER_ROCK_ZBOUNDARY - ROCK_DGSEM_ZBOUNDARY), \
        f"{_name}={_z} falls in the DG/ghost cap, not in a SEM domain"



def write_uniform_model_file(n_elem, vp):
    fd, path = tempfile.mkstemp(suffix=".txt", text=True)
    with os.fdopen(fd, "w") as f:
        f.write("Model Vp element\n")
        f.write(f"{n_elem}\n")
        for _ in range(n_elem):
            f.write(f"{vp}\n")
    return path


def face_dofs(order, k_loc):
    return np.array([i + j*(order+1) + k_loc*(order+1)**2
                      for j in range(order+1) for i in range(order+1)])


# -----------------------------------------------------------------------------
# Snapshot grid: ONE shared dense (x, z) grid for the whole global domain, at
# ORDER_MAX resolution. Every region (bot SEM/DG, mid pMin/pMax, top DG/SEM) is
# resampled onto it via true tensor-product GLL Lagrange reconstruction (see
# dense_grid_rows in bilayer_mesh_common.py) instead of each keeping its own
# native per-order z-row spacing -- that native spacing is what produced the
# "strata" look in gnuplot's matrix/image (row density jumps at every order
# change: SEM/DG order2 vs order6, or DG cap vs SEM bulk). One shared grid means
# every row is dz = LZ/(NZ_DENSE-1) apart everywhere, so the image reads like a
# single uniform-order SEM solution.
#
# Regions are described by (order, ex, dxe, dze, z0, sign, ez_count, ez_offset):
# global z = z0 + sign * local_z, ez_offset = how many local elements to skip
# before this region's own [0, ez_count) (matches the ranges the previous
# per-region dg_region_rows/sem_region_rows calls used).
# -----------------------------------------------------------------------------
NZ_DENSE = ORDER_MAX * EZ + 1
NX_DENSE = ORDER_MAX * EX + 1
DENSE_X = np.linspace(0.0, LX, NX_DENSE)
DENSE_Z = np.linspace(0.0, LZ, NZ_DENSE)
DXE = LX / EX
DZE = LZ / EZ  # == LZ_top/EZ_top == LZ_mid/EZ_mid == LZ_bot/EZ_bot (25 m here)

# Per-mesh (z0, sign) mapping local z=0 (element index 0, BEFORE any ez_offset) to global z:
# global z = z0_mesh + sign * local_z. bot's local z runs opposite to global (see the
# LZ_bot - z conversions used for its source/receivers above); mid and top run the same way.
BOT_Z0, BOT_SIGN = WATER_ROCK_ZBOUNDARY, -1
MID_Z0 = WATER_ROCK_ZBOUNDARY - Z_ELEM_INTERFACE * DZE  # global z at mid mesh's local z=0
MID_SIGN = 1
TOP_Z0, TOP_SIGN = WATER_ROCK_ZBOUNDARY, 1


def region_z0(mesh_z0, sign, ez_offset):
    """z0 for dense_grid_rows when the region starts `ez_offset` elements into the mesh
    (dense_grid_rows/get_local_block's ez_i is 0-based *within the region*, not the mesh,
    so the region's own z0 must be shifted by the skipped elements' physical extent).
    """
    return mesh_z0 + sign * ez_offset * DZE


def main():
    kokkos.initialize()
    try:
        run()
    finally:
        kokkos.finalize()


def run():
    t_start = time.perf_counter()
    timer = PhaseTimer(sync=SYNC_TIMERS)
    timer.tic()
    memspace, layout = detect_default_memspace()

    # -------------------------------------------------------------------
    # Models 
    # -------------------------------------------------------------------
    vp_path = write_uniform_model_file(N_ELEMENTS_top, VEL_WATER)
    try:
        model_top = CartesianStructBuilderTop(
            EX, LX, EY, LY, EZ_top, LZ_top,
            is_model_on_nodes=False, is_elastic=False,
            # free_surface_on_top=False: the uniform/reference and pure-DG-p-adaptive scripts
            # build their model via build_bilayer_model() with an EMPTY boundaries_t array ->
            # isFreeSurface() false everywhere -> absorbing top. The DG kernels never consult
            # isFreeSurface at all, so an all-absorbing setup is the only boundary condition
            # every config of the comparison can share. With True here, this run alone had a
            # reflecting surface -> strong source/receiver ghost (delayed inflated 2nd lobe +
            # trailing dip) absent from the reference, wrongly read as a coupling artifact.
            model_file=vp_path).get_model(free_surface_on_top=False)
    finally:
        os.remove(vp_path)

    model_mid, nodes_x, nodes_y, nodes_z, _model_keepalive = build_bilayer_model(
        ORDER_MAX, EX, EY, EZ_mid, LX, LY, LZ_mid,
        VEL_WATER, VEL_ROCK, Z_ELEM_INTERFACE, 
        memspace, layout)
    
    vp_path = write_uniform_model_file(N_ELEMENTS_bot, VEL_ROCK)
    try:
        model_bot = CartesianStructBuilderBot(
            EX, LX, EY, LY, EZ_bot, LZ_bot, 
            is_model_on_nodes=False, is_elastic=False, 
            model_file=vp_path).get_model(free_surface_on_top=False)
    finally:
        os.remove(vp_path)

    print(f"Top model built: {N_ELEMENTS_top} elements, {model_top.get_number_of_nodes()} nodes")
    print(f"Middle model built: {N_ELEMENTS_mid} elements, {model_mid.get_number_of_nodes()} nodes")
    print(f"Bottom model built: {N_ELEMENTS_bot} elements, {model_bot.get_number_of_nodes()} nodes")


    # -------------------------------------------------------------------
    # Solvers
    # -------------------------------------------------------------------
    solver_top = solver.create_solver(
        method_type=solver.MethodType.DGSEM,
        implem_type=solver.ImplemType.MAKUTU,
        mesh_type=solver.MeshType.STRUCT,
        model_location=solver.ModelLocationType.ONELEMENTS,
        physic_type=solver.PhysicType.ACOUSTIC,
        order=ORDER_MAX
    )
    solver_top.set_z_boundary(WATER_DGSEM_ZBOUNDARY)
    solver_mid = solver.create_solver(
        method_type=solver.MethodType.DGPADAPTIVE,
        implem_type=solver.ImplemType.MAKUTU,
        mesh_type=solver.MeshType.UNSTRUCT,
        model_location=solver.ModelLocationType.ONELEMENTS,
        physic_type=solver.PhysicType.ACOUSTIC,
        order=ORDER_MAX,
        order_min=ORDER_MIN
    )
    solver_mid.set_z_boundary(PADAPTIVE_ZBOUNDARY)
    solver_bot = solver.create_solver(
        method_type=solver.MethodType.DGSEM,
        implem_type=solver.ImplemType.MAKUTU,
        mesh_type=solver.MeshType.STRUCT,
        model_location=solver.ModelLocationType.ONELEMENTS,
        physic_type=solver.PhysicType.ACOUSTIC,
        order=ORDER_MIN
    )
    solver_bot.set_z_boundary(ROCK_DGSEM_ZBOUNDARY)
    print("Solver created")

    solver_top.compute_fe_init(model_top)
    solver_mid.compute_fe_init(model_mid)
    solver_bot.compute_fe_init(model_bot)
    print("Solver initialized")


    # -------------------------------------------------------------------
    # Wavefield
    # -------------------------------------------------------------------
    # TOP wavefield
    n_node_top = model_top.get_number_of_nodes()
    pn_top_sem_prev = kk_zeros((n_node_top,), kokkos.float32, memspace, layout)
    pn_top_sem_curr = kk_zeros((n_node_top,), kokkos.float32, memspace, layout)
    pn_top_dg_prev = kk_zeros((N_ELEMENTS_top, N_DOF_MAX), kokkos.float32, memspace, layout)
    pn_top_dg_curr = kk_zeros((N_ELEMENTS_top, N_DOF_MAX), kokkos.float32, memspace, layout)

    # MIDDLE wavefield
    pn_mid_pmax_prev = kk_zeros((N_ELEMENTS_mid, N_DOF_MAX), kokkos.float32, memspace, layout)
    pn_mid_pmax_curr = kk_zeros((N_ELEMENTS_mid, N_DOF_MAX), kokkos.float32, memspace, layout)
    pn_mid_pmin_prev = kk_zeros((N_ELEMENTS_mid, N_DOF_MIN), kokkos.float32, memspace, layout)
    pn_mid_pmin_curr = kk_zeros((N_ELEMENTS_mid, N_DOF_MIN), kokkos.float32, memspace, layout)

    # BOTTOM wavefield
    n_node_bot = model_bot.get_number_of_nodes()
    pn_bot_dg_prev = kk_zeros((N_ELEMENTS_bot, N_DOF_MIN), kokkos.float32, memspace, layout)
    pn_bot_dg_curr = kk_zeros((N_ELEMENTS_bot, N_DOF_MIN), kokkos.float32, memspace, layout)
    pn_bot_sem_prev = kk_zeros((n_node_bot,), kokkos.float32, memspace, layout)
    pn_bot_sem_curr = kk_zeros((n_node_bot,), kokkos.float32, memspace, layout)
    print("Wavefield buffers allocated ")


    # -------------------------------------------------------------------
    # Source: Ricker in the rock(SEM_low) domain 
    # port of SolverUtils::evaluateRicker order=2
    # -------------------------------------------------------------------
    t = np.arange(N_SAMPLES, dtype=np.float32) * DT
    TPEAK = 0.2
    lam = (np.pi * F0) ** 2
    tau = t - TPEAK
    ricker = 2 * lam * (2 * lam * tau ** 2 - 1) * np.exp(-lam * tau ** 2)
    ricker = np.where((t <= -0.9 * TPEAK) | (t >= 2.9 * TPEAK), 0.0, ricker)

    # TOP source (ricker or 0)
    rhs_top_dg_term = kk_zeros((N_SRC, N_SAMPLES), kokkos.float32, memspace, layout)
    rhs_top_sem_term = kk_zeros((N_SRC, N_SAMPLES), kokkos.float32, memspace, layout)
    rhs_top_element = kk_zeros((N_SRC,), kokkos.int32, memspace, layout)
    rhs_top_weights = kk_zeros((N_SRC, N_DOF_MAX), kokkos.float32, memspace, layout)

    # MIDDLE source (0)
    rhs_mid_pmin_term = kk_zeros((0, 0), kokkos.float32, memspace, layout)
    rhs_mid_pmax_term = kk_zeros((0, 0), kokkos.float32, memspace, layout)
    rhs_mid_element = kk_zeros((0,), kokkos.int32, memspace, layout)
    rhs_mid_pmin_weights = kk_zeros((0, 0), kokkos.float32, memspace, layout)
    rhs_mid_pmax_weights = kk_zeros((0, 0), kokkos.float32, memspace, layout)

    # BOTTOM source (ricker or 0)
    rhs_bot_dg_term = kk_zeros((N_SRC, N_SAMPLES), kokkos.float32, memspace, layout)
    rhs_bot_sem_term = kk_zeros((N_SRC, N_SAMPLES), kokkos.float32, memspace, layout)
    rhs_bot_element = kk_zeros((N_SRC,), kokkos.int32, memspace, layout)
    rhs_bot_weights = kk_zeros((N_SRC, N_DOF_MIN), kokkos.float32, memspace, layout)


    if SRC_COORD[2] > WATER_ROCK_ZBOUNDARY:
        src_term, src_element, src_weights = rhs_top_sem_term, rhs_top_element, rhs_top_weights
        local_src_z, dz_src = SRC_COORD[2] - WATER_ROCK_ZBOUNDARY, LZ_top / EZ_top
        src_order = ORDER_MAX
    else:
        # NOTE: the bottom mesh's local z runs REVERSED vs global z; the zeta coordinate
        # fed to compute_rhs_weights below is expressed in that reversed local frame, which
        # mirrors the weights along z. Fine while the source stays in the top/water domain
        # (current scenario); revisit before placing the source in the bottom mesh.
        src_term, src_element, src_weights = rhs_bot_sem_term, rhs_bot_element, rhs_bot_weights
        local_src_z, dz_src = LZ_bot - SRC_COORD[2], LZ_bot / EZ_bot
        src_order = ORDER_MIN

    # True tensorised GLL nodal-basis weights at the physical source position -- NOT a delta
    # on local dof 0: that hack pinned the source to the element's first node (up to one
    # element-size position error, e.g. 20 m in y here) and made the injected amplitude
    # depend on the order, which showed up as a time shift + amplitude mismatch vs the
    # uniform reference (same fix as bilayer_uniform_solver.py / compute_rhs_weights).
    dxe_src, dye_src = LX / EX, LY / EY
    ex_src = int(SRC_COORD[0] / dxe_src)
    ey_src = int(SRC_COORD[1] / dye_src)
    ez_src = int(local_src_z / dz_src)

    np.array(src_term, copy=False)[0, :] = ricker
    np.array(src_element, copy=False)[0] = ex_src + ey_src * EX + ez_src * EX * EY
    src_elem_origin = (ex_src * dxe_src, ey_src * dye_src, ez_src * dz_src)
    np.array(src_weights, copy=False)[0, :] = compute_rhs_weights(
        src_order, (SRC_COORD[0], SRC_COORD[1], local_src_z), src_elem_origin,
        dxe_src, dye_src, dz_src)


    # -------------------------------------------------------------------
    # Data strucure
    # -------------------------------------------------------------------
    # TOP data
    wavefield_top = solver.DGSEMWavefieldAcoustic(pn_top_dg_prev, pn_top_dg_curr, pn_top_sem_prev, pn_top_sem_curr)
    rhs_top = solver.DGSEMRhsAcoustic(rhs_top_dg_term, rhs_top_sem_term, rhs_top_element, rhs_top_weights)
    data_top = solver.DGSEMsolverData(wavefield_top, rhs_top)

    # MIDDLE data
    wavefield_mid = solver.DGPAdaptiveWavefieldAcoustic(pn_mid_pmin_prev, pn_mid_pmin_curr, pn_mid_pmax_prev, pn_mid_pmax_curr)
    rhs_mid = solver.DGPAdaptiveRhsAcoustic(rhs_mid_pmin_term, rhs_mid_pmax_term, rhs_mid_element, rhs_mid_pmin_weights, rhs_mid_pmax_weights)
    data_mid = solver.DGPAdaptiveSolverData(wavefield_mid, rhs_mid)

    # BOTTOM data
    wavefield_bot = solver.DGSEMWavefieldAcoustic(pn_bot_dg_prev, pn_bot_dg_curr, pn_bot_sem_prev, pn_bot_sem_curr)
    rhs_bot = solver.DGSEMRhsAcoustic(rhs_bot_dg_term, rhs_bot_sem_term, rhs_bot_element, rhs_bot_weights)
    data_bot = solver.DGSEMsolverData(wavefield_bot, rhs_bot)


    # -------------------------------------------------------------------
    # Receiver lines: nearest global SEM nodes along RCV_X_COORDS, one line
    # in each SEM domain (water/top, rock/bottom).
    # -------------------------------------------------------------------
    ny_nodes_bot = ORDER_MIN * EY + 1

    # TOP receiver line (water, SEM_high)
    rcv_local_z_top = RCV_Z_TOP - WATER_ROCK_ZBOUNDARY
    dz_top = LZ_top / (ORDER_MAX * EZ_top)
    nx_nodes = ORDER_MAX * EX + 1
    ny_nodes = ORDER_MAX * EY + 1
    dx_top = LX / (ORDER_MAX * EX)
    dy_top = LY / (ORDER_MAX * EY)
    j_top = round(RCV_Y / dy_top)
    k_top = round(rcv_local_z_top / dz_top)
    rcv_nodes_top = np.array([
        round(x / dx_top) + j_top * nx_nodes + k_top * nx_nodes * ny_nodes
        for x in RCV_X_COORDS
    ])

    # BOTTOM receiver line (rock, SEM_low) -- local z runs reversed vs global,
    # same convention as the source placement above.
    rcv_local_z_bot = LZ_bot - RCV_Z_BOT
    dz_bot = LZ_bot / (ORDER_MIN * EZ_bot)
    nx_nodes_bot = ORDER_MIN * EX + 1
    dx_bot = LX / (ORDER_MIN * EX)
    dy_bot = LY / (ORDER_MIN * EY)
    j_bot = round(RCV_Y / dy_bot)
    k_bot = round(rcv_local_z_bot / dz_bot)
    rcv_nodes_bot = np.array([
        round(x / dx_bot) + j_bot * nx_nodes_bot + k_bot * nx_nodes_bot * ny_nodes_bot
        for x in RCV_X_COORDS
    ])


    # -------------------------------------------------------------------
    # Time step loop
    # -------------------------------------------------------------------
    face_high_max = face_dofs(ORDER_MAX, ORDER_MAX)
    face_low_max  = face_dofs(ORDER_MAX, 0)
    face_high_min = face_dofs(ORDER_MIN, ORDER_MIN)
    face_low_min  = face_dofs(ORDER_MIN, 0)

    layer = np.arange(EX*EY)
    
    top_ghost = layer
    top_real  = layer + 1*EX*EY

    bot_ghost = layer
    bot_real  = layer + 1*EX*EY

    mid_pmin_ghost = layer
    mid_pmin_real  = layer + 1*EX*EY
    mid_pmax_real  = layer + 2*EX*EY
    mid_pmax_ghost = layer + 3*EX*EY

    # -------------------------------------------------------------------
    # Device-side index arrays + trace buffers: the whole time loop runs without
    # numpy touching any wavefield buffer (ghost exchange, receiver gather and
    # |p|_max diagnostics all execute as device kernels through cupy on zero-copy
    # wrappers of the UVM buffers). The wavefields' UVM pages therefore stay
    # device-resident all run long instead of ping-ponging host<->device at every
    # step (was ~17% of wall-clock, see the "gather" row in PhaseTimer's report).
    # -------------------------------------------------------------------
    cp_top_ghost, cp_top_real = cp.asarray(top_ghost), cp.asarray(top_real)
    cp_bot_ghost, cp_bot_real = cp.asarray(bot_ghost), cp.asarray(bot_real)
    cp_mid_pmin_ghost, cp_mid_pmin_real = cp.asarray(mid_pmin_ghost), cp.asarray(mid_pmin_real)
    cp_mid_pmax_ghost, cp_mid_pmax_real = cp.asarray(mid_pmax_ghost), cp.asarray(mid_pmax_real)
    cp_face_high_max, cp_face_low_max = cp.asarray(face_high_max), cp.asarray(face_low_max)
    cp_face_high_min, cp_face_low_min = cp.asarray(face_high_min), cp.asarray(face_low_min)
    cp_rcv_top, cp_rcv_bot = cp.asarray(rcv_nodes_top), cp.asarray(rcv_nodes_bot)

    # Traces accumulate device-side in kokkos-owned UVM buffers (read back once at the
    # end through numpy); cupy writes them through zero-copy wrappers.
    rcv_trace_top_kk = kk_zeros((N_RCV, N_SAMPLES), kokkos.float32, memspace, layout)
    rcv_trace_bot_kk = kk_zeros((N_RCV, N_SAMPLES), kokkos.float32, memspace, layout)
    rcv_trace_top_cp = as_cupy(rcv_trace_top_kk)
    rcv_trace_bot_cp = as_cupy(rcv_trace_bot_kk)

    def cp_pmax(view):
        return float(cp.abs(as_cupy(view)).max())

    timer.toc("setup")

    for it in range(N_SAMPLES):
        timer.tic()
        solver_top.compute_one_step(DT, it, data_top)
        solver_mid.compute_one_step(DT, it, data_mid)
        solver_bot.compute_one_step(DT, it, data_bot)
        timer.toc("compute")

        # Swap FIRST, then exchange ghosts on the post-swap current buffers (p^{n+1}): these
        # are exactly the buffers the next step's flux kernels read. Exchanging before the
        # swap (as done previously) landed the neighbor values in what became the *previous*
        # buffer, so the seam ran with a one-step lag -- an O(dt) transparency defect at both
        # mesh-to-mesh interfaces that showed up as a spurious partial reflection.
        data_top.swap_wavefields()
        data_mid.swap_wavefields()
        data_bot.swap_wavefields()
        timer.toc("swap")

        # Device-side ghost exchange: cupy fancy-indexing kernels on zero-copy wrappers of
        # the UVM buffers (launch cost only when SYNC_TIMERS=False; execution overlaps into
        # the next step's compute).
        top_dg = as_cupy(wavefield_top.get_dg_current_field(0))
        bot_dg = as_cupy(wavefield_bot.get_dg_current_field(0))
        mid_pmax = as_cupy(wavefield_mid.get_pmax_current_field(0))
        mid_pmin = as_cupy(wavefield_mid.get_pmin_current_field(0))
        top_dg[cp.ix_(cp_top_ghost, cp_face_high_max)] = mid_pmax[cp.ix_(cp_mid_pmax_real, cp_face_high_max)]
        mid_pmax[cp.ix_(cp_mid_pmax_ghost, cp_face_low_max)] = top_dg[cp.ix_(cp_top_real, cp_face_low_max)]
        bot_dg[cp.ix_(cp_bot_ghost, cp_face_high_min)] = mid_pmin[cp.ix_(cp_mid_pmin_real, cp_face_low_min)]
        mid_pmin[cp.ix_(cp_mid_pmin_ghost, cp_face_high_min)] = bot_dg[cp.ix_(cp_bot_real, cp_face_low_min)]
        timer.toc("exchange")

        # Receiver traces accumulate on the device; read back once after the loop.
        rcv_trace_top_cp[:, it] = as_cupy(wavefield_top.get_sem_previous_field(0))[cp_rcv_top]
        rcv_trace_bot_cp[:, it] = as_cupy(wavefield_bot.get_sem_previous_field(0))[cp_rcv_bot]
        # CuPy launches on the legacy default stream; Kokkos' default instance stream may
        # differ, so drain before the next compute_one_step reads the exchanged ghosts. This
        # sync is load-bearing for correctness (not just timing), so it stays regardless of
        # SYNC_TIMERS -- PhaseTimer's own fence right after is then a cheap no-op.
        cp.cuda.Stream.null.synchronize()
        timer.toc("gather")

        if it % PRINT_INTERVAL == 0:
            print(f"step {it:5d}  |p|_max top.dg={cp_pmax(wavefield_top.get_dg_current_field(0)):.3e}"
                  f"  top.sem={cp_pmax(wavefield_top.get_sem_previous_field(0)):.3e}"
                  f"  mid.pmin={cp_pmax(wavefield_mid.get_pmin_current_field(0)):.3e}"
                  f"  mid.pmax={cp_pmax(wavefield_mid.get_pmax_current_field(0)):.3e}"
                  f"  bot.dg={cp_pmax(wavefield_bot.get_dg_current_field(0)):.3e}"
                  f"  bot.sem={cp_pmax(wavefield_bot.get_sem_previous_field(0)):.3e}")
        timer.toc("diag")

        if it % SNAP_INTERVAL == 0:
            # Single x-z slice (mid-Y), bottom-to-top in GLOBAL z: bot's SEM+DG, mid's
            # pMin+pMax, top's DG+SEM, each resampled by true tensor-product GLL Lagrange
            # reconstruction onto the ONE shared dense (DENSE_X, DENSE_Z) grid (see
            # dense_grid_rows in bilayer_mesh_common.py) instead of each region's own
            # native per-order z-row spacing -- that's what used to render as "strata"
            # (row density jumps at every SEM/DG or order2/order6 boundary). Every row
            # below is dz = LZ/(NZ_DENSE-1) apart everywhere, no reversal needed: each
            # dense_grid_rows call already returns rows in ascending global z.
            # Only place inside the loop where numpy touches the wavefields (full host
            # read, amortized over SNAP_INTERVAL steps). The per-step stream sync above
            # already drained the exchange/gather kernels; compute_one_step fences the
            # Kokkos stream at its own end, so host reads are safe here.
            pn_top_dg_curr_np = np.array(wavefield_top.get_dg_current_field(0), copy=False)
            pn_bot_dg_curr_np = np.array(wavefield_bot.get_dg_current_field(0), copy=False)
            pn_mid_pmax_curr_np = np.array(wavefield_mid.get_pmax_current_field(0), copy=False)
            pn_mid_pmin_curr_np = np.array(wavefield_mid.get_pmin_current_field(0), copy=False)
            sem_top_prev_np = np.array(wavefield_top.get_sem_previous_field(0), copy=False)
            sem_bot_prev_np = np.array(wavefield_bot.get_sem_previous_field(0), copy=False)
            timer.toc("snap_read")

            # Each DG cap is 2 elements (1 ghost + 1 truly-solved, see WATER/ROCK_DGSEM_ZBOUNDARY
            # above): the ghost only mirrors the neighboring mesh's own real boundary element, so
            # drawing it too double-draws that physical slab (visible as an extra "layer" at every
            # bot<->mid and mid<->top seam). ez_offset=1/count=1 below keeps only the truly-solved
            # element; local index 0 (the ghost) is skipped in every DG-cap call.
            ey_mid = EY // 2
            bot_sem_rows = dense_grid_rows(
                sem_local_block_getter(sem_bot_prev_np, ORDER_MIN, EX, EY, ny_nodes_bot, ny_nodes_bot // 2, 2),
                ORDER_MIN, EX, DXE, DZE, region_z0(BOT_Z0, BOT_SIGN, 2), BOT_SIGN, EZ_bot - 2, DENSE_X, DENSE_Z)
            bot_dg_rows = dense_grid_rows(
                dg_local_block_getter(pn_bot_dg_curr_np, ORDER_MIN, EX, EY, ey_mid, 1),
                ORDER_MIN, EX, DXE, DZE, region_z0(BOT_Z0, BOT_SIGN, 1), BOT_SIGN, 1, DENSE_X, DENSE_Z)
            mid_pmin_rows = dense_grid_rows(
                dg_local_block_getter(pn_mid_pmin_curr_np, ORDER_MIN, EX, EY, ey_mid, 1),
                ORDER_MIN, EX, DXE, DZE, region_z0(MID_Z0, MID_SIGN, 1), MID_SIGN, 1, DENSE_X, DENSE_Z)
            mid_pmax_rows = dense_grid_rows(
                dg_local_block_getter(pn_mid_pmax_curr_np, ORDER_MAX, EX, EY, ey_mid, Z_ELEM_INTERFACE),
                ORDER_MAX, EX, DXE, DZE, region_z0(MID_Z0, MID_SIGN, Z_ELEM_INTERFACE), MID_SIGN, 1, DENSE_X, DENSE_Z)
            top_dg_rows = dense_grid_rows(
                dg_local_block_getter(pn_top_dg_curr_np, ORDER_MAX, EX, EY, ey_mid, 1),
                ORDER_MAX, EX, DXE, DZE, region_z0(TOP_Z0, TOP_SIGN, 1), TOP_SIGN, 1, DENSE_X, DENSE_Z)
            top_sem_rows = dense_grid_rows(
                sem_local_block_getter(sem_top_prev_np, ORDER_MAX, EX, EY, ny_nodes, ny_nodes // 2, 2),
                ORDER_MAX, EX, DXE, DZE, region_z0(TOP_Z0, TOP_SIGN, 2), TOP_SIGN, EZ_top - 2, DENSE_X, DENSE_Z)

            np.savetxt(f"slice_{it:05d}.dat", np.array(
                bot_sem_rows + bot_dg_rows + mid_pmin_rows + mid_pmax_rows + top_dg_rows + top_sem_rows))
            timer.toc("snap_write")

    # Drain the cupy stream before reading the accumulated traces from the host
    # (the Kokkos stream is already fenced by the last compute_one_step).
    cp.cuda.Stream.null.synchronize()
    rcv_trace_top = np.array(rcv_trace_top_kk, copy=False)
    rcv_trace_bot = np.array(rcv_trace_bot_kk, copy=False)

    np.savetxt("bilayer_sem_padaptive_receiver_trace_top.txt",
               np.column_stack([t] + [rcv_trace_top[r] for r in range(N_RCV)]),
               header="time " + " ".join(f"pressure_rcv{r}" for r in range(N_RCV)) +
                      " (bilayer p-adaptive, receiver line, water/top mesh)")
    np.savetxt("bilayer_sem_padaptive_receiver_trace_bot.txt",
               np.column_stack([t] + [rcv_trace_bot[r] for r in range(N_RCV)]),
               header="time " + " ".join(f"pressure_rcv{r}" for r in range(N_RCV)) +
                      " (bilayer p-adaptive, receiver line, rock/bottom mesh)")
    print("Wrote bilayer_sem_padaptive_receiver_trace_{top,bot}.txt")
    timer.toc("io_final")

    # -------------------------------------------------------------------
    # Cost report: wall-clock time and total DOF count (approximate, local
    # per-element DOF summed per region -- consistent across runs, not a
    # unique-global-node count).
    # -------------------------------------------------------------------
    n_elem_mid_pmin = EX * EY * Z_ELEM_INTERFACE
    n_elem_mid_pmax = EX * EY * (EZ_mid - Z_ELEM_INTERFACE)
    total_dof = (N_DOF_MAX * N_ELEMENTS_top
                 + N_DOF_MIN * n_elem_mid_pmin + N_DOF_MAX * n_elem_mid_pmax
                 + N_DOF_MIN * N_ELEMENTS_bot)
    elapsed = time.perf_counter() - t_start
    timer.report(elapsed, extra=f"dof={total_dof}  steps/s={N_SAMPLES / elapsed:.2f}")


if __name__ == "__main__":
    main()