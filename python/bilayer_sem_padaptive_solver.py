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
from numpy.polynomial import legendre
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))

os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OMP_THREAD_LIMIT", "6")
os.environ.setdefault("KOKKOS_NUM_THREADS", "6")

import kokkos  # noqa: E402  (from pykokkos-base, built with the TPLs)

from pyfuntides import model, solver  # noqa: E402

from bilayer_mesh_common import (  # noqa: E402
    detect_default_memspace, kk_zeros, build_bilayer_model,
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


_GLL_CACHE = {}


def _gll(order):
    """Gauss-Lobatto-Legendre nodes (on [-1, 1]) and barycentric interpolation weights for a
    degree-`order` element -- matches the solver's internal 1D nodal basis (endpoints plus the
    interior roots of P'_order). Cached per order since called every row/every snapshot.
    """
    if order in _GLL_CACHE:
        return _GLL_CACHE[order]
    if order == 1:
        nodes = np.array([-1.0, 1.0])
    else:
        interior = legendre.Legendre.basis(order).deriv().roots().real
        nodes = np.sort(np.concatenate(([-1.0], interior, [1.0])))
    weights = np.array([1.0 / np.prod(nodes[j] - np.delete(nodes, j)) for j in range(len(nodes))])
    _GLL_CACHE[order] = (nodes, weights)
    return nodes, weights


def _lagrange_interp_1d(nodes, weights, values, x):
    """Barycentric Lagrange interpolation of `values` (sampled at `nodes`) at point `x`."""
    diff = x - nodes
    hit = np.flatnonzero(np.abs(diff) < 1e-12)
    if hit.size:
        return values[hit[0]]
    terms = weights / diff
    return float(np.dot(terms, values) / np.sum(terms))


def _dense_x_frac(ix, ex, dense_order):
    """Physical x position of dense-grid index `ix`, in element units (integer part = element
    index, fractional part = position within that element), plus its reference coordinate
    in [-1, 1]. Clamps the last point to the right edge of the last element.
    """
    if ix == dense_order * ex:
        return ex - 1, 1.0
    frac = ix / dense_order
    ix_e = int(frac)
    return ix_e, 2.0 * (frac - ix_e) - 1.0


def dg_region_rows(field_curr, order, ex, ey, ez_start, ez_count, dense_order):
    """X-z rows (mid-Y, one row per (z-element, z-dof)) of a per-element DG/DGPAdaptive field
    region, resampled onto the dense_order x-grid. When order != dense_order, reconstructs the
    true polynomial via Lagrange interpolation on the region's own GLL nodes (not
    nearest-neighbor) so the transition between p-adaptive regions renders smoothly.
    """
    n1d = order + 1
    ey_mid, ib_mid = ey // 2, order // 2
    nx_dense = dense_order * ex + 1
    identity = order == dense_order
    nodes, weights = _gll(order)
    rows = []
    for ez_i in range(ez_start, ez_start + ez_count):
        for iz_dof in range(n1d):
            row = np.empty(nx_dense, dtype=field_curr.dtype)
            dof0 = ib_mid * n1d + iz_dof * n1d * n1d
            for ix in range(nx_dense):
                ix_e, xi = _dense_x_frac(ix, ex, dense_order)
                elem = ix_e + ey_mid * ex + ez_i * ex * ey
                if identity:
                    ix_d = order if ix == nx_dense - 1 else ix % dense_order
                    row[ix] = field_curr[elem, ix_d + dof0]
                else:
                    local_vals = field_curr[elem, dof0:dof0 + n1d]
                    row[ix] = _lagrange_interp_1d(nodes, weights, local_vals, xi)
            rows.append(row)
    return rows


def sem_region_rows(sem_curr, order, ex, ey, ny, iz_start, iz_end, dense_order):
    """X-z rows (mid-Y, one row per z-node) of a nodal SEM field region, resampled onto the
    dense_order x-grid via true Lagrange interpolation when order != dense_order (see
    dg_region_rows), so the seam with the neighboring DG region lines up smoothly.
    """
    nx_native = order * ex + 1
    nx_dense = dense_order * ex + 1
    iy_mid = ny // 2
    identity = order == dense_order
    nodes, weights = _gll(order)
    rows = []
    for iz in range(iz_start, iz_end):
        row = np.empty(nx_dense, dtype=sem_curr.dtype)
        for ix in range(nx_dense):
            ix_e, xi = _dense_x_frac(ix, ex, dense_order)
            if identity:
                native_ix = ix_e * order + (order if ix == nx_dense - 1 else ix % dense_order)
                row[ix] = sem_curr[native_ix + iy_mid * nx_native + iz * nx_native * ny]
            else:
                base = ix_e * order
                local_vals = np.array([sem_curr[base + k + iy_mid * nx_native + iz * nx_native * ny]
                                        for k in range(order + 1)])
                row[ix] = _lagrange_interp_1d(nodes, weights, local_vals, xi)
        rows.append(row)
    return rows


def main():
    kokkos.initialize()
    try:
        run()
    finally:
        kokkos.finalize()


def run():
    t_start = time.perf_counter()
    memspace, layout = detect_default_memspace()

    # -------------------------------------------------------------------
    # Models 
    # -------------------------------------------------------------------
    vp_path = write_uniform_model_file(N_ELEMENTS_top, VEL_WATER)
    try:
        model_top = CartesianStructBuilderTop(
            EX, LX, EY, LY, EZ_top, LZ_top, 
            is_model_on_nodes=False, is_elastic=False, 
            model_file=vp_path).get_model(free_surface_on_top=True)
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
    else:
        src_term, src_element, src_weights = rhs_bot_sem_term, rhs_bot_element, rhs_bot_weights
        local_src_z, dz_src = LZ_bot - SRC_COORD[2], LZ_bot / EZ_bot

    np.array(src_term, copy=False)[0, :] = ricker
    np.array(src_element, copy=False)[0] = (
        int(SRC_COORD[0] / (LX / EX)) + 
        int(SRC_COORD[1] / (LY / EY)) * EX +
        int(local_src_z / dz_src) * EX * EY
    )
    np.array(src_weights, copy=False)[0, 0] = 1.0


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
    rcv_trace_top = np.zeros((N_RCV, N_SAMPLES), dtype=np.float32)

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
    rcv_trace_bot = np.zeros((N_RCV, N_SAMPLES), dtype=np.float32)


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

    t_setup = time.perf_counter() - t_start
    t_compute = 0.0
    t_sync = 0.0
    t_other = 0.0

    for it in range(N_SAMPLES):
        _t0 = time.perf_counter()
        solver_top.compute_one_step(DT, it, data_top)
        solver_mid.compute_one_step(DT, it, data_mid)
        solver_bot.compute_one_step(DT, it, data_bot)
        _t1 = time.perf_counter()
        t_compute += _t1 - _t0

        # Swap FIRST, then exchange ghosts on the post-swap current buffers (p^{n+1}): these
        # are exactly the buffers the next step's flux kernels read. Exchanging before the
        # swap (as done previously) landed the neighbor values in what became the *previous*
        # buffer, so the seam ran with a one-step lag -- an O(dt) transparency defect at both
        # mesh-to-mesh interfaces that showed up as a spurious partial reflection.
        data_top.swap_wavefields()
        data_mid.swap_wavefields()
        data_bot.swap_wavefields()

        pn_top_dg_curr_np = np.array(wavefield_top.get_dg_current_field(0), copy=False)
        pn_bot_dg_curr_np = np.array(wavefield_bot.get_dg_current_field(0), copy=False)
        pn_mid_pmax_curr_np = np.array(wavefield_mid.get_pmax_current_field(0), copy=False)
        pn_mid_pmin_curr_np = np.array(wavefield_mid.get_pmin_current_field(0), copy=False)

        pn_top_dg_curr_np[np.ix_(top_ghost, face_high_max)]      = pn_mid_pmax_curr_np[np.ix_(mid_pmax_real, face_high_max)]
        pn_mid_pmax_curr_np[np.ix_(mid_pmax_ghost, face_low_max)] = pn_top_dg_curr_np[np.ix_(top_real, face_low_max)]

        pn_bot_dg_curr_np[np.ix_(bot_ghost, face_high_min)]       = pn_mid_pmin_curr_np[np.ix_(mid_pmin_real, face_low_min)]
        pn_mid_pmin_curr_np[np.ix_(mid_pmin_ghost, face_high_min)] = pn_bot_dg_curr_np[np.ix_(bot_real, face_low_min)]
        _t2 = time.perf_counter()
        t_sync += _t2 - _t1

        sem_top_prev_np = np.array(wavefield_top.get_sem_previous_field(0), copy=False)
        sem_bot_prev_np = np.array(wavefield_bot.get_sem_previous_field(0), copy=False)
        rcv_trace_top[:, it] = sem_top_prev_np[rcv_nodes_top]
        rcv_trace_bot[:, it] = sem_bot_prev_np[rcv_nodes_bot]
        _t3 = time.perf_counter()
        t_other += _t3 - _t2

        if it % PRINT_INTERVAL == 0:
            print(f"step {it:5d}  |p|_max top.dg={np.abs(pn_top_dg_curr_np).max():.3e}"
                  f"  top.sem={np.abs(sem_top_prev_np).max():.3e}"
                  f"  mid.pmin={np.abs(pn_mid_pmin_curr_np).max():.3e}"
                  f"  mid.pmax={np.abs(pn_mid_pmax_curr_np).max():.3e}"
                  f"  bot.dg={np.abs(pn_bot_dg_curr_np).max():.3e}"
                  f"  bot.sem={np.abs(sem_bot_prev_np).max():.3e}")

        if it % SNAP_INTERVAL == 0:
            # Single x-z slice (mid-Y), bottom-to-top in GLOBAL z: bot's SEM+DG, mid's
            # pMin+pMax, top's DG+SEM -- all resampled onto the dense ORDER_MAX x-grid so
            # columns line up and the three regions concatenate into one valid matrix.
            # model_bot's local z runs opposite to global z (see build_layered_model /
            # the LZ_bot - z conversions above), so its rows are reversed here.
            bot_rows = (dg_region_rows(pn_bot_dg_curr_np, ORDER_MIN, EX, EY, 0, 2, ORDER_MAX)
                        + sem_region_rows(sem_bot_prev_np, ORDER_MIN, EX, EY, ny_nodes_bot,
                                          2 * ORDER_MIN, ORDER_MIN * EZ_bot + 1, ORDER_MAX))
            bot_rows = bot_rows[::-1]
            mid_rows = (dg_region_rows(pn_mid_pmin_curr_np, ORDER_MIN, EX, EY, 0, Z_ELEM_INTERFACE, ORDER_MAX)
                        + dg_region_rows(pn_mid_pmax_curr_np, ORDER_MAX, EX, EY, Z_ELEM_INTERFACE,
                                         EZ_mid - Z_ELEM_INTERFACE, ORDER_MAX))
            top_rows = (dg_region_rows(pn_top_dg_curr_np, ORDER_MAX, EX, EY, 0, 2, ORDER_MAX)
                        + sem_region_rows(sem_top_prev_np, ORDER_MAX, EX, EY, ny_nodes,
                                          2 * ORDER_MAX, ORDER_MAX * EZ_top + 1, ORDER_MAX))
            np.savetxt(f"slice_{it:05d}.dat", np.array(bot_rows + mid_rows + top_rows))


    np.savetxt("bilayer_sem_padaptive_receiver_trace_top.txt",
               np.column_stack([t] + [rcv_trace_top[r] for r in range(N_RCV)]),
               header="time " + " ".join(f"pressure_rcv{r}" for r in range(N_RCV)) +
                      " (bilayer p-adaptive, receiver line, water/top mesh)")
    np.savetxt("bilayer_sem_padaptive_receiver_trace_bot.txt",
               np.column_stack([t] + [rcv_trace_bot[r] for r in range(N_RCV)]),
               header="time " + " ".join(f"pressure_rcv{r}" for r in range(N_RCV)) +
                      " (bilayer p-adaptive, receiver line, rock/bottom mesh)")
    print("Wrote bilayer_sem_padaptive_receiver_trace_{top,bot}.txt")

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
    print(f"wall_clock={elapsed:.2f}s  total_dof={total_dof}"
          f"  t_setup={t_setup:.2f}s  t_compute={t_compute:.2f}s"
          f"  t_sync={t_sync:.2f}s  t_other={t_other:.2f}s")


if __name__ == "__main__":
    main()