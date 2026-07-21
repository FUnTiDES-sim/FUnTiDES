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
from pathlib import Path

import numpy as np
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))

os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OMP_THREAD_LIMIT", "6")
os.environ.setdefault("KOKKOS_NUM_THREADS", "6")

import kokkos  # noqa: E402  (from pykokkos-base, built with the TPLs)

from pyfuntides import model, solver  # noqa: E402

ModelUnstruct = model.ModelUnstruct_f32_i32
ModelUnstructData = model.ModelUnstructData_f32_i32

# =============================================================================
# Parameters
# =============================================================================
EX, EY, EZ = 100, 100, 100               # elements per direction
LX, LY, LZ = 2000.0, 2000.0, 2000.0   # domain size [m]
WATER_ROCK_ZBOUNDARY = 1000.0   # water_rock interface in the global mesh (top mesh + middle mesh + bottom mesh)
VEL_WATER = 2000.0
VEL_ROCK = 2000.0

DT = 0.0001
N_SAMPLES = 10000
F0 = 10.0
PRINT_INTERVAL = 100     # stdout |p|_max diagnostics every N steps
SNAP_INTERVAL = 100      # x-z slice snapshot every N steps (gnuplot: plot 'file' matrix with image)
SRC_COORD = (1000.0, 1000.0, 1400.0)   # rock region
RCV_COORD = (1500.0, 1000.0, 1700.0)   # water region
N_SRC = 1


ORDER_MIN = 2   # fast region, bottom(rock)
ORDER_MAX = 3   # slow region, top(water)
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

assert (RCV_COORD[2] > WATER_ROCK_ZBOUNDARY + WATER_DGSEM_ZBOUNDARY or
        RCV_COORD[2] < WATER_ROCK_ZBOUNDARY - ROCK_DGSEM_ZBOUNDARY), \
    f"RCV_COORD z={RCV_COORD[2]} falls in the DG/ghost cap, not in a SEM domain"



def detect_default_memspace():
    """Inspect the pybind11-generated docstring to know if this build is CPU or CUDA."""
    doc = solver.DGSEMWavefieldAcoustic.__init__.__doc__
    if doc and ("CudaUVMSpace" in doc or "CudaSpace" in doc):
        return kokkos.CudaUVMSpace, kokkos.LayoutLeft
    return kokkos.HostSpace, kokkos.LayoutRight


def kk_zeros(shape, dtype, memspace, layout):
    """Allocate a Kokkos array (portable CPU/GPU) and zero it via a numpy view."""
    arr = kokkos.array(list(shape), dtype=dtype, space=memspace, layout=layout)
    np.array(arr, copy=False)[:] = 0
    return arr


def build_bilayer_model(order, ex, ey, ez, lx, ly, lz, vel_top, vel_bottom,
                         z_elem_interface, memspace, layout):
    """Build a water/rock bilayer ModelUnstruct directly in Python (no C++ builder).

    Adapted from examples/fe/solver_cartesian_bilayer.py:create_bilayer_model_data_from_params.
    Returns (model_obj, nodes_coords_x, nodes_coords_y, nodes_coords_z) -- the coordinate
    arrays are kept so callers can locate the nearest global node for receivers.
    """
    n_elem = ex * ey * ez

    mesh_model_vp = np.empty((n_elem,), dtype=np.float32)
    mesh_model_rho = np.ones((n_elem,), dtype=np.float32)  # acoustic: rho unused, keep at 1.0
    for k in range(ez):
        for j in range(ey):
            for i in range(ex):
                e = i + j * ex + k * ex * ey
                mesh_model_vp[e] = vel_bottom if k < z_elem_interface else vel_top

    nx_nodes, ny_nodes, nz_nodes = order * ex + 1, order * ey + 1, order * ez + 1
    n_node = nx_nodes * ny_nodes * nz_nodes
    dx_node, dy_node, dz_node = lx / (order * ex), ly / (order * ey), lz / (order * ez)

    nodes_coords_x = np.empty(n_node, dtype=np.float32)
    nodes_coords_y = np.empty(n_node, dtype=np.float32)
    nodes_coords_z = np.empty(n_node, dtype=np.float32)
    for k in range(nz_nodes):
        for j in range(ny_nodes):
            for i in range(nx_nodes):
                idx = i + j * nx_nodes + k * nx_nodes * ny_nodes
                nodes_coords_x[idx] = i * dx_node
                nodes_coords_y[idx] = j * dy_node
                nodes_coords_z[idx] = k * dz_node

    n_gll_per_elem = (order + 1) ** 3
    elem_to_nodes = np.empty((n_elem, n_gll_per_elem), dtype=np.int32)
    for ez_i in range(ez):
        for ey_i in range(ey):
            for ex_i in range(ex):
                e = ex_i + ey_i * ex + ez_i * ex * ey
                for k_loc in range(order + 1):
                    for j_loc in range(order + 1):
                        for i_loc in range(order + 1):
                            i_glob = ex_i * order + i_loc
                            j_glob = ey_i * order + j_loc
                            k_glob = ez_i * order + k_loc
                            node_idx = i_glob + j_glob * nx_nodes + k_glob * nx_nodes * ny_nodes
                            local_idx = i_loc + j_loc * (order + 1) + k_loc * (order + 1) ** 2
                            elem_to_nodes[e, local_idx] = node_idx

    kk_elem_to_nodes = kokkos.array([n_elem, n_gll_per_elem], dtype=kokkos.int32, space=memspace, layout=layout)
    np.array(kk_elem_to_nodes, copy=False)[:] = elem_to_nodes

    kk_nodes_coords_x = kk_zeros((n_node,), kokkos.float32, memspace, layout)
    kk_nodes_coords_y = kk_zeros((n_node,), kokkos.float32, memspace, layout)
    kk_nodes_coords_z = kk_zeros((n_node,), kokkos.float32, memspace, layout)
    np.array(kk_nodes_coords_x, copy=False)[:] = nodes_coords_x
    np.array(kk_nodes_coords_y, copy=False)[:] = nodes_coords_y
    np.array(kk_nodes_coords_z, copy=False)[:] = nodes_coords_z

    kk_model_vp_element = kk_zeros((n_elem,), kokkos.float32, memspace, layout)
    kk_model_rho_element = kk_zeros((n_elem,), kokkos.float32, memspace, layout)
    np.array(kk_model_vp_element, copy=False)[:] = mesh_model_vp
    np.array(kk_model_rho_element, copy=False)[:] = mesh_model_rho

    # Acoustic, model-on-elements: node-side property arrays and every elastic/anisotropic
    # field are unused by the solver -- dummy 1-element placeholders, per
    # solver_cartesian_bilayer.py's proven pattern.
    kk_dummy_1d = kk_zeros((1,), kokkos.float32, memspace, layout)
    kk_empty_1d_int = kokkos.array([0], dtype=kokkos.int32, space=memspace, layout=layout)
    kk_dummy_3d = kk_zeros((1, 1, 1), kokkos.float32, memspace, layout)

    model_data = ModelUnstructData(
        order, n_elem, n_node, lx, ly, lz,
        True,   # is_model_on_nodes (matches solver_cartesian_bilayer.py; solver still reads
                # the element-side arrays via ModelLocationType.ONELEMENTS)
        False,  # is_elastic
        kk_elem_to_nodes,
        kk_nodes_coords_x, kk_nodes_coords_y, kk_nodes_coords_z,
        kk_dummy_1d, kk_model_vp_element,      # vp: node (dummy), element (real)
        kk_dummy_1d, kk_model_rho_element,     # rho: node (dummy), element (real)
        kk_dummy_1d, kk_dummy_1d,              # vs: node, element
        kk_dummy_1d, kk_dummy_1d,              # delta: node, element
        kk_dummy_1d, kk_dummy_1d,              # epsilon: node, element
        kk_dummy_1d, kk_dummy_1d,              # gamma: node, element
        kk_dummy_1d, kk_dummy_1d,              # theta: node, element
        kk_dummy_1d, kk_dummy_1d,              # phi: node, element
        kk_dummy_3d,                           # C_tensor_element
        kk_empty_1d_int,                       # boundaries_t
    )

    model_obj = ModelUnstruct(model_data)
    # CRITICAL: replicates what CartesianUnstructBuilder::getModel() does in C++.
    # Skipping this leaves internal face/mesh structures uninitialised.
    model_obj.build_face_connectivity()

    # Keep the kokkos arrays alive for as long as model_obj is used: ModelUnstructData's
    # bindings wrap these buffers via python_view_type_t (zero-copy), which may not extend
    # their lifetime independently of the Python objects that created them. Letting them be
    # garbage-collected once this function returns risks a use-after-free the first time the
    # C++ side actually reads through them (e.g. DGSEMsolver::computeFEInit -> globalNodeIndex).
    keepalive = (kk_elem_to_nodes, kk_nodes_coords_x, kk_nodes_coords_y, kk_nodes_coords_z,
                 kk_model_vp_element, kk_model_rho_element, kk_dummy_1d, kk_empty_1d_int, kk_dummy_3d)

    return model_obj, nodes_coords_x, nodes_coords_y, nodes_coords_z, keepalive


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


def _dense_ix_map(ix, ex, order, dense_order):
    """Map x-node index `ix` on the dense_order grid to (element_x, local_dof) in a region's
    OWN order. Nearest-neighbor when order != dense_order (visual check only, not
    quantitatively accurate) -- same trick as sem_proxy.cc's pMin/pMax x-resampling.
    Identity when order == dense_order.
    """
    nx_dense = dense_order * ex + 1
    ix_e = ex - 1 if ix == nx_dense - 1 else ix // dense_order
    ix_d_dense = dense_order if ix == nx_dense - 1 else ix % dense_order
    ix_d = min(order, round(ix_d_dense / dense_order * order))
    return ix_e, ix_d


def dg_region_rows(field_curr, order, ex, ey, ez_start, ez_count, dense_order):
    """X-z rows (mid-Y, one row per (z-element, z-dof)) of a per-element DG/DGPAdaptive field
    region, resampled onto the dense_order x-grid. Port of sem_proxy.cc's slice_dgsem_xz /
    slice_dgpadaptive_xz per-region loops.
    """
    n1d = order + 1
    ey_mid, ib_mid = ey // 2, order // 2
    nx_dense = dense_order * ex + 1
    rows = []
    for ez_i in range(ez_start, ez_start + ez_count):
        for iz_dof in range(n1d):
            row = np.empty(nx_dense, dtype=field_curr.dtype)
            for ix in range(nx_dense):
                ix_e, ix_d = _dense_ix_map(ix, ex, order, dense_order)
                elem = ix_e + ey_mid * ex + ez_i * ex * ey
                dof = ix_d + ib_mid * n1d + iz_dof * n1d * n1d
                row[ix] = field_curr[elem, dof]
            rows.append(row)
    return rows


def sem_region_rows(sem_curr, order, ex, ey, ny, iz_start, iz_end, dense_order):
    """X-z rows (mid-Y, one row per z-node) of a nodal SEM field region, resampled onto the
    dense_order x-grid (same per-element nearest-neighbor mapping as dg_region_rows, so the
    seam with the neighboring DG region lines up column-for-column).
    """
    nx_native = order * ex + 1
    nx_dense = dense_order * ex + 1
    iy_mid = ny // 2
    rows = []
    for iz in range(iz_start, iz_end):
        row = np.empty(nx_dense, dtype=sem_curr.dtype)
        for ix in range(nx_dense):
            ix_e, ix_d = _dense_ix_map(ix, ex, order, dense_order)
            native_ix = ix_e * order + ix_d
            row[ix] = sem_curr[native_ix + iy_mid * nx_native + iz * nx_native * ny]
        rows.append(row)
    return rows


def main():
    kokkos.initialize()
    try:
        run()
    finally:
        kokkos.finalize()


def run():
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
    # Receiver: nearest global SEM_high node to RCV_COORD (water region).
    # -------------------------------------------------------------------
    if RCV_COORD[2] > WATER_ROCK_ZBOUNDARY:
        order_rcv = ORDER_MAX
        rcv_local_z = RCV_COORD[2] - WATER_ROCK_ZBOUNDARY
        dz = LZ_top / (ORDER_MAX * EZ_top)
        wavefield_rcv = wavefield_top
    else:
        order_rcv = ORDER_MIN
        rcv_local_z = LZ_bot - RCV_COORD[2]
        dz = LZ_bot / (ORDER_MIN * EZ_bot)
        wavefield_rcv = wavefield_bot
        
    nx_nodes = order_rcv * EX + 1
    ny_nodes = order_rcv * EY + 1
    ny_nodes_bot = ORDER_MIN * EY + 1
    dx = LX / (order_rcv * EX)
    dy = LY / (order_rcv * EY)
    
    i = round(RCV_COORD[0] / dx)
    j = round(RCV_COORD[1] / dy)
    k = round(rcv_local_z / dz)   
    rcv_node = i + j * nx_nodes + k * nx_nodes * ny_nodes
    rcv_trace = np.zeros(N_SAMPLES, dtype=np.float32)


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

    for it in range(N_SAMPLES):
        solver_top.compute_one_step(DT, it, data_top)
        solver_mid.compute_one_step(DT, it, data_mid)
        solver_bot.compute_one_step(DT, it, data_bot)

        pn_top_dg_curr_np = np.array(wavefield_top.get_dg_current_field(0), copy=False)
        pn_bot_dg_curr_np = np.array(wavefield_bot.get_dg_current_field(0), copy=False)
        pn_mid_pmax_curr_np = np.array(wavefield_mid.get_pmax_current_field(0), copy=False)
        pn_mid_pmin_curr_np = np.array(wavefield_mid.get_pmin_current_field(0), copy=False)

        pn_top_dg_curr_np[np.ix_(top_ghost, face_high_max)]      = pn_mid_pmax_curr_np[np.ix_(mid_pmax_real, face_high_max)]
        pn_mid_pmax_curr_np[np.ix_(mid_pmax_ghost, face_low_max)] = pn_top_dg_curr_np[np.ix_(top_real, face_low_max)]

        pn_bot_dg_curr_np[np.ix_(bot_ghost, face_high_min)]       = pn_mid_pmin_curr_np[np.ix_(mid_pmin_real, face_low_min)]
        pn_mid_pmin_curr_np[np.ix_(mid_pmin_ghost, face_high_min)] = pn_bot_dg_curr_np[np.ix_(bot_real, face_low_min)]

        data_top.swap_wavefields()
        data_mid.swap_wavefields()
        data_bot.swap_wavefields()

        sem_rcv_prev_np = np.array(wavefield_rcv.get_sem_previous_field(0), copy=False)
        rcv_trace[it] = sem_rcv_prev_np[rcv_node]

        sem_top_prev_np = np.array(wavefield_top.get_sem_previous_field(0), copy=False)
        sem_bot_prev_np = np.array(wavefield_bot.get_sem_previous_field(0), copy=False)
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


    np.savetxt("bilayer_sem_padaptive_receiver_trace.txt",
               np.column_stack([t, rcv_trace]),
               header="time pressure_at_receiver (bilayer p-adative SEM, python orchestration)")
    print("Wrote bilayer_sem_padaptive_receiver_trace.txt")


if __name__ == "__main__":
    main()