"""
Runs one validation case with the pyfuntides solver and records the
receiver trace.

Kokkos must be initialized by the caller, and finalized only after this
module's functions have returned (every Kokkos view is local to them).
"""

import kokkos
import numpy as np
import pyfuntides.model as Model
import pyfuntides.solver as Solver

import cases
from acquisition import element_index, element_nodes, interpolation_weights
from analytical_acoustic_homogeneous import evaluate_ricker

STRUCT_BUILDERS = {
    1: Model.CartesianStructBuilder_f32_i32_O1,
    2: Model.CartesianStructBuilder_f32_i32_O2,
    3: Model.CartesianStructBuilder_f32_i32_O3,
}


def default_memspace():
    """Memory space and layout matching the views the bindings were compiled for."""
    doc = Solver.WavefieldAcoustic.__init__.__doc__ or ""
    if "CudaUVMSpace" in doc or "CudaSpace" in doc:
        return kokkos.CudaUVMSpace, kokkos.LayoutLeft
    return kokkos.HostSpace, kokkos.LayoutRight


def build_model(case):
    n = case.elements_per_direction
    size = cases.DOMAIN_SIZE
    if case.mesh == "cartesian":
        builder = STRUCT_BUILDERS[case.order](n, size, n, size, n, size, case.on_nodes, case.is_elastic)
    else:
        params = Model.CartesianParams_f32_i32()
        params.ex, params.ey, params.ez = n, n, n
        params.lx, params.ly, params.lz = size, size, size
        params.order = case.order
        params.is_model_on_nodes = case.on_nodes
        params.is_elastic = case.is_elastic
        builder = Model.CartesianUnstructBuilder_f32_i32(params)
    return builder.get_model()


METHOD_TYPES = {
    "sem": Solver.MethodType.SEM,
    "dg": Solver.MethodType.DG,
    "dgsem": Solver.MethodType.DGSEM,
    "dgpadaptive": Solver.MethodType.DGPADAPTIVE,
}


def build_solver(case):
    solver = Solver.create_solver(
        METHOD_TYPES[case.method],
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.STRUCT if case.mesh == "cartesian" else Solver.MeshType.UNSTRUCT,
        Solver.ModelLocationType.ONNODES if case.on_nodes else Solver.ModelLocationType.ONELEMENTS,
        Solver.PhysicType.ELASTIC if case.is_elastic else Solver.PhysicType.ACOUSTIC,
        case.order,
        case.order_min,
    )
    if case.is_elastic:
        solver.set_anisotropy_type(Model.AnisotropyType.ISO)
    return solver


def material_at_source(model, case, element):
    """(vp, vs, rho) of the homogeneous medium, read back from the model; vs is 0 for acoustic runs."""
    if case.on_nodes:
        node = int(element_nodes(model, element, case.order)[0])
        vp, rho = model.get_model_vp_on_node(node), model.get_model_rho_on_node(node)
        vs = model.get_model_vs_on_node(node) if case.is_elastic else 0.0
    else:
        vp, rho = model.get_model_vp_on_element(element), model.get_model_rho_on_element(element)
        vs = model.get_model_vs_on_element(element) if case.is_elastic else 0.0
    return float(vp), float(vs), float(rho)


def _allocator(memspace, layout):
    """Returns alloc(shape, values=None, dtype=float32): a Kokkos array filled with values, or zeros."""

    def alloc(shape, values=None, dtype=kokkos.float32):
        view = kokkos.array(list(shape), dtype=dtype, space=memspace, layout=layout)
        np.array(view, copy=False)[...] = 0 if values is None else values
        return view

    return alloc


def _locate(model, case, point):
    """(element, interpolation weights) of a point, in the basis of the sub-domain holding it."""
    n_elements = (case.elements_per_direction,) * 3
    element = element_index(point, n_elements, cases.DOMAIN_SIZE)
    return element, interpolation_weights(model, element, case.order, point, case.basis_order(point))


def run(case, memspace=None, layout=None):
    """
    Simulate one case.

    Returns (time, trace, (vp, vs, rho)): trace has one column per entry of
    case.components, and row k holds the field at time (k + 1) * dt, i.e. after
    the k-th update.
    """
    if memspace is None:
        memspace, layout = default_memspace()
    model = build_model(case)
    solver = build_solver(case)
    solver.compute_fe_init(model, [0.0, 0.0, 0.0], False, 0.015)
    source = _locate(model, case, case.source)
    receiver = _locate(model, case, case.receiver)
    material = material_at_source(model, case, source[0])
    # Same sampling as SourceTimeFunction::computeSourceTerm: sample k is the source at time k * dt.
    wavelet = evaluate_ricker(np.arange(case.n_samples) * case.dt, cases.F0, cases.RICKER_ORDER, cases.TPEAK)

    time_loop = {
        "sem": _time_loop_sem,
        "dg": _time_loop_dg,
        "dgsem": _time_loop_dgsem,
        "dgpadaptive": _time_loop_dgpadaptive,
    }[case.method]
    trace = time_loop(case, model, solver, source, receiver, wavelet, _allocator(memspace, layout))
    time = (np.arange(case.n_samples) + 1) * case.dt
    return time, trace, material


def _time_loop_sem(case, model, solver, source, receiver, wavelet, alloc):
    n_fields = len(case.components)
    rhs_element = alloc([1], [source[0]], kokkos.int32)
    rhs_weights = alloc([1, len(source[1])], source[1])
    rhs_terms = [alloc([1, case.n_samples], wavelet) for _ in range(n_fields)]
    fields = [alloc([model.get_number_of_nodes()]) for _ in range(2 * n_fields)]
    if case.is_elastic:
        wavefield = Solver.WavefieldElastic(*fields)
        rhs = Solver.RhsElastic(*rhs_terms, rhs_element, rhs_weights)
        data = Solver.SEMsolverDataElastic(wavefield, rhs)
    else:
        wavefield = Solver.WavefieldAcoustic(*fields)
        rhs = Solver.RhsAcoustic(rhs_terms[0], rhs_element, rhs_weights)
        data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    rcv_nodes = element_nodes(model, receiver[0], case.order)
    trace = np.zeros((case.n_samples, n_fields))
    for k in range(case.n_samples):
        solver.compute_forces(case.dt, k, data)
        solver.update_solution_forward(case.dt, data)
        # The update writes the new field into the previous-field buffer: read it once swapped.
        data.swap_wavefields()
        for c in range(n_fields):
            field = np.array(data.get_current_field(c), copy=False)
            trace[k, c] = field[rcv_nodes] @ receiver[1]
    return trace


def _alternating_loop(case, solver, data, buffers, read):
    """
    Time loop of the DG-based solvers, whose data objects have no field accessor.

    The update writes the new field into the buffer that was "previous" before the step and the
    swap exchanges the two, so the new field alternates between buffers[0] and buffers[1].
    The coupled solvers implement the step only as compute_one_step: their compute_forces and
    update_solution_forward do nothing.
    """
    trace = np.zeros((case.n_samples, 1))
    for k in range(case.n_samples):
        if case.is_coupled:
            solver.compute_one_step(case.dt, k, data)
        else:
            solver.compute_forces(case.dt, k, data)
            solver.update_solution_forward(case.dt, data)
        data.swap_wavefields()
        trace[k, 0] = read(buffers[k % 2])
    return trace


def _dg_reader(receiver):
    element, weights = receiver
    return lambda field: field[element, :] @ weights


def _time_loop_dg(case, model, solver, source, receiver, wavelet, alloc):
    rhs_element = alloc([1], [source[0]], kokkos.int32)
    rhs_weights = alloc([1, len(source[1])], source[1])
    rhs_term = alloc([1, case.n_samples], wavelet)
    shape = [model.get_number_of_elements(), model.get_number_of_points_per_element()]
    prev, curr = alloc(shape), alloc(shape)
    wavefield = Solver.DGWavefieldAcoustic(prev, curr)
    rhs = Solver.RhsAcoustic(rhs_term, rhs_element, rhs_weights)
    data = Solver.DGsolverDataAcoustic(wavefield, rhs)
    buffers = (np.array(prev, copy=False), np.array(curr, copy=False))
    return _alternating_loop(case, solver, data, buffers, _dg_reader(receiver))


def _time_loop_dgsem(case, model, solver, source, receiver, wavelet, alloc):
    src_in_dg = case.in_lower_domain(case.source)
    rhs_element = alloc([1], [source[0]], kokkos.int32)
    rhs_weights = alloc([1, len(source[1])], source[1])
    dg_term = alloc([1, case.n_samples], wavelet if src_in_dg else None)
    sem_term = alloc([1, case.n_samples], None if src_in_dg else wavelet)
    dg_shape = [model.get_number_of_elements(), model.get_number_of_points_per_element()]
    dg_prev, dg_curr = alloc(dg_shape), alloc(dg_shape)
    sem_prev, sem_curr = alloc([model.get_number_of_nodes()]), alloc([model.get_number_of_nodes()])
    wavefield = Solver.DGSEMWavefieldAcoustic(dg_prev, dg_curr, sem_prev, sem_curr)
    rhs = Solver.DGSEMRhsAcoustic(dg_term, sem_term, rhs_element, rhs_weights)
    data = Solver.DGSEMsolverData(wavefield, rhs)

    if case.in_lower_domain(case.receiver):
        buffers = (np.array(dg_prev, copy=False), np.array(dg_curr, copy=False))
        read = _dg_reader(receiver)
    else:
        buffers = (np.array(sem_prev, copy=False), np.array(sem_curr, copy=False))
        rcv_nodes = element_nodes(model, receiver[0], case.order)
        read = lambda field: field[rcv_nodes] @ receiver[1]
    return _alternating_loop(case, solver, data, buffers, read)


def _time_loop_dgpadaptive(case, model, solver, source, receiver, wavelet, alloc):
    n_elements = model.get_number_of_elements()
    n_min, n_max = (case.order_min + 1) ** 3, (case.order + 1) ** 3
    src_in_pmin = case.in_lower_domain(case.source)
    rhs_element = alloc([1], [source[0]], kokkos.int32)
    pmin_weights = alloc([1, n_min], source[1] if src_in_pmin else None)
    pmax_weights = alloc([1, n_max], None if src_in_pmin else source[1])
    pmin_term = alloc([1, case.n_samples], wavelet if src_in_pmin else None)
    pmax_term = alloc([1, case.n_samples], None if src_in_pmin else wavelet)
    pmin_prev, pmin_curr = alloc([n_elements, n_min]), alloc([n_elements, n_min])
    pmax_prev, pmax_curr = alloc([n_elements, n_max]), alloc([n_elements, n_max])
    wavefield = Solver.DGPAdaptiveWavefieldAcoustic(pmin_prev, pmin_curr, pmax_prev, pmax_curr)
    rhs = Solver.DGPAdaptiveRhsAcoustic(pmin_term, pmax_term, rhs_element, pmin_weights, pmax_weights)
    data = Solver.DGPAdaptiveSolverData(wavefield, rhs)

    if case.in_lower_domain(case.receiver):
        buffers = (np.array(pmin_prev, copy=False), np.array(pmin_curr, copy=False))
    else:
        buffers = (np.array(pmax_prev, copy=False), np.array(pmax_curr, copy=False))
    return _alternating_loop(case, solver, data, buffers, _dg_reader(receiver))
