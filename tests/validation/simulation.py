"""
Runs one validation case with the pyfuntides SEM solver and records the
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


def build_solver(case):
    solver = Solver.create_solver(
        Solver.MethodType.SEM,
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.STRUCT if case.mesh == "cartesian" else Solver.MeshType.UNSTRUCT,
        Solver.ModelLocationType.ONNODES if case.on_nodes else Solver.ModelLocationType.ONELEMENTS,
        Solver.PhysicType.ELASTIC if case.is_elastic else Solver.PhysicType.ACOUSTIC,
        case.order,
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


def _zeros(shape, dtype, space, layout):
    view = kokkos.array(shape, dtype=dtype, space=space, layout=layout)
    np.array(view, copy=False)[...] = 0
    return view


def run(case, memspace=None, layout=None):
    """
    Simulate one case.

    Returns (time, trace, (vp, vs, rho)): trace has one column per entry of
    case.components, and row k holds the field at time (k + 1) * dt, i.e. after
    the k-th update.
    """
    if memspace is None:
        memspace, layout = default_memspace()
    n_elements = (case.elements_per_direction,) * 3
    n_samples = case.n_samples
    n_fields = len(case.components)

    model = build_model(case)
    solver = build_solver(case)
    solver.compute_fe_init(model, [0.0, 0.0, 0.0], False, 0.015)

    src_element = element_index(cases.SOURCE, n_elements, cases.DOMAIN_SIZE)
    rcv_element = element_index(cases.RECEIVER, n_elements, cases.DOMAIN_SIZE)
    rcv_nodes = element_nodes(model, rcv_element, case.order)
    rcv_weights = interpolation_weights(model, rcv_element, case.order, cases.RECEIVER)
    material = material_at_source(model, case, src_element)

    n_points = model.get_number_of_points_per_element()
    rhs_element = _zeros([1], kokkos.int32, memspace, layout)
    np.array(rhs_element, copy=False)[0] = src_element
    rhs_weights = _zeros([1, n_points], kokkos.float32, memspace, layout)
    np.array(rhs_weights, copy=False)[0, :] = interpolation_weights(model, src_element, case.order, cases.SOURCE)

    # Same sampling as SourceTimeFunction::computeSourceTerm: sample k is the source at time k * dt.
    wavelet = evaluate_ricker(np.arange(n_samples) * cases.DT, cases.F0, cases.RICKER_ORDER, cases.TPEAK)
    rhs_terms = []
    for _ in range(n_fields):
        term = _zeros([1, n_samples], kokkos.float32, memspace, layout)
        np.array(term, copy=False)[0, :] = wavelet
        rhs_terms.append(term)

    n_nodes = model.get_number_of_nodes()
    fields = [_zeros([n_nodes], kokkos.float32, memspace, layout) for _ in range(2 * n_fields)]
    if case.is_elastic:
        wavefield = Solver.WavefieldElastic(*fields)
        rhs = Solver.RhsElastic(*rhs_terms, rhs_element, rhs_weights)
        data = Solver.SEMsolverDataElastic(wavefield, rhs)
    else:
        wavefield = Solver.WavefieldAcoustic(*fields)
        rhs = Solver.RhsAcoustic(rhs_terms[0], rhs_element, rhs_weights)
        data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    trace = np.zeros((n_samples, n_fields))
    for k in range(n_samples):
        solver.compute_forces(cases.DT, k, data)
        solver.update_solution_forward(cases.DT, data)
        # The update writes the new field into the previous-field buffer: read it once swapped.
        data.swap_wavefields()
        for c in range(n_fields):
            field = np.array(data.get_current_field(c), copy=False)
            trace[k, c] = field[rcv_nodes] @ rcv_weights

    time = (np.arange(n_samples) + 1) * cases.DT
    return time, trace, material
