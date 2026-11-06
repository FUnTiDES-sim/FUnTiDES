"""2-D acoustic SEM solver through pyfuntides (built with COMPILE_SEM_2D)."""

import kokkos
import numpy as np
import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest

pytestmark = pytest.mark.skipif(
    not hasattr(Model, "CartesianUnstruct2DBuilder_f32_i32"),
    reason="FUnTiDES built without COMPILE_SEM_2D",
)

VP = 1500.0  # default vp of the 2-D builder (rho = 1)


def array(shape, dtype=kokkos.float32):
    """Zero Kokkos array in the memory space and layout of the FUnTiDES build."""
    space = kokkos.DefaultMemorySpace
    layout = kokkos.LayoutRight if space == kokkos.HostSpace else kokkos.LayoutLeft
    kk = kokkos.array(list(shape), dtype=dtype, space=space, layout=layout)
    np.array(kk, copy=False)[...] = 0
    return kk


def build_model(order, ex, ez, lx, lz, free_surface=True):
    params = Model.CartesianParams_f32_i32(order, ex, 1, ez, lx, 1.0, lz, False, False)
    return Model.CartesianUnstruct2DBuilder_f32_i32(params).get_model(free_surface)


def make_solver(order, model):
    solver = Solver.create_solver(
        Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.UNSTRUCT,
        Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ACOUSTIC, order, dim=2)
    solver.compute_fe_init(model)
    return solver


@pytest.mark.parametrize("order", [1, 2, 3])
def test_model_and_mass(order):
    model = build_model(order, 8, 4, 2000.0, 1000.0)
    assert model.get_dimension() == 2
    assert model.get_order() == order
    assert model.get_number_of_points_per_element() == (order + 1) ** 2
    assert model.get_number_of_nodes() == (8 * order + 1) * (4 * order + 1)
    solver = make_solver(order, model)
    mass = np.array(solver.get_mass_matrix(), copy=False)
    assert mass.sum() == pytest.approx(2000.0 * 1000.0 / VP**2, rel=1e-4)


def test_point_source_shot():
    """Ricker source at the centre: left/right symmetric arrivals at the expected
    time, free-surface nodes at zero."""
    order, ex, ez, L, H = 3, 21, 11, 2100.0, 1100.0
    model = build_model(order, ex, ez, L, H, free_surface=True)
    solver = make_solver(order, model)
    n_node = model.get_number_of_nodes()
    npe = model.get_number_of_points_per_element()

    dt = 0.2 * model.get_min_spacing() / VP
    n_steps = int(0.6 / dt)
    f0 = 10.0
    t = np.arange(n_steps) * dt - 1.0 / f0
    a = (np.pi * f0) ** 2
    ricker = (1 - 2 * a * t**2) * np.exp(-a * t**2)

    src_elem = ex // 2 + (ez // 2) * ex
    term, elem, weights = array([1, n_steps]), array([1], kokkos.int32), array([1, npe])
    np.array(term, copy=False)[0, :] = ricker
    np.array(elem, copy=False)[0] = src_elem
    np.array(weights, copy=False)[0, :] = 1.0 / npe
    prev, cur = array([n_node]), array([n_node])
    data = Solver.SEMsolverDataAcoustic(Solver.WavefieldAcoustic(prev, cur),
                                        Solver.RhsAcoustic(term, elem, weights))

    # Receivers: a node about d left of the source centre, and its mirror node
    x = np.array([model.node_coord(i, 0) for i in range(n_node)])
    z = np.array([model.node_coord(i, 2) for i in range(n_node)])
    xs, zs = 1050.0, 550.0
    left = int(np.argmin((x - (xs - 400.0)) ** 2 + (z - zs) ** 2))
    right = int(np.argmin((x - (2 * xs - x[left])) ** 2 + (z - z[left]) ** 2))
    assert abs(x[right] - (2 * xs - x[left])) < 1e-3 and z[right] == z[left]
    d = np.hypot(x[left] - xs, z[left] - zs)
    traces = np.zeros((n_steps, 2))
    for s in range(n_steps):
        solver.compute_one_step(dt, s, data)
        data.swap_wavefields()
        p = np.array(data.get_current_field(0), copy=False)
        traces[s] = p[left], p[right]

    assert np.isfinite(traces).all() and np.abs(traces).max() > 0
    np.testing.assert_allclose(traces[:, 0], traces[:, 1], atol=1e-5 * np.abs(traces).max())
    # 2-D Green's function: the peak follows the arrival d / vp + 1 / f0 by a fraction of a period
    t_peak = np.argmax(np.abs(traces[:, 0])) * dt
    assert d / VP + 1.0 / f0 - 0.005 < t_peak < d / VP + 1.0 / f0 + 0.05

    p = np.array(data.get_current_field(0), copy=False)
    top = np.isclose(z, H)
    assert top.any() and np.all(p[top] == 0.0)


def test_unsupported_configurations():
    model = build_model(2, 4, 4, 400.0, 400.0)
    with pytest.raises(RuntimeError):
        Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.STRUCT,
                             Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ACOUSTIC, 2, dim=2)
    with pytest.raises(RuntimeError):
        Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.UNSTRUCT,
                             Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ELASTIC, 2, dim=2)
    with pytest.raises(RuntimeError):
        Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.UNSTRUCT,
                             Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ACOUSTIC, 2, dim=5)
    # A 3-D solver rejects the 2-D mesh, and the 2-D solver a mesh of another order
    solver3d = Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.UNSTRUCT,
                                    Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ACOUSTIC, 2)
    solver2d = Solver.create_solver(Solver.MethodType.SEM, Solver.ImplemType.MAKUTU, Solver.MeshType.UNSTRUCT,
                                    Solver.ModelLocationType.ONELEMENTS, Solver.PhysicType.ACOUSTIC, 3, dim=2)
    with pytest.raises(RuntimeError):
        solver2d.compute_fe_init(model)
    with pytest.raises(RuntimeError):
        solver3d.compute_fe_init(model)


def test_model_data_dimension_attribute():
    names = [n for n in dir(Model.ModelUnstructData_f32_i32) if not n.startswith("_")]
    assert "dim" in names
