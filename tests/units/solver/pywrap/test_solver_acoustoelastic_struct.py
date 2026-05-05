"""Smoke tests for the acousto-elastic solver on a structured Cartesian mesh.

Each test builds a 10x10x10 two-layer mesh (fluid above, solid below),
runs one time step with a source in the fluid domain, and asserts that
no exception is raised.

Two-layer geometry:
  - Total height: 1500 m  (lz = 1500, hz = 150 m per element)
  - Interface:    750 m   (acousto_elastic_boundary_z = 750)
  - Source:       900 m   (z >= 750  => fluid domain)
  - Receiver:     375 m   (z <  750  => solid domain)
"""

import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest
import solver_utils as Utils


class AEStructData:
    """Geometry and discretisation parameters for the two-layer test case."""

    domain_size = 1500.0
    boundary_z = 750.0   # fluid above, solid below

    def __init__(self, order):
        self.order = order
        self.ex = self.ey = self.ez = 10
        self.hx = self.domain_size / self.ex
        self.hy = self.domain_size / self.ey
        self.hz = self.domain_size / self.ez
        self.nx = self.ex * order + 1
        self.ny = self.ey * order + 1
        self.nz = self.ez * order + 1
        self.n_dof = self.nx * self.ny * self.nz


# (order, builder_class, is_model_on_nodes)
_BUILDER = {
    (1, False): Model.CartesianStructBuilder_f32_i32_O1,
    (1, True):  Model.CartesianStructBuilder_f32_i32_O1,
    (2, False): Model.CartesianStructBuilder_f32_i32_O2,
    (2, True):  Model.CartesianStructBuilder_f32_i32_O2,
    (3, False): Model.CartesianStructBuilder_f32_i32_O3,
    (3, True):  Model.CartesianStructBuilder_f32_i32_O3,
}

test_cases = [
    (1, False),
    (1, True),
    (2, False),
    (2, True),
    (3, False),
    (3, True),
]


@pytest.fixture
def ae_struct(request):
    order, on_nodes = request.param
    sd = AEStructData(order)
    builder_cls = _BUILDER[(order, on_nodes)]
    builder = builder_cls(
        sd.ex, sd.domain_size, sd.ey, sd.domain_size, sd.ez, sd.domain_size,
        on_nodes,   # is_model_on_nodes
        True,       # is_elastic (ignored for AE — material set by AE branch)
        is_acousto_elastic=True,
        acousto_elastic_boundary_z=sd.boundary_z,
    )
    return sd, builder, on_nodes


class TestAcoustoElasticSolverStruct:
    @pytest.mark.parametrize("ae_struct", test_cases, indirect=True)
    def test_one_step(self, ae_struct):
        sd, builder, on_nodes = ae_struct

        n_rhs = 1
        dt = 0.001
        time_sample = 1
        n_time_steps = 1
        f0 = 5.0

        model = builder.get_model()

        model_location = (
            Solver.ModelLocationType.ONNODES
            if on_nodes
            else Solver.ModelLocationType.ONELEMENTS
        )

        solver = Solver.create_solver(
            Solver.MethodType.SEM,
            Solver.ImplemType.MAKUTU,
            Solver.MeshType.STRUCT,
            model_location,
            Solver.PhysicType.ACOUSTOELASTIC,
            sd.order,
        )

        solver.compute_fe_init(model)

        # Acoustic (pressure) wavefield
        kk_pn_prev, _, kk_pn_curr, _ = Utils.allocate_pressure(sd.n_dof)
        # Elastic (displacement) wavefield
        kk_uxn_prev, _, kk_uxn_curr, _ = Utils.allocate_displacementx(sd.n_dof)
        kk_uyn_prev, _, kk_uyn_curr, _ = Utils.allocate_displacementy(sd.n_dof)
        kk_uzn_prev, _, kk_uzn_curr, _ = Utils.allocate_displacementz(sd.n_dof)

        # Source element in the fluid domain (upper half of the mesh)
        kk_rhs_element, _ = Utils.allocate_rhs_element(
            n_rhs, sd.ex, sd.ey, sd.ez
        )
        kk_rhs_weights, _ = Utils.allocate_rhs_weight(n_rhs, model)

        # Acoustic source term (fluid domain); elastic terms are zero
        kk_acoustic_term, _ = Utils.allocate_rhs_term(
            n_rhs, n_time_steps, dt, f0
        )
        kk_elastic_termx, _ = Utils.allocate_rhs_term(
            n_rhs, n_time_steps, dt=0.0, f0=1.0
        )
        kk_elastic_termy, _ = Utils.allocate_rhs_term(
            n_rhs, n_time_steps, dt=0.0, f0=1.0
        )
        kk_elastic_termz, _ = Utils.allocate_rhs_term(
            n_rhs, n_time_steps, dt=0.0, f0=1.0
        )

        wavefield = Solver.WavefieldAcoustoElastic(
            kk_pn_prev,  kk_pn_curr,
            kk_uxn_prev, kk_uxn_curr,
            kk_uyn_prev, kk_uyn_curr,
            kk_uzn_prev, kk_uzn_curr,
        )
        rhs = Solver.RhsAcoustoElastic(
            kk_acoustic_term, kk_rhs_element, kk_rhs_weights,
            kk_elastic_termx, kk_elastic_termy, kk_elastic_termz,
        )
        data = Solver.SEMsolverDataAcoustoElastic(wavefield, rhs)

        # Should not raise
        solver.compute_one_step(dt, time_sample, data)
