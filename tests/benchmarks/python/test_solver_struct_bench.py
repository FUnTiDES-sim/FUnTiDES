import pyproxys.model as Model
import pyproxys.solver as Solver
import pytest
import solver_utils as Utils
import benchmark_groups as Groups


class StructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 100
        self.domain_size = 2000
        self.hx = self.domain_size / self.ex
        self.hy = self.domain_size / self.ey
        self.hz = self.domain_size / self.ez
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz


@pytest.fixture
def struct(request):
    order, builder_cls, on_nodes = request.param

    sd = StructData(order)

    builder = builder_cls(sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, on_nodes)

    return sd, builder


test_cases = [
    # f32, i32 cases (only ones supported by solver so far)
    (1, Model.CartesianStructBuilder_f32_i32_O1, True),
    (1, Model.CartesianStructBuilder_f32_i32_O1, False),
    (2, Model.CartesianStructBuilder_f32_i32_O2, True),
    (2, Model.CartesianStructBuilder_f32_i32_O2, False),
    (3, Model.CartesianStructBuilder_f32_i32_O3, True),
    (3, Model.CartesianStructBuilder_f32_i32_O3, False),
]


class TestSolverStruct:
    @pytest.mark.benchmark(group=Groups.BenchmarkGroup.COMPUTE_FE_INIT.name)
    @pytest.mark.parametrize("struct", test_cases, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.GEOS, Solver.ImplemType.SHIVA]
    )
    def test_solver_fe_init(self, struct, implem, benchmark):
        sd, builder = struct

        model = builder.get_model()

        # TODO remove when we reactivate SHIVA
        if implem == Solver.ImplemType.SHIVA:
            return

        solver = Solver.create_solver(
            Solver.MethodType.SEM, implem, Solver.MeshType.STRUCT, sd.order
        )

        benchmark(solver.compute_fe_init, model)

    @pytest.mark.benchmark(group=Groups.BenchmarkGroup.COMPUTE_ONE_STEP.name)
    @pytest.mark.parametrize("struct", test_cases, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.GEOS, Solver.ImplemType.SHIVA]
    )
    def test_solver_one_step(self, struct, implem, benchmark):
        sd, builder = struct
        n_rhs = 2
        dt = 0.001
        time_sample = 1
        n_time_steps = 1
        f0 = 5.0

        model = builder.get_model()

        # TODO remove when we reactivate SHIVA
        if implem == Solver.ImplemType.SHIVA:
            return

        solver = Solver.create_solver(
            Solver.MethodType.SEM, implem, Solver.MeshType.STRUCT, sd.order
        )

        solver.compute_fe_init(model)

        kk_pnGlobal, _ = Utils.allocate_pressure(sd.n_dof)
        kk_RHSElement, _ = Utils.allocate_rhs_element(n_rhs, sd.ex, sd.ey, sd.ez)
        kk_RHSWeights, _ = Utils.allocate_rhs_weight(n_rhs, model)
        kk_RHSTerm, _ = Utils.allocate_rhs_term(n_rhs, n_time_steps, dt, f0)

        data = Solver.SEMsolverData(
            0, 1, kk_RHSTerm, kk_pnGlobal, kk_RHSElement, kk_RHSWeights
        )

        benchmark(solver.compute_one_step, dt, time_sample, data)
