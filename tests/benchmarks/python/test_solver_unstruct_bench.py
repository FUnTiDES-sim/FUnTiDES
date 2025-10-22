import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest
import solver_utils as Utils
import benchmark_groups as Groups


class UnstructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 100
        self.lx = self.ly = self.lz = 2000
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz


@pytest.fixture
def unstruct(request):
    order, param_cls, builder_cls, on_nodes = request.param

    sd = UnstructData(order)

    params = param_cls()
    params.ex, params.ey, params.ez = sd.ex, sd.ey, sd.ez
    params.lx, params.ly, params.lz = sd.lx, sd.ly, sd.lz
    params.order = order
    params.is_model_on_nodes = on_nodes

    builder = builder_cls(params)

    return sd, params, builder


test_cases = [
    # f32, i32 cases (only ones supported by solver so far)
    (1, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, True),
    (1, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, False),
    (2, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, True),
    (2, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, False),
    (3, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, True),
    (3, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, False),
]


class TestSolverUnstruct:
    @pytest.mark.benchmark(group=Groups.BenchmarkGroup.COMPUTE_FE_INIT.name)
    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.MAKUTU, Solver.ImplemType.SHIVA]
    )
    def test_solver_fe_init(self, unstruct, implem, benchmark):
        sd, _, builder = unstruct

        model = builder.get_model()

        # TODO remove when we reactivate SHIVA
        if implem == Solver.ImplemType.SHIVA:
            return

        solver = Solver.create_solver(
            Solver.MethodType.SEM, implem, Solver.MeshType.UNSTRUCT, sd.order
        )

        benchmark(solver.compute_fe_init, model)

    @pytest.mark.benchmark(group=Groups.BenchmarkGroup.COMPUTE_ONE_STEP.name)
    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.MAKUTU, Solver.ImplemType.SHIVA]
    )
    def test_solver_one_step(self, unstruct, implem, benchmark):
        sd, _, builder = unstruct
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
            Solver.MethodType.SEM, implem, Solver.MeshType.UNSTRUCT, sd.order
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
