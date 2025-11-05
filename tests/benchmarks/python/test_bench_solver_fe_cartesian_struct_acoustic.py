import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest
import solver_utils as Utils
import benchmark_groups as Groups
from data_structures import StructData


@pytest.fixture
def struct(request):
    order, builder_cls, on_nodes, is_elastic = request.param

    sd = StructData(order)

    builder = builder_cls(
        sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, on_nodes, is_elastic
    )

    return sd, builder, on_nodes, is_elastic


test_cases = [
    # f32, i32 cases (only ones supported by solver so far)
    (1, Model.CartesianStructBuilder_f32_i32_O1, True, False),
    (1, Model.CartesianStructBuilder_f32_i32_O1, False, False),
    (2, Model.CartesianStructBuilder_f32_i32_O2, True, False),
    (2, Model.CartesianStructBuilder_f32_i32_O2, False, False),
    (3, Model.CartesianStructBuilder_f32_i32_O3, True, False),
    (3, Model.CartesianStructBuilder_f32_i32_O3, False, False),
]


class TestSolverStructAcoustic:
    @pytest.mark.benchmark(group=Groups.BenchmarkGroup.COMPUTE_FE_INIT.name)
    @pytest.mark.parametrize("struct", test_cases, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.MAKUTU, Solver.ImplemType.SHIVA]
    )
    def test_solver_fe_init(self, struct, implem, benchmark):
        sd, builder, on_nodes, is_elastic = struct

        model = builder.get_model()

        # TODO remove when we reactivate SHIVA
        if implem == Solver.ImplemType.SHIVA:
            return

        model_location = (
            Solver.ModelLocationType.ONNODES
            if on_nodes
            else Solver.ModelLocationType.ONELEMENTS
        )

        physic_type = Solver.PhysicType.ACOUSTIC

        solver = Solver.create_solver(
            Solver.MethodType.SEM,
            implem,
            Solver.MeshType.STRUCT,
            model_location,
            physic_type,
            sd.order,
        )

        benchmark(solver.compute_fe_init, model)

    @pytest.mark.benchmark(group=Groups.BenchmarkGroup.COMPUTE_ONE_STEP.name)
    @pytest.mark.parametrize("struct", test_cases, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.MAKUTU, Solver.ImplemType.SHIVA]
    )
    def test_solver_one_step(self, struct, implem, benchmark):
        sd, builder, on_nodes, is_elastic = struct
        n_rhs = 2
        dt = 0.001
        time_sample = 1
        n_time_steps = 1500
        f0 = 5.0

        model = builder.get_model()

        # TODO remove when we reactivate SHIVA
        if implem == Solver.ImplemType.SHIVA:
            return

        model_location = (
            Solver.ModelLocationType.ONNODES
            if on_nodes
            else Solver.ModelLocationType.ONELEMENTS
        )

        physic_type = (
            Solver.PhysicType.ACOUSTIC if not is_elastic else Solver.PhysicType.ELASTIC
        )

        solver = Solver.create_solver(
            Solver.MethodType.SEM,
            implem,
            Solver.MeshType.STRUCT,
            model_location,
            physic_type,
            sd.order,
        )

        solver.compute_fe_init(model)

        kk_pnGlobal, _ = Utils.allocate_pressure(sd.n_dof)
        kk_RHSElement, _ = Utils.allocate_rhs_element(n_rhs, sd.ex, sd.ey, sd.ez)
        kk_RHSWeights, _ = Utils.allocate_rhs_weight(n_rhs, model)
        kk_RHSTerm, _ = Utils.allocate_rhs_term(n_rhs, n_time_steps, dt, f0)

        data = Solver.SEMsolverDataAcoustic(
            0, 1, kk_RHSTerm, kk_pnGlobal, kk_RHSElement, kk_RHSWeights
        )

        benchmark(solver.compute_one_step, dt, time_sample, data)
