import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest
import solver_utils as Utils


class UnstructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 10
        self.lx = self.ly = self.lz = 1500
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz


@pytest.fixture
def unstruct(request):
    order, param_cls, builder_cls, on_nodes, is_elastic = request.param

    sd = UnstructData(order)

    params = param_cls()

    params.ex, params.ey, params.ez = sd.ex, sd.ey, sd.ez
    params.lx, params.ly, params.lz = sd.lx, sd.ly, sd.lz
    params.order = order
    params.is_model_on_nodes = on_nodes
    params.is_elastic = is_elastic
    builder = builder_cls(params)
    return sd, params, builder, on_nodes, is_elastic


test_cases_struct = [
    # f32, i32 cases (only ones supported by solver so far)
    (
        1,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        True,
        True,
    ),
    (
        1,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        False,
        True,
    ),
    (
        1,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        True,
        False,
    ),
    (
        1,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        False,
        False,
    ),
    (
        2,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        True,
        True,
    ),
    (
        2,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        False,
        True,
    ),
    (
        2,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        True,
        False,
    ),
    (
        2,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        False,
        False,
    ),
    (
        3,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        True,
        True,
    ),
    (
        3,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        False,
        True,
    ),
    (
        3,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        True,
        False,
    ),
    (
        3,
        Model.CartesianParams_f32_i32,
        Model.CartesianUnstructBuilder_f32_i32,
        False,
        False,
    ),
]


class TestSolverUnstruct:
    @pytest.mark.parametrize("unstruct", test_cases_struct, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.MAKUTU, Solver.ImplemType.SHIVA]
    )
    def test_solver_one_step(self, unstruct, implem):
        sd, _, builder, is_model_on_nodes, is_elastic = (
            unstruct  # Déstructure avec le booléen
        )
        n_rhs = 2
        dt = 0.001
        time_sample = 1
        n_time_steps = 1
        f0 = 5.0
        model = builder.get_model()

        # TODO remove when we reactivate SHIVA
        if implem == Solver.ImplemType.SHIVA:
            return

        # Convertir le booléen en ModelLocationType
        model_location = (
            Solver.ModelLocationType.ONNODES
            if is_model_on_nodes
            else Solver.ModelLocationType.ONELEMENTS
        )

        physic_type = (
            Solver.PhysicType.ACOUSTIC if not is_elastic else Solver.PhysicType.ELASTIC
        )

        solver = Solver.create_solver(
            Solver.MethodType.SEM,
            implem,
            Solver.MeshType.UNSTRUCT,
            model_location,
            physic_type,
            sd.order,
        )

        solver.compute_fe_init(model)
        if physic_type == Solver.PhysicType.ACOUSTIC:
            kk_pnGlobal, _ = Utils.allocate_pressure(sd.n_dof)
            kk_RHSElement, _ = Utils.allocate_rhs_element(n_rhs, sd.ex, sd.ey, sd.ez)
            kk_RHSWeights, _ = Utils.allocate_rhs_weight(n_rhs, model)
            kk_RHSTerm, _ = Utils.allocate_rhs_term(n_rhs, n_time_steps, dt, f0)
            data = Solver.SEMsolverDataAcoustic(
                0, 1, kk_RHSTerm, kk_pnGlobal, kk_RHSElement, kk_RHSWeights
            )
        else:
            kk_uxnGlobal, _ = Utils.allocate_displacementx(sd.n_dof)
            kk_uynGlobal, _ = Utils.allocate_displacementy(sd.n_dof)
            kk_uznGlobal, _ = Utils.allocate_displacementz(sd.n_dof)
            kk_RHSElement, _ = Utils.allocate_rhs_element(n_rhs, sd.ex, sd.ey, sd.ez)
            kk_RHSWeights, _ = Utils.allocate_rhs_weight(n_rhs, model)
            kk_RHSTermx, _ = Utils.allocate_rhs_term(n_rhs, n_time_steps, dt, f0)
            kk_RHSTermy, _ = Utils.allocate_rhs_term(n_rhs, n_time_steps, dt, f0)
            kk_RHSTermz, _ = Utils.allocate_rhs_term(n_rhs, n_time_steps, dt, f0)
            data = Solver.SEMsolverDataElastic(
                0,
                1,
                kk_RHSTermx,
                kk_RHSTermy,
                kk_RHSTermz,
                kk_uxnGlobal,
                kk_uynGlobal,
                kk_uznGlobal,
                kk_RHSElement,
                kk_RHSWeights,
            )
        solver.compute_one_step(dt, time_sample, data)
