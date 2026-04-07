import pyfuntides.gradient as Gradient
import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest
import gradient_utils as Utils


class UnstructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 10
        self.lx = self.ly = self.lz = 1500
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz
        self.n_elements = self.ex * self.ey * self.ez


@pytest.fixture
def unstruct(request):
    order, param_cls, builder_cls, on_nodes = request.param
    sd = UnstructData(order)
    params = param_cls()
    params.ex, params.ey, params.ez = sd.ex, sd.ey, sd.ez
    params.lx, params.ly, params.lz = sd.lx, sd.ly, sd.lz
    params.order = order
    params.is_model_on_nodes = on_nodes
    params.is_elastic = True
    builder = builder_cls(params)
    return sd, builder, on_nodes


test_cases = [
    (1, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, True),
    (1, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, False),
    (2, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, True),
    (2, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, False),
    (3, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, True),
    (3, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32, False),
]


class TestDifferentiatorElasticUnstruct:
    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    @pytest.mark.parametrize("implem", [Solver.ImplemType.MAKUTU])
    def test_differentiator_compute(self, unstruct, implem):
        sd, builder, is_model_on_nodes = unstruct
        model = builder.get_model()

        model_location = (
            Solver.ModelLocationType.ONNODES
            if is_model_on_nodes
            else Solver.ModelLocationType.ONELEMENTS
        )

        differentiator = Gradient.create_differentiator(
            implem,
            Solver.MeshType.UNSTRUCT,
            model_location,
            Solver.PhysicType.ELASTIC,
            sd.order,
        )

        n_nodes = sd.n_dof
        n_grad = sd.n_dof if is_model_on_nodes else sd.n_elements

        kk_ux_fwd, _ = Utils.allocate_field(n_nodes)
        kk_uy_fwd, _ = Utils.allocate_field(n_nodes)
        kk_uz_fwd, _ = Utils.allocate_field(n_nodes)
        kk_ux_adj, _ = Utils.allocate_field(n_nodes)
        kk_uy_adj, _ = Utils.allocate_field(n_nodes)
        kk_uz_adj, _ = Utils.allocate_field(n_nodes)
        kk_ux_dt2, _ = Utils.allocate_field(n_nodes)
        kk_uy_dt2, _ = Utils.allocate_field(n_nodes)
        kk_uz_dt2, _ = Utils.allocate_field(n_nodes)
        kk_grad_rho, _ = Utils.allocate_field(n_grad)
        kk_grad_lambda, _ = Utils.allocate_field(n_grad)
        kk_grad_mu, _ = Utils.allocate_field(n_grad)

        dt = 0.001
        fwd = Gradient.WavefieldViewForwardElastic(
            kk_ux_fwd, kk_uy_fwd, kk_uz_fwd
        )
        bwd = Gradient.WavefieldViewBackwardElastic(
            kk_ux_adj, kk_uy_adj, kk_uz_adj,
            kk_ux_dt2, kk_uy_dt2, kk_uz_dt2,
        )
        grad = Gradient.GradientElastic(
            kk_grad_rho, kk_grad_lambda, kk_grad_mu
        )
        data = Gradient.GradientDataElastic(fwd, bwd, grad)

        differentiator.compute(model, data, dt)
