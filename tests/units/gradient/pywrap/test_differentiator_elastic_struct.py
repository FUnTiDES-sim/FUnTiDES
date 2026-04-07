import pyfuntides.gradient as Gradient
import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest
import gradient_utils as Utils


class StructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 10
        self.domain_size = 1500
        self.hx = self.domain_size / self.ex
        self.hy = self.domain_size / self.ey
        self.hz = self.domain_size / self.ez
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz
        self.n_elements = self.ex * self.ey * self.ez


@pytest.fixture
def struct(request):
    order, builder_cls, on_nodes = request.param
    sd = StructData(order)
    builder = builder_cls(
        sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, on_nodes, True
    )
    return sd, builder, on_nodes


test_cases = [
    (1, Model.CartesianStructBuilder_f32_i32_O1, True),
    (1, Model.CartesianStructBuilder_f32_i32_O1, False),
    (2, Model.CartesianStructBuilder_f32_i32_O2, True),
    (2, Model.CartesianStructBuilder_f32_i32_O2, False),
    (3, Model.CartesianStructBuilder_f32_i32_O3, True),
    (3, Model.CartesianStructBuilder_f32_i32_O3, False),
]


class TestDifferentiatorElasticStruct:
    @pytest.mark.parametrize("struct", test_cases, indirect=True)
    @pytest.mark.parametrize("implem", [Solver.ImplemType.MAKUTU])
    def test_differentiator_compute(self, struct, implem):
        sd, builder, is_model_on_nodes = struct
        model = builder.get_model()

        model_location = (
            Solver.ModelLocationType.ONNODES
            if is_model_on_nodes
            else Solver.ModelLocationType.ONELEMENTS
        )

        differentiator = Gradient.create_differentiator(
            implem,
            Solver.MeshType.STRUCT,
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
