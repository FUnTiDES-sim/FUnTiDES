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
        sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, on_nodes, False
    )

    return sd, builder, on_nodes


test_cases = [
    # f32, i32 cases (only ones supported by the differentiator so far)
    (1, Model.CartesianStructBuilder_f32_i32_O1, True),
    (1, Model.CartesianStructBuilder_f32_i32_O1, False),
    (2, Model.CartesianStructBuilder_f32_i32_O2, True),
    (2, Model.CartesianStructBuilder_f32_i32_O2, False),
    (3, Model.CartesianStructBuilder_f32_i32_O3, True),
    (3, Model.CartesianStructBuilder_f32_i32_O3, False),
]


class TestDifferentiatorAcousticStruct:
    @pytest.mark.parametrize("struct", test_cases, indirect=True)
    @pytest.mark.parametrize(
        "implem", [Solver.ImplemType.MAKUTU]
    )
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
            Solver.PhysicType.ACOUSTIC,
            sd.order,
        )

        n_nodes = sd.n_dof
        n_grad = sd.n_dof if is_model_on_nodes else sd.n_elements

        kk_pn, _ = Utils.allocate_field(n_nodes)
        kk_qn, _ = Utils.allocate_field(n_nodes)
        kk_qn_prev, _ = Utils.allocate_field(n_nodes)
        kk_qn_prev_prev, _ = Utils.allocate_field(n_nodes)
        kk_grad_kappa, _ = Utils.allocate_field(n_grad)
        kk_grad_buoyancy, _ = Utils.allocate_field(n_grad)

        dt = 0.001
        fwd = Gradient.WavefieldViewForwardAcoustic(kk_pn)
        bwd = Gradient.WavefieldViewBackwardAcoustic(kk_qn, kk_qn_prev, kk_qn_prev_prev)
        grad = Gradient.GradientAcoustic(kk_grad_kappa, kk_grad_buoyancy)
        data = Gradient.GradientDataAcoustic(fwd, bwd, grad)

        differentiator.compute(model, data, dt)
