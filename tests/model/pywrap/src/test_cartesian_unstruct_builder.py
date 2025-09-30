import pytest
import pyproxys.model as Model


class UnstructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 10
        self.lx = self.ly = self.lz = 1500
        self.order = order


@pytest.fixture
def unstruct(request):
    order, param_cls, builder_cls = request.param

    sd = UnstructData(order)

    params = param_cls()
    params.ex, params.ey, params.ez = sd.ex, sd.ey, sd.ez
    params.lx, params.ly, params.lz = sd.lx, sd.ly, sd.lz
    params.order = order

    builder = builder_cls(params)

    return sd, params, builder


test_cases = [
    # f32, i32 cases
    (1, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32),
    (2, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32),
    (3, Model.CartesianParams_f32_i32, Model.CartesianUnstructBuilder_f32_i32),
    # f64, i32 cases
    (1, Model.CartesianParams_f64_i32, Model.CartesianUnstructBuilder_f64_i32),
    (2, Model.CartesianParams_f64_i32, Model.CartesianUnstructBuilder_f64_i32),
    (3, Model.CartesianParams_f64_i32, Model.CartesianUnstructBuilder_f64_i32),
    # f32, i64 cases
    (1, Model.CartesianParams_f32_i64, Model.CartesianUnstructBuilder_f32_i64),
    (2, Model.CartesianParams_f32_i64, Model.CartesianUnstructBuilder_f32_i64),
    (3, Model.CartesianParams_f32_i64, Model.CartesianUnstructBuilder_f32_i64),
    # f64, i64 cases
    (1, Model.CartesianParams_f64_i64, Model.CartesianUnstructBuilder_f64_i64),
    (2, Model.CartesianParams_f64_i64, Model.CartesianUnstructBuilder_f64_i64),
    (3, Model.CartesianParams_f64_i64, Model.CartesianUnstructBuilder_f64_i64),
]


class TestCartesianUnstructBuilder:
    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_cartesian_unstruct_builder(self, unstruct):
        _, _, builder = unstruct
        model = builder.get_model()
        assert model is not None
