import pyfuntides.model as Model
import pytest


class StructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 10
        self.domain_size = 1500
        self.hx = self.domain_size / self.ex
        self.hy = self.domain_size / self.ey
        self.hz = self.domain_size / self.ez
        self.order = order


@pytest.fixture
def struct(request):
    order, builder_cls, on_nodes, is_elastic = request.param

    sd = StructData(order)

    builder = builder_cls(sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, on_nodes, is_elastic)

    return sd, builder


test_cases = [
    # f32, i32 cases
    (1, Model.CartesianStructBuilder_f32_i32_O1, True, True),
    (1, Model.CartesianStructBuilder_f32_i32_O1, False, True),
    (1, Model.CartesianStructBuilder_f32_i32_O1, True, False),
    (1, Model.CartesianStructBuilder_f32_i32_O1, False, False),
    (2, Model.CartesianStructBuilder_f32_i32_O2, True, True),
    (2, Model.CartesianStructBuilder_f32_i32_O2, False, True),
    (2, Model.CartesianStructBuilder_f32_i32_O2, True, False),
    (2, Model.CartesianStructBuilder_f32_i32_O2, False, False),
    (3, Model.CartesianStructBuilder_f32_i32_O3, True, True),
    (3, Model.CartesianStructBuilder_f32_i32_O3, False, True),
    (3, Model.CartesianStructBuilder_f32_i32_O3, True, False),
    (3, Model.CartesianStructBuilder_f32_i32_O3, False, False),
    # f64, i32 cases
    (1, Model.CartesianStructBuilder_f64_i32_O1, True, True),
    (1, Model.CartesianStructBuilder_f64_i32_O1, False, True),
    (1, Model.CartesianStructBuilder_f64_i32_O1, True, False),
    (1, Model.CartesianStructBuilder_f64_i32_O1, False, False),
    (2, Model.CartesianStructBuilder_f64_i32_O2, True, True),
    (2, Model.CartesianStructBuilder_f64_i32_O2, False, True),
    (2, Model.CartesianStructBuilder_f64_i32_O2, True, False),
    (2, Model.CartesianStructBuilder_f64_i32_O2, False, False),
    (3, Model.CartesianStructBuilder_f64_i32_O3, True, True),
    (3, Model.CartesianStructBuilder_f64_i32_O3, False, True),
    (3, Model.CartesianStructBuilder_f64_i32_O3, True, False),
    (3, Model.CartesianStructBuilder_f64_i32_O3, False, False),
    # f32, i64 cases
    (1, Model.CartesianStructBuilder_f32_i64_O1, True, True),
    (1, Model.CartesianStructBuilder_f32_i64_O1, False, True),
    (1, Model.CartesianStructBuilder_f32_i64_O1, True, False),
    (1, Model.CartesianStructBuilder_f32_i64_O1, False, False),
    (2, Model.CartesianStructBuilder_f32_i64_O2, True, True),
    (2, Model.CartesianStructBuilder_f32_i64_O2, False, True),
    (2, Model.CartesianStructBuilder_f32_i64_O2, True, False),
    (2, Model.CartesianStructBuilder_f32_i64_O2, False, False),
    (3, Model.CartesianStructBuilder_f32_i64_O3, True, True),
    (3, Model.CartesianStructBuilder_f32_i64_O3, False, True),
    (3, Model.CartesianStructBuilder_f32_i64_O3, True, False),
    (3, Model.CartesianStructBuilder_f32_i64_O3, False, False),
    # f64, i64 cases
    (1, Model.CartesianStructBuilder_f64_i64_O1, True, True),
    (1, Model.CartesianStructBuilder_f64_i64_O1, False, True),
    (1, Model.CartesianStructBuilder_f64_i64_O1, True, False),
    (1, Model.CartesianStructBuilder_f64_i64_O1, False, False),
    (2, Model.CartesianStructBuilder_f64_i64_O2, True, True),
    (2, Model.CartesianStructBuilder_f64_i64_O2, False, True),
    (2, Model.CartesianStructBuilder_f64_i64_O2, True, False),
    (2, Model.CartesianStructBuilder_f64_i64_O2, False, False),
    (3, Model.CartesianStructBuilder_f64_i64_O3, True, True),
    (3, Model.CartesianStructBuilder_f64_i64_O3, False, True),
    (3, Model.CartesianStructBuilder_f64_i64_O3, True, False),
    (3, Model.CartesianStructBuilder_f64_i64_O3, False, False),
]


class TestCartesianStructBuilder:
    @pytest.mark.parametrize("struct", test_cases, indirect=True)
    def test_cartesian_struct_builder(self, struct):
        _, builder = struct
        model = builder.get_model()
        assert model is not None
