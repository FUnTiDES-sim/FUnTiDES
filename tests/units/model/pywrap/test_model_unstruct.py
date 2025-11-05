import kokkos
import numpy as np
import pyfuntides.model as Model
import pytest


class UnstructData:
    """
    A cube divided into 4 elements
    +----+----+
    | E0 | E1 |
    +----+----+
    | E2 | E3 |
    +----+----+
    """

    def __init__(self, order):
        self.lx = self.ly = self.lz = 1500
        self.order = order
        self.n_elements_per_dim = 2
        self.n_elements = self.n_elements_per_dim**3
        self.n_nodes = (self.n_elements_per_dim * order + 1) ** 3
        self.memspace = kokkos.DefaultMemorySpace
        self.layout = kokkos.LayoutRight
        # TODO so far model only accepts f32, i32 for kokkos arrays
        float_type = kokkos.float32
        scalar_type = kokkos.int32
        # coords
        self.kk_nodes_coords_x = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_nodes_coords_y = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_nodes_coords_z = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_global_node_index = kokkos.array(
            [self.n_elements, self.n_nodes],
            dtype=scalar_type,
            space=self.memspace,
            layout=self.layout,
        )
        # model
        self.kk_model_vp_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_vp_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_rho_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_rho_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_vs_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_vs_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_delta_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_delta_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_epsilon_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_epsilon_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_gamma_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_gamma_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_theta_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_theta_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_phi_node = kokkos.array(
            [self.n_nodes], dtype=float_type, space=self.memspace, layout=self.layout
        )
        self.kk_model_phi_element = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )

        # boundaries
        self.kk_boundaries = kokkos.array(
            [self.n_elements], dtype=float_type, space=self.memspace, layout=self.layout
        )

    def generate_global_coordinates(self):
        n_points = 2 * self.order + 1
        i, j, k = np.meshgrid(
            np.arange(n_points), np.arange(n_points), np.arange(n_points), indexing="ij"
        )
        i = i.ravel()
        j = j.ravel()
        k = k.ravel()

        dx = self.lx / (n_points - 1)
        dy = self.ly / (n_points - 1)
        dz = self.lz / (n_points - 1)

        nodes_coords_x = np.array(self.kk_nodes_coords_x, copy=False)
        nodes_coords_y = np.array(self.kk_nodes_coords_y, copy=False)
        nodes_coords_z = np.array(self.kk_nodes_coords_z, copy=False)

        for idx in range(self.n_nodes):
            nodes_coords_x[idx] = i[idx] * dx
            nodes_coords_y[idx] = j[idx] * dy
            nodes_coords_z[idx] = k[idx] * dz

    def generate_global_node_index_map(self):
        n_points = self.order + 1
        ex, ey, ez = np.meshgrid(
            np.arange(self.n_elements_per_dim),
            np.arange(self.n_elements_per_dim),
            np.arange(self.n_elements_per_dim),
            indexing="ij",
        )

        ex = ex.ravel()
        ey = ey.ravel()
        ez = ez.ravel()

        global_node_index = np.array(self.kk_global_node_index, copy=False)

        for e in range(self.n_elements):
            for i in range(n_points):
                for j in range(n_points):
                    for k in range(n_points):
                        local_index = i + j * n_points + k * n_points * n_points
                        global_index = (
                            (ex[e] * self.order + i)
                            + (ey[e] * self.order + j) * (2 * self.order + 1)
                            + (ez[e] * self.order + k) * (2 * self.order + 1) ** 2
                        )
                        global_node_index[e, local_index] = global_index

    def fill_model(self):
        model_vp_node = np.array(self.kk_model_vp_node, copy=False)
        model_vp_element = np.array(self.kk_model_vp_element, copy=False)
        model_rho_node = np.array(self.kk_model_rho_node, copy=False)
        model_rho_element = np.array(self.kk_model_rho_element, copy=False)
        model_vs_node = np.array(self.kk_model_vs_node, copy=False)
        model_vs_element = np.array(self.kk_model_vs_element, copy=False)
        model_delta_node = np.array(self.kk_model_delta_node, copy=False)
        model_delta_element = np.array(self.kk_model_delta_element, copy=False)
        model_epsilon_node = np.array(self.kk_model_epsilon_node, copy=False)
        model_epsilon_element = np.array(self.kk_model_epsilon_element, copy=False)
        model_gamma_node = np.array(self.kk_model_gamma_node, copy=False)
        model_gamma_element = np.array(self.kk_model_gamma_element, copy=False)
        model_phi_node = np.array(self.kk_model_phi_node, copy=False)
        model_phi_element = np.array(self.kk_model_phi_element, copy=False)
        model_theta_node = np.array(self.kk_model_theta_node, copy=False)
        model_theta_element = np.array(self.kk_model_theta_element, copy=False)
        model_vp_node[:] = 3000.0
        model_vp_element[:] = 3000.0
        # set different mid value for min/max testing
        model_vp_node[self.n_nodes // 2] = 3600.0
        model_vp_element[self.n_elements // 2] = 3500.0
        model_rho_node[:] = 2500.0
        model_rho_element[:] = 2500.0
        model_vs_node[:] = 755.0
        model_vs_element[:] = 755.0
        model_delta_node[:] = 0.1
        model_delta_element[:] = 0.1
        model_epsilon_node[:] = 0.2
        model_epsilon_element[:] = 0.2
        model_gamma_node[:] = 0.15
        model_gamma_element[:] = 0.15
        model_theta_node[:] = 30.0
        model_theta_element[:] = 30.0
        model_phi_node[:] = 45.0
        model_phi_element[:] = 45.0

    def fill_boundaries(self):
        boundaries = np.array(self.kk_boundaries, copy=False)
        boundaries[:] = 1.0


@pytest.fixture
def unstruct(request):
    order, param_cls, model_cls, on_nodes, is_elastic = request.param

    ud = UnstructData(order)
    ud.generate_global_coordinates()
    ud.generate_global_node_index_map()
    ud.fill_model()
    ud.fill_boundaries()

    params = param_cls(
        ud.order,
        ud.n_elements,
        ud.n_nodes,
        ud.lx,
        ud.ly,
        ud.lz,
        on_nodes,
        is_elastic,
        ud.kk_global_node_index,
        ud.kk_nodes_coords_x,
        ud.kk_nodes_coords_y,
        ud.kk_nodes_coords_z,
        ud.kk_model_vp_node,
        ud.kk_model_vp_element,
        ud.kk_model_rho_node,
        ud.kk_model_rho_element,
        ud.kk_model_vs_node,
        ud.kk_model_vs_element,
        ud.kk_model_delta_node,
        ud.kk_model_delta_element,
        ud.kk_model_epsilon_node,
        ud.kk_model_epsilon_element,
        ud.kk_model_gamma_node,
        ud.kk_model_gamma_element,
        ud.kk_model_theta_node,
        ud.kk_model_theta_element,
        ud.kk_model_phi_node,
        ud.kk_model_phi_element,
        ud.kk_boundaries,
    )

    return ud, model_cls, params


test_cases = [
    # f32, i32 cases
    (1, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, True, True),
    (1, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, True, False),
    (1, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, False, True),
    (1, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, False, False),
    (2, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, True, True),
    (2, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, True, False),
    (2, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, False, True),
    (2, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, False, False),
    (3, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, True, True),
    (3, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, True, False),
    (3, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, False, True),
    (3, Model.ModelUnstructData_f32_i32, Model.ModelUnstruct_f32_i32, False, False),
    # f64, i32 cases
    (1, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, True, True),
    (1, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, True, False),
    (1, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, False, True),
    (1, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, False, False),
    (2, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, True, True),
    (2, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, True, False),
    (2, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, False, True),
    (2, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, False, False),
    (3, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, True, True),
    (3, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, True, False),
    (3, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, False, True),
    (3, Model.ModelUnstructData_f64_i32, Model.ModelUnstruct_f64_i32, False, False),
    ## f32, i64 cases
    (1, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, True, True),
    (1, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, True, False),
    (1, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, False, True),
    (1, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, False, False),
    (2, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, True, True),
    (2, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, True, False),
    (2, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, False, True),
    (2, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, False, False),
    (3, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, True, True),
    (3, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, True, False),
    (3, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, False, True),
    (3, Model.ModelUnstructData_f32_i64, Model.ModelUnstruct_f32_i64, False, False),
    ## f64, i64 cases
    (1, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, True, True),
    (1, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, True, False),
    (1, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, False, True),
    (1, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, False, False),
    (2, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, True, True),
    (2, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, True, False),
    (2, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, False, True),
    (2, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, False, False),
    (3, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, True, True),
    (3, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, True, False),
    (3, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, False, True),
    (3, Model.ModelUnstructData_f64_i64, Model.ModelUnstruct_f64_i64, False, False),
]


class TestModelUnstruct:
    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_model_unstruct(self, unstruct):
        _, model_cls, params = unstruct
        model = model_cls(params)
        assert model is not None

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_node_coord(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        # Test coordinates for first node
        assert model.node_coord(0, 0) == 0.0
        assert model.node_coord(0, 1) == 0.0
        assert model.node_coord(0, 2) == 0.0

        # Test coordinates for middle node
        mid_node_index = data.n_nodes // 2
        assert model.node_coord(mid_node_index, 0) == 750.0
        assert model.node_coord(mid_node_index, 1) == 750.0
        assert model.node_coord(mid_node_index, 2) == 750.0

        # Test coordinates for last node
        last_node_index = data.n_nodes - 1
        assert model.node_coord(last_node_index, 0) == 1500.0
        assert model.node_coord(last_node_index, 1) == 1500.0
        assert model.node_coord(last_node_index, 2) == 1500.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_global_node_index(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        stride = 2 * data.order + 1

        # Test first element node indices
        assert model.global_node_index(0, 0, 0, 0) == 0
        assert model.global_node_index(0, 1, 0, 0) == 1
        assert model.global_node_index(0, 0, 1, 0) == stride

        # If order is 2 or higher, test mid-node indices
        if data.order >= 2:
            expected_mid = 1 + 1 * stride + 1 * stride * stride
            expected_corner = 2 + 2 * stride + 2 * stride * stride
            assert model.global_node_index(0, 1, 1, 1) == expected_mid
            assert model.global_node_index(0, 2, 2, 2) == expected_corner

        # Test last element node indices
        last_elem = data.n_elements - 1
        offset = data.order
        assert (
            model.global_node_index(last_elem, 0, 0, 0)
            == offset + offset * stride + offset * stride * stride
        )
        assert (
            model.global_node_index(last_elem, data.order, data.order, data.order)
            == data.n_nodes - 1
        )

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_vp_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            if n == data.n_nodes // 2:
                assert model.get_model_vp_on_node(n) == 3600.0
            else:
                assert model.get_model_vp_on_node(n) == 3000.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_vp_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            if e == data.n_elements // 2:
                assert model.get_model_vp_on_element(e) == 3500.0
            else:
                assert model.get_model_vp_on_element(e) == 3000.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_rho_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            assert model.get_model_rho_on_node(n) == 2500.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_rho_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            assert model.get_model_rho_on_element(e) == 2500.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_vs_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            assert model.get_model_vs_on_node(n) == 755.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_vs_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            assert model.get_model_vs_on_element(e) == 755.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_delta_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            assert model.get_model_delta_on_node(n) == pytest.approx(0.1)

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_delta_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            assert model.get_model_delta_on_element(e) == pytest.approx(0.1)

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_epsilon_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            assert model.get_model_epsilon_on_node(n) == pytest.approx(0.2)

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_epsilon_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            assert model.get_model_epsilon_on_element(e) == pytest.approx(0.2)

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_gamma_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            assert model.get_model_gamma_on_node(n) == pytest.approx(0.15)

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_gamma_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            assert model.get_model_gamma_on_element(e) == pytest.approx(0.15)

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_phi_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            assert model.get_model_phi_on_node(n) == 45.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_phi_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            assert model.get_model_phi_on_element(e) == 45.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_theta_on_node(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for n in range(data.n_nodes):
            assert model.get_model_theta_on_node(n) == 30.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_model_theta_on_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        for e in range(data.n_elements):
            assert model.get_model_theta_on_element(e) == 30.0

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_number_of_elements(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        assert model.get_number_of_elements() == data.n_elements

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_number_of_nodes(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        assert model.get_number_of_nodes() == data.n_nodes

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_number_of_points_per_element(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        assert model.get_number_of_points_per_element() == (data.order + 1) ** 3

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_order(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        assert model.get_order() == data.order

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_domain_size(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        assert model.domain_size(0) == data.lx
        assert model.domain_size(1) == data.ly
        assert model.domain_size(2) == data.lz

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_min_spacing(self, unstruct):
        data, model_cls, params = unstruct
        model = model_cls(params)

        n_points = 2 * data.order + 1
        dx = data.lx / (n_points - 1)
        dy = data.ly / (n_points - 1)
        dz = data.lz / (n_points - 1)
        min_d = min(dx, dy, dz)

        assert model.get_min_spacing() == min_d

    @pytest.mark.parametrize("unstruct", test_cases, indirect=True)
    def test_get_max_speed(self, unstruct):
        _, model_cls, params = unstruct
        model = model_cls(params)

        assert model.get_max_speed() == 3600.0
