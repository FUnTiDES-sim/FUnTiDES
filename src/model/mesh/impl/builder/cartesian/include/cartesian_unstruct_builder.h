#ifndef SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_
#define SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_

#include <builder.h>
#include <model_unstruct.h>

#include "cartesian_params.h"
#include "sep_params.h"

namespace model
{
template <typename FloatType, typename ScalarType>
class CartesianUnstructBuilder : public ModelBuilderBase<FloatType, ScalarType>
{
 public:
  using ModelBuilderBase<FloatType, ScalarType>::MAX_ORDER;

  CartesianUnstructBuilder() {}

  CartesianUnstructBuilder(const CartesianParams<FloatType, ScalarType>& p)
      : ex_(p.ex),
        ey_(p.ey),
        ez_(p.ez),
        lx_(p.lx),
        ly_(p.ly),
        lz_(p.lz),
        order_(p.order),
        isModelOnNodes_(p.isModelOnNodes),
        isElastic_(p.isElastic)

  {
    initGlobalNodeList();
    initNodesCoords();
    initModelsUniform();
  }

  CartesianUnstructBuilder(const SepParams<FloatType, ScalarType>& p,
                           const int order)
      : ex_(p.getN1()), ey_(p.getN2()), ez_(p.getN3()), order_(order)
  {
    lx_ = p.getD1() * ex_;
    ly_ = p.getD2() * ey_;
    lz_ = p.getD3() * ez_;

    isModelOnNodes_ = false;
    isElastic_ = false;

    initGlobalNodeList();
    initNodesCoords();
    initModelsFromSep(p);
    // printModelsToFiles();
  }

  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel()
      const override
  {
    model::ModelUnstructData<FloatType, ScalarType> modelData;

    modelData.order_ = order_;
    modelData.n_element_ = ex_ * ey_ * ez_;
    modelData.n_node_ =
        (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);
    modelData.lx_ = lx_;
    modelData.ly_ = ly_;
    modelData.lz_ = lz_;

    modelData.global_node_index_ = global_node_index_;
    modelData.nodes_coords_x_ = nodes_coords_x_;
    modelData.nodes_coords_y_ = nodes_coords_y_;
    modelData.nodes_coords_z_ = nodes_coords_z_;

    modelData.isModelOnNodes_ = isModelOnNodes_;
    modelData.isElastic_ = isElastic_;

    modelData.model_vp_node_ = model_vp_node_;
    modelData.model_rho_node_ = model_rho_node_;
    modelData.model_vs_node_ = model_vs_node_;
    modelData.model_delta_node_ = model_delta_node_;
    modelData.model_epsilon_node_ = model_epsilon_node_;
    modelData.model_gamma_node_ = model_gamma_node_;
    modelData.model_theta_node_ = model_theta_node_;
    modelData.model_phi_node_ = model_phi_node_;
    modelData.model_vp_element_ = model_vp_element_;
    modelData.model_rho_element_ = model_rho_element_;
    modelData.model_vs_element_ = model_vs_element_;
    modelData.model_delta_element_ = model_delta_element_;
    modelData.model_epsilon_element_ = model_epsilon_element_;
    modelData.model_gamma_element_ = model_gamma_element_;
    modelData.model_theta_element_ = model_theta_element_;
    modelData.model_phi_element_ = model_phi_element_;

    auto model = std::make_shared<model::ModelUnstruct<FloatType, ScalarType>>(
        modelData);

    if (isElastic_ && !isModelOnNodes_)
    {
      model->initElasticityTensors();
    }

    return model;
  }

  ~CartesianUnstructBuilder() = default;

 private:
  ScalarType ex_, ey_, ez_;
  FloatType lx_, ly_, lz_;
  int order_;
  bool isModelOnNodes_;
  bool isElastic_;

  ARRAY_INT_VIEW global_node_index_;
  VECTOR_REAL_VIEW nodes_coords_x_;
  VECTOR_REAL_VIEW nodes_coords_y_;
  VECTOR_REAL_VIEW nodes_coords_z_;
  // Models view
  VECTOR_REAL_VIEW model_vp_node_;
  VECTOR_REAL_VIEW model_vp_element_;
  VECTOR_REAL_VIEW model_rho_node_;
  VECTOR_REAL_VIEW model_rho_element_;
  VECTOR_REAL_VIEW model_vs_node_;
  VECTOR_REAL_VIEW model_vs_element_;
  VECTOR_REAL_VIEW model_delta_node_;
  VECTOR_REAL_VIEW model_delta_element_;
  VECTOR_REAL_VIEW model_epsilon_node_;
  VECTOR_REAL_VIEW model_epsilon_element_;
  VECTOR_REAL_VIEW model_gamma_node_;
  VECTOR_REAL_VIEW model_gamma_element_;
  VECTOR_REAL_VIEW model_theta_node_;
  VECTOR_REAL_VIEW model_theta_element_;
  VECTOR_REAL_VIEW model_phi_node_;
  VECTOR_REAL_VIEW model_phi_element_;

  VECTOR_REAL_VIEW boundaries_t_;

  void initGlobalNodeList()
  {
    int nodes_x = order_ + 1;
    int nodes_y = order_ + 1;
    int nodes_z = order_ + 1;
    int total_nodes = nodes_x * nodes_y * nodes_z;
    global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(
        ex_ * ey_ * ez_, total_nodes, "global node index");
    int nx = ex_ * order_ + 1;  // Total nodes in x direction
    int ny = ey_ * order_ + 1;  // Total nodes in y direction
    int nz = ez_ * order_ + 1;  // Total nodes in z direction

    for (int k = 0; k < ez_; k++)
    {
      for (int j = 0; j < ey_; j++)
      {
        for (int i = 0; i < ex_; i++)
        {
          int elementNum = i + j * ex_ + k * ex_ * ey_;
          // Corrected offset calculation
          int offset = i * order_ + j * order_ * nx + k * order_ * nx * ny;

          for (int m = 0; m < order_ + 1; m++)
          {  // z-direction
            for (int n = 0; n < order_ + 1; n++)
            {  // y-direction
              for (int l = 0; l < order_ + 1; l++)
              {  // x-direction
                int dofLocal =
                    l + n * (order_ + 1) + m * (order_ + 1) * (order_ + 1);
                int dofGlobal = offset + l + n * nx + m * nx * ny;
                global_node_index_(elementNum, dofLocal) = dofGlobal;
              }
            }
          }
        }
      }
    }
  }

  void getCoordInOneDirection(const int& h, const int& n_element, float* coord)
  {
    float xi[MAX_ORDER + 1];

    switch (order_)
    {
      case 1:
        xi[0] = -1.f;
        xi[1] = 1.f;
        break;
      case 2:
        xi[0] = -1.f;
        xi[1] = 0.f;
        xi[2] = 1.f;
        break;
      case 3: {
        static constexpr float sqrt5 = 2.2360679774997897f;
        xi[0] = -1.0f;
        xi[1] = -1.f / sqrt5;
        xi[2] = 1.f / sqrt5;
        xi[3] = 1.f;
        break;
      }
      case 4: {
        static constexpr float sqrt3_7 = 0.6546536707079771f;
        xi[0] = -1.0f;
        xi[1] = -sqrt3_7;
        xi[2] = 0.0f;
        xi[3] = sqrt3_7;
        xi[4] = 1.0f;
        break;
      }
      case 5: {
        static constexpr float sqrt__7_plus_2sqrt7__ = 3.50592393273573196f;
        static constexpr float sqrt__7_mins_2sqrt7__ = 1.30709501485960033f;
        static constexpr float sqrt_inv21 = 0.218217890235992381f;
        xi[0] = -1.0f;
        xi[1] = -sqrt_inv21 * sqrt__7_plus_2sqrt7__;
        xi[2] = -sqrt_inv21 * sqrt__7_mins_2sqrt7__;
        xi[3] = sqrt_inv21 * sqrt__7_mins_2sqrt7__;
        xi[4] = sqrt_inv21 * sqrt__7_plus_2sqrt7__;
        xi[5] = 1.0f;
        break;
      }
      default:
        break;
    }

    int i = n_element;
    float x0 = i * h;
    float x1 = (i + 1) * h;
    float b = (x1 + x0) / 2.f;
    float a = b - x0;

    for (int j = 0; j < order_ + 1; j++)
    {
      coord[j] = a * xi[j] + b;
    }
  }

  void initNodesCoords()
  {
    int nodes_x = ex_ * order_ + 1;
    int nodes_y = ey_ * order_ + 1;
    int nodes_z = ez_ * order_ + 1;
    int total_nodes = nodes_x * nodes_y * nodes_z;

    // Init the structure within mesh
    nodes_coords_x_ =
        allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords x");
    nodes_coords_y_ =
        allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords y");
    nodes_coords_z_ =
        allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords z");

    float coord_x[MAX_ORDER + 1];
    float coord_y[MAX_ORDER + 1];
    float coord_z[MAX_ORDER + 1];

    auto hx = lx_ / ex_;
    auto hy = ly_ / ey_;
    auto hz = lz_ / ez_;

    for (int n = 0; n < ez_; n++)
    {
      getCoordInOneDirection(hz, n, coord_z);
      for (int m = 0; m < ey_; m++)
      {
        getCoordInOneDirection(hy, m, coord_y);
        for (int l = 0; l < ex_; l++)
        {
          getCoordInOneDirection(hx, l, coord_x);

          for (int k = 0; k < order_ + 1; k++)
          {
            for (int j = 0; j < order_ + 1; j++)
            {
              for (int i = 0; i < order_ + 1; i++)
              {
                int global_i = l * order_ + i;
                int global_j = m * order_ + j;
                int global_k = n * order_ + k;

                int global_node_index = global_i + global_j * nodes_x +
                                        global_k * nodes_x * nodes_y;

                if (global_i < nodes_x && global_j < nodes_y &&
                    global_k < nodes_z)
                {
                  nodes_coords_x_(global_node_index) = coord_x[i];
                  nodes_coords_y_(global_node_index) = coord_y[j];
                  nodes_coords_z_(global_node_index) = coord_z[k];
                }
              }
            }
          }
        }
      }
    }
  }

  /**
   * @brief Copy data from std::vector to VECTOR_REAL_VIEW
   *
   * @param src Source vector (from SEP file)
   * @param dest Destination view (already allocated)
   * @param count Number of elements to copy
   */
  void copyToView(const std::vector<float>& src, VECTOR_REAL_VIEW& dest,
                  size_t count)
  {
    for (size_t i = 0; i < count; ++i)
    {
      dest(i) = src[i];
    }
  }

  /**
   * @brief Initialize models from SEP parameter file
   *
   * Reads binary velocity data from the SEP file and stores it in
   * model_vp_element_ or model_vp_node_ depending on isModelOnNodes_.
   *
   * @param p SEP parameters containing file path and dimensions
   */
  void initModelsFromSep(const SepParams<FloatType, ScalarType>& p)
  {
    const int n_element = ex_ * ey_ * ez_;

    // Read velocity model from SEP file
    std::vector<float> vp_data = p.template readBinaryData<float>();

    if (isModelOnNodes_)
    {
      throw std::runtime_error(
          "[Cartesian Builder] Cannot init model on nodes from sep files.");
    }
    else
    {
      if (static_cast<size_t>(n_element) != vp_data.size())
      {
        throw std::runtime_error("SEP data size (" +
                                 std::to_string(vp_data.size()) +
                                 ") does not match expected element count (" +
                                 std::to_string(n_element) + ")");
      }

      model_vp_element_ =
          allocateVector<VECTOR_REAL_VIEW>(n_element, "model vp elem");
      model_rho_element_ =
          allocateVector<VECTOR_REAL_VIEW>(n_element, "model rho elem");

      copyToView(vp_data, model_vp_element_, n_element);

      // Initialize rho with uniform value
      for (int i = 0; i < n_element; i++)
      {
        model_rho_element_(i) = 1.0f;
      }

      if (isElastic_)
      {
        model_vs_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model vs element");
        model_delta_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model delta element");
        model_gamma_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model gamma element");
        model_epsilon_element_ = allocateVector<VECTOR_REAL_VIEW>(
            n_element, "model epsilon element");
        model_theta_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model theta element");
        model_phi_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model phi element");

        for (int i = 0; i < n_element; i++)
        {
          model_vs_element_(i) = 755.0f;
          model_delta_element_(i) = 0.2f;
          model_epsilon_element_(i) = 0.3f;
          model_gamma_element_(i) = 0.08f;
          model_theta_element_(i) = 30.0f;
          model_phi_element_(i) = 45.0f;
        }
      }
    }
  }

  /**
   * @brief Initialize models with uniform values
   *
   * Creates velocity and density models with constant values throughout
   * the domain. Used when no external model file is provided.
   */
  void initModelsUniform()
  {
    int n_element = ex_ * ey_ * ez_;
    int n_node = (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);

    if (isModelOnNodes_)
    {
      model_rho_node_ =
          allocateVector<VECTOR_REAL_VIEW>(n_node, "model rho node");
      model_vp_node_ =
          allocateVector<VECTOR_REAL_VIEW>(n_node, "model vp node");

      for (int i = 0; i < n_node; i++)
      {
        model_rho_node_(i) = 1.0f;
        model_vp_node_(i) = 1500.0f;
      }

      if (isElastic_)
      {
        model_vs_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model vs node");
        model_delta_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model delta node");
        model_gamma_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model gamma node");
        model_epsilon_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model epsilon node");
        model_theta_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model theta node");
        model_phi_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model phi node");

        for (int i = 0; i < n_node; i++)
        {
          model_vs_node_(i) = 755.0f;
          model_delta_node_(i) = 0.2f;
          model_epsilon_node_(i) = 0.3f;
          model_gamma_node_(i) = 0.08f;
          model_theta_node_(i) = 30.0f;
          model_phi_node_(i) = 45.0f;
        }
      }
    }
    else
    {
      model_rho_element_ =
          allocateVector<VECTOR_REAL_VIEW>(n_element, "model rho elem");
      model_vp_element_ =
          allocateVector<VECTOR_REAL_VIEW>(n_element, "model vp elem");

      for (int i = 0; i < n_element; i++)
      {
        model_rho_element_(i) = 1.0f;
        model_vp_element_(i) = 1500.0f;
      }

      if (isElastic_)
      {
        model_vs_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model vs element");
        model_delta_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model delta element");
        model_gamma_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model gamma element");
        model_epsilon_element_ = allocateVector<VECTOR_REAL_VIEW>(
            n_element, "model epsilon element");
        model_theta_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model theta element");
        model_phi_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_element, "model phi element");

        for (int i = 0; i < n_element; i++)
        {
          model_vs_element_(i) = 755.0f;
          model_delta_element_(i) = 0.2f;
          model_epsilon_element_(i) = 0.3f;
          model_gamma_element_(i) = 0.08f;
          model_theta_element_(i) = 30.0f;
          model_phi_element_(i) = 45.0f;
        }
      }
    }
  }

  void printModelsToFiles() const
  {
    // Print vp_element
    if (model_vp_element_.extent(0) > 0)
    {
      auto host_vp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                         model_vp_element_);
      std::ofstream file_vp_elem("model_vp_element.txt");
      for (size_t i = 0; i < host_vp.extent(0); i++)
      {
        file_vp_elem << i << " " << host_vp(i) << "\n";
      }
      file_vp_elem.close();
    }

    // Print vp_node
    if (model_vp_node_.extent(0) > 0)
    {
      auto host_vp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                         model_vp_node_);
      std::ofstream file_vp_node("model_vp_node.txt");
      for (size_t i = 0; i < host_vp.extent(0); i++)
      {
        file_vp_node << i << " " << host_vp(i) << "\n";
      }
      file_vp_node.close();
    }

    // Print rho_element
    if (model_rho_element_.extent(0) > 0)
    {
      auto host_rho = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                          model_rho_element_);
      std::ofstream file_rho_elem("model_rho_element.txt");
      for (size_t i = 0; i < host_rho.extent(0); i++)
      {
        file_rho_elem << i << " " << host_rho(i) << "\n";
      }
      file_rho_elem.close();
    }

    // Print rho_node
    if (model_rho_node_.extent(0) > 0)
    {
      auto host_rho = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                          model_rho_node_);
      std::ofstream file_rho_node("model_rho_node.txt");
      for (size_t i = 0; i < host_rho.extent(0); i++)
      {
        file_rho_node << i << " " << host_rho(i) << "\n";
      }
      file_rho_node.close();
    }

    // Print nodes coordinates
    if (nodes_coords_x_.extent(0) > 0)
    {
      auto host_x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                        nodes_coords_x_);
      auto host_y = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                        nodes_coords_y_);
      auto host_z = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                        nodes_coords_z_);

      std::ofstream file_coords("nodes_coords.txt");
      for (size_t i = 0; i < host_x.extent(0); i++)
      {
        file_coords << i << " " << host_x(i) << " " << host_y(i) << " "
                    << host_z(i) << "\n";
      }
      file_coords.close();
    }
  }
};

}  // namespace model
#endif  // SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_
