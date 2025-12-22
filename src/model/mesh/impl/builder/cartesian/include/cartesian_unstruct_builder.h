#ifndef SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_
#define SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_

#include <builder.h>
#include <model_unstruct.h>

#include "cartesian_params.h"

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
    if (isModelOnNodes_)
    {
      // High-Order Mode: Connectivity contains all GLL points
      initGlobalNodeList();
      initNodesCoords();
    }
    else
    {
      // Optimized Mode: Connectivity contains only 8 corners
      initGlobalCornerNodeList();
      initCornerNodesCoords();
    }
    initModels();
  }

  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel()
      const override
  {
    model::ModelUnstructData<FloatType, ScalarType> modelData;

    // Mathematical GLL order remains constant
    modelData.order_ = order_;
    modelData.n_element_ = ex_ * ey_ * ez_;

    // n_node_ represents the actual number of coordinates/properties allocated
    if (isModelOnNodes_)
    {
      modelData.n_node_ =
          (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);
    }
    else
    {
      modelData.n_node_ = (ex_ + 1) * (ey_ + 1) * (ez_ + 1);
    }

    modelData.lx_ = lx_;
    modelData.ly_ = ly_;
    modelData.lz_ = lz_;

    modelData.global_node_index_ = global_node_index_;
    modelData.nodes_coords_x_ = nodes_coords_x_;
    modelData.nodes_coords_y_ = nodes_coords_y_;
    modelData.nodes_coords_z_ = nodes_coords_z_;

    modelData.isModelOnNodes_ = isModelOnNodes_;
    modelData.isElastic_ = isElastic_;

    // Map property views
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

 private:
  ScalarType ex_, ey_, ez_;
  FloatType lx_, ly_, lz_;
  int order_;
  bool isModelOnNodes_;
  bool isElastic_;

  ARRAY_INT_VIEW global_node_index_;
  VECTOR_REAL_VIEW nodes_coords_x_, nodes_coords_y_, nodes_coords_z_;
  VECTOR_REAL_VIEW model_vp_node_, model_vp_element_;
  VECTOR_REAL_VIEW model_rho_node_, model_rho_element_;
  VECTOR_REAL_VIEW model_vs_node_, model_vs_element_;
  VECTOR_REAL_VIEW model_delta_node_, model_delta_element_;
  VECTOR_REAL_VIEW model_epsilon_node_, model_epsilon_element_;
  VECTOR_REAL_VIEW model_gamma_node_, model_gamma_element_;
  VECTOR_REAL_VIEW model_theta_node_, model_theta_element_;
  VECTOR_REAL_VIEW model_phi_node_, model_phi_element_;

  // ==========================================================================
  // OPTIMIZED CORNER-ONLY INITIALIZATION
  // ==========================================================================
  void initGlobalCornerNodeList()
  {
    global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(ex_ * ey_ * ez_, 8,
                                                         "global corner index");

    int nx = ex_ + 1;
    int ny = ey_ + 1;

    for (int k = 0; k < ez_; k++)
    {
      for (int j = 0; j < ey_; j++)
      {
        for (int i = 0; i < ex_; i++)
        {
          int elementNum = i + j * ex_ + k * ex_ * ey_;
          for (int m = 0; m < 2; m++)
          {
            for (int n = 0; n < 2; n++)
            {
              for (int l = 0; l < 2; l++)
              {
                int localIdx = l + n * 2 + m * 4;
                int globalIdx = (i + l) + (j + n) * nx + (k + m) * nx * ny;
                global_node_index_(elementNum, localIdx) = globalIdx;
              }
            }
          }
        }
      }
    }
  }

  void initCornerNodesCoords()
  {
    int nx = ex_ + 1;
    int ny = ey_ + 1;
    int nz = ez_ + 1;
    int total_nodes = nx * ny * nz;

    nodes_coords_x_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "coords x");
    nodes_coords_y_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "coords y");
    nodes_coords_z_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "coords z");

    auto dx = lx_ / static_cast<double>(ex_);
    auto dy = ly_ / static_cast<double>(ey_);
    auto dz = lz_ / static_cast<double>(ez_);

    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int i = 0; i < nx; i++)
        {
          int idx = i + j * nx + k * nx * ny;
          nodes_coords_x_(idx) = i * dx;
          nodes_coords_y_(idx) = j * dy;
          nodes_coords_z_(idx) = k * dz;
        }
      }
    }
  }

  // ==========================================================================
  // HIGH-ORDER GLL INITIALIZATION
  // ==========================================================================
  void initGlobalNodeList()
  {
    int n_local = (order_ + 1) * (order_ + 1) * (order_ + 1);
    global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(
        ex_ * ey_ * ez_, n_local, "global node index");

    int nx = ex_ * order_ + 1;
    int ny = ey_ * order_ + 1;

    for (int k = 0; k < ez_; k++)
    {
      for (int j = 0; j < ey_; j++)
      {
        for (int i = 0; i < ex_; i++)
        {
          int elementNum = i + j * ex_ + k * ex_ * ey_;
          int offset = i * order_ + j * order_ * nx + k * order_ * nx * ny;
          for (int m = 0; m <= order_; m++)
          {
            for (int n = 0; n <= order_; n++)
            {
              for (int l = 0; l <= order_; l++)
              {
                int localIdx =
                    l + n * (order_ + 1) + m * (order_ + 1) * (order_ + 1);
                int globalIdx = offset + l + n * nx + m * nx * ny;
                global_node_index_(elementNum, localIdx) = globalIdx;
              }
            }
          }
        }
      }
    }
  }

  void initNodesCoords()
  {
    int nx = ex_ * order_ + 1;
    int ny = ey_ * order_ + 1;
    int nz = ez_ * order_ + 1;
    int total_nodes = nx * ny * nz;

    nodes_coords_x_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "coords x");
    nodes_coords_y_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "coords y");
    nodes_coords_z_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "coords z");

    float cx[MAX_ORDER + 1], cy[MAX_ORDER + 1], cz[MAX_ORDER + 1];
    auto hx = lx_ / ex_;
    auto hy = ly_ / ey_;
    auto hz = lz_ / ez_;

    for (int n = 0; n < ez_; n++)
    {
      getCoordInOneDirection(hz, n, cz);
      for (int m = 0; m < ey_; m++)
      {
        getCoordInOneDirection(hy, m, cy);
        for (int l = 0; l < ex_; l++)
        {
          getCoordInOneDirection(hx, l, cx);
          for (int k = 0; k <= order_; k++)
          {
            for (int j = 0; j <= order_; j++)
            {
              for (int i = 0; i <= order_; i++)
              {
                int gi = l * order_ + i;
                int gj = m * order_ + j;
                int gk = n * order_ + k;
                int idx = gi + gj * nx + gk * nx * ny;
                nodes_coords_x_(idx) = cx[i];
                nodes_coords_y_(idx) = cy[j];
                nodes_coords_z_(idx) = cz[k];
              }
            }
          }
        }
      }
    }
  }

  void getCoordInOneDirection(const FloatType h, const int n_element,
                              float* coord)
  {
    float xi[MAX_ORDER + 1];
    switch (order_)
    {
      case 1:
        xi[0] = -1.0f;
        xi[1] = 1.0f;
        break;
      case 2:
        xi[0] = -1.0f;
        xi[1] = 0.0f;
        xi[2] = 1.0f;
        break;
      case 3: {
        static constexpr float s5 = 0.4472135955f;
        xi[0] = -1.0f;
        xi[1] = -s5;
        xi[2] = s5;
        xi[3] = 1.0f;
        break;
      }
      case 4: {
        static constexpr float s3_7 = 0.6546536707f;
        xi[0] = -1.0f;
        xi[1] = -s3_7;
        xi[2] = 0.0f;
        xi[3] = s3_7;
        xi[4] = 1.0f;
        break;
      }
      case 5: {
        static constexpr float s21 = 0.2182178902f;
        static constexpr float s7p = 3.5059239327f;
        static constexpr float s7m = 1.3070950148f;
        xi[0] = -1.0f;
        xi[1] = -s21 * s7p;
        xi[2] = -s21 * s7m;
        xi[3] = s21 * s7m;
        xi[4] = s21 * s7p;
        xi[5] = 1.0f;
        break;
      }
    }
    float x0 = static_cast<float>(n_element) * h;
    float x1 = static_cast<float>(n_element + 1) * h;
    float mid = (x1 + x0) * 0.5f;
    float hw = (x1 - x0) * 0.5f;
    for (int j = 0; j <= order_; j++) coord[j] = hw * xi[j] + mid;
  }

  void initModels()
  {
    int n_el = ex_ * ey_ * ez_;
    if (isModelOnNodes_)
    {
      int n_no = (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);
      model_rho_node_ = allocateVector<VECTOR_REAL_VIEW>(n_no, "rho node");
      model_vp_node_ = allocateVector<VECTOR_REAL_VIEW>(n_no, "vp node");
      for (int i = 0; i < n_no; ++i)
      {
        model_rho_node_[i] = 1.0;
        model_vp_node_[i] = 1500.0;
      }
      if (isElastic_)
      {
        model_vs_node_ = allocateVector<VECTOR_REAL_VIEW>(n_no, "vs node");
        model_delta_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_no, "delta node");
        model_epsilon_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_no, "eps node");
        model_gamma_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_no, "gamma node");
        model_theta_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_no, "theta node");
        model_phi_node_ = allocateVector<VECTOR_REAL_VIEW>(n_no, "phi node");
        for (int i = 0; i < n_no; ++i) model_vs_node_[i] = 755.0;
      }
    }
    else
    {
      model_rho_element_ = allocateVector<VECTOR_REAL_VIEW>(n_el, "rho elem");
      model_vp_element_ = allocateVector<VECTOR_REAL_VIEW>(n_el, "vp elem");
      for (int i = 0; i < n_el; ++i)
      {
        model_rho_element_[i] = 1.0;
        model_vp_element_[i] = 1500.0;
      }
      if (isElastic_)
      {
        model_vs_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_el, "vs element");
        model_delta_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_el, "delta element");
        model_epsilon_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_el, "eps element");
        model_gamma_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_el, "gamma element");
        model_theta_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_el, "theta element");
        model_phi_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_el, "phi element");
        for (int i = 0; i < n_el; ++i) model_vs_element_[i] = 755.0;
      }
    }
  }
};
}  // namespace model
#endif
