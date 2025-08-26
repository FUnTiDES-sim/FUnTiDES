#ifndef SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_
#define SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_

#include "cartesian_params.h"
#include "dataType.hpp"
#include "commonMacros.hpp"
#include <model.hpp>

namespace model_builder {
  template <typename FloatType, typename ScalarType, int Order>
  class CartesianUnstructBuilder : ModelBuilderBase<FloatType, ScalarType>{
   public:
    CartesianUnstructBuilder() { }

    CartesianUnstructBuilder (const CartesianParams & p) {
      ex_ = p.ex;
      ey_ = p.ey;
      ez_ = p.ez;

      lx_ = p.lx;
      ly_ = p.ly;
      lz_ = p.lz;

      initGlobalNodeList();
      initNodesCoords();
    }

    model::ModelUnstruct<FloatType, ScalarType> getModel() {
      model::ModelUnstructData<FloatType, ScalarType> modelData;

      modelData.order_ = order;
      modelData.n_element_ = ex_ * ey_ * ez_;
      modelData.n_node_ = 0;
      modelData.lx_ = lx_;
      modelData.ly_ = ly_;
      modelData.lz_ = lz_;
    }

    ~CartesianUnstructBuilder() = default;

  private:
    ScalarType ex_, ey_, ez_;
    FloatType lx_, ly_, lz_;

    ARRAY_INT_VIEW  global_node_index_;
    VECTOR_REAL_VIEW nodes_coords_x_;
    VECTOR_REAL_VIEW nodes_coords_y_;
    VECTOR_REAL_VIEW nodes_coords_z_;

    void initGlobalNodeList() {
      int nodes_x = Order + 1;
      int nodes_y = Order + 1;
      int nodes_z = Order + 1;
      int total_nodes = nodes_x * nodes_y * nodes_z;
      global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(ex_*ey_*ez_,
                                                          total_nodes,
                                                          "global node index");
      int nx = ex_ * Order + 1;  // Total nodes in x direction
      int ny = ey_ * Order + 1;  // Total nodes in y direction
      int nz = ez_ * Order + 1;  // Total nodes in z direction

      for (int k = 0; k < ez_; k++) {
          for (int j = 0; j < ey_; j++) {
              for (int i = 0; i < ex_; i++) {
                  int elementNum = i + j * ex_ + k * ex_ * ey_;
                  // Corrected offset calculation
                  int offset = i * Order + j * Order * nx + k * Order * nx * ny;

                  for (int m = 0; m < Order + 1; m++) {      // z-direction
                      for (int n = 0; n < Order + 1; n++) {  // y-direction
                          for (int l = 0; l < Order + 1; l++) {  // x-direction
                              int dofLocal = l + n * (Order + 1) + m * (Order + 1) * (Order + 1);
                              int dofGlobal = offset + l + n * nx + m * nx * ny;
                              global_node_index_(elementNum, dofLocal) = dofGlobal;
                          }
                      }
                  }
              }
          }
      }
  }

    void getCoordInOneDirection(const int& h, const int& n_element, float* coord) {
      float xi[Order+1];

      switch (Order) {
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

      for (int j = 0; j < Order + 1; j++) {
        coord[j] = a * xi[j] + b;
      }
    }

    void initNodesCoords() {
      int nodes_x = ex_ * Order + 1;
      int nodes_y = ey_ * Order + 1;
      int nodes_z = ez_ * Order + 1;
      int total_nodes = nodes_x * nodes_y * nodes_z;

      // Init the structure within mesh
      nodes_coords_x_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords x");
      nodes_coords_y_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords y");
      nodes_coords_z_ = allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords z");

      float coord_x[Order+1];
      float coord_y[Order+1];
      float coord_z[Order+1];

      auto hx = lx_ / ex_;
      auto hy = ly_ / ey_;
      auto hz = lz_ / ez_;

      for (int n = 0; n < ez_; n++) {
        getCoordInOneDirection(hz, n, coord_z);
        for (int m = 0; m < ey_; m++) {
          getCoordInOneDirection(hy, m, coord_y);
          for (int l = 0; l < ex_; l++) {
            getCoordInOneDirection(hx, l, coord_x);

            for (int k = 0; k < Order + 1; k++) {
              for (int j = 0; j < Order + 1; j++) {
                for (int i = 0; i < Order + 1; i++) {
                  int global_i = l * Order + i;
                  int global_j = m * Order + j;
                  int global_k = n * Order + k;

                  int global_node_index = global_i
                                          + global_j * nodes_x
                                          + global_k * nodes_x * nodes_y;

                  if (global_i < nodes_x && global_j < nodes_y && global_k < nodes_z) {
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

  };
}
#endif  // SRC_MODEL_CARTESIANMESH_INCLUDE_CARTESIAN_UNSTRUCT_MESH_H_
