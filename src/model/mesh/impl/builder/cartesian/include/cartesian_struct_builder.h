#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_

#pragma once

#include <builder.h>
#include <model_struct.h>

namespace model
{
template <typename FloatType, typename ScalarType, int Order>
class CartesianStructBuilder : public ModelBuilderBase<FloatType, ScalarType>
{
 public:
  CartesianStructBuilder(ScalarType ex, FloatType lx, ScalarType ey,
                         FloatType ly, ScalarType ez, FloatType lz,
                         bool isModelOnNodes, bool isElastic,
                         FloatType ox = 0.0, FloatType oy = 0.0,
                         FloatType oz = 0.0,
                         FloatType global_lx = -1.0, FloatType global_ly = -1.0,
                         FloatType global_lz = -1.0, FloatType global_ox = 0.0,
                         FloatType global_oy = 0.0, FloatType global_oz = 0.0,
                         bool isAcoustoElastic = false,
                         FloatType acoustoElasticBoundaryZ =
                             static_cast<FloatType>(0))
      : ex_(ex),
        ey_(ey),
        ez_(ez),
        lx_(lx),
        ly_(ly),
        lz_(lz),
        isModelOnNodes_(isModelOnNodes),
        isElastic_(isElastic),
        ox_(ox),
        oy_(oy),
        oz_(oz),
        global_lx_(global_lx < 0 ? lx : global_lx),
        global_ly_(global_ly < 0 ? ly : global_ly),
        global_lz_(global_lz < 0 ? lz : global_lz),
        global_ox_(global_ox),
        global_oy_(global_oy),
        global_oz_(global_oz),
        isAcoustoElastic_(isAcoustoElastic),
        acoustoElasticBoundaryZ_(acoustoElasticBoundaryZ)
  {
  }

  ~CartesianStructBuilder() = default;

  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel()
      const override
  {
    model::ModelStructData<FloatType, ScalarType> data;
    data.ex_ = ex_;
    data.ey_ = ey_;
    data.ez_ = ez_;
    data.dx_ = lx_;
    data.dy_ = ly_;
    data.dz_ = lz_;
    data.ox_ = ox_;
    data.oy_ = oy_;
    data.oz_ = oz_;
    data.isModelOnNodes_ = isModelOnNodes_;
    data.isElastic_ = isElastic_;

    auto temp_model = model::ModelStruct<FloatType, ScalarType, Order>(data);

    const int n_node = temp_model.getNumberOfNodes();
    FloatType tol = temp_model.getMinSpacing() * 1e-4;

    FloatType x_min = global_ox_, x_max = global_ox_ + global_lx_;
    FloatType y_min = global_oy_, y_max = global_oy_ + global_ly_;
    FloatType z_min = global_oz_, z_max = global_oz_ + global_lz_;

    auto boundaries_t =
        allocateVector<VECTOR_REAL_VIEW>(n_node, "boundaries_t");

    for (int n = 0; n < n_node; ++n)
    {
      FloatType x = temp_model.nodeCoord(n, 0);
      FloatType y = temp_model.nodeCoord(n, 1);
      FloatType z = temp_model.nodeCoord(n, 2);

      bool at_xmin = (fabs(x - x_min) < tol);
      bool at_xmax = (fabs(x - x_max) < tol);
      bool at_ymin = (fabs(y - y_min) < tol);
      bool at_ymax = (fabs(y - y_max) < tol);
      bool at_zmin = (fabs(z - z_min) < tol);
      bool at_zmax = (fabs(z - z_max) < tol);

      bool on_boundary =
          at_xmin || at_xmax || at_ymin || at_ymax || at_zmin || at_zmax;

      if (!on_boundary)
        boundaries_t(n) = static_cast<FloatType>(BoundaryFlag::InteriorNode);
      else if (at_zmax)
        boundaries_t(n) = static_cast<FloatType>(BoundaryFlag::Surface);
      else
        boundaries_t(n) = static_cast<FloatType>(BoundaryFlag::Damping);
    }

    data.boundaries_t_ = boundaries_t;

    // Bicouche model for acoustoelastic: fluid layer (z >= boundary) vs solid.
    // vs=0 in the fluid → TagElements classifies it as acoustic.
    if (isAcoustoElastic_)
    {
      if (isModelOnNodes_)
      {
        data.model_vp_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model_vp_node");
        data.model_vs_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model_vs_node");
        data.model_rho_node_ =
            allocateVector<VECTOR_REAL_VIEW>(n_node, "model_rho_node");

        for (int n = 0; n < n_node; ++n)
        {
          bool const is_fluid =
              (temp_model.nodeCoord(n, 2) >= acoustoElasticBoundaryZ_);
          data.model_vp_node_[n] =
              is_fluid ? static_cast<FloatType>(1500) : static_cast<FloatType>(3000);
          data.model_vs_node_[n] =
              is_fluid ? static_cast<FloatType>(0) : static_cast<FloatType>(1500);
          data.model_rho_node_[n] =
              is_fluid ? static_cast<FloatType>(1000) : static_cast<FloatType>(2000);
        }
      }
      else
      {
        int const n_elem = ex_ * ey_ * ez_;
        FloatType const hz = lz_ / ez_;

        data.model_vp_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_elem, "model_vp_elem");
        data.model_vs_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_elem, "model_vs_elem");
        data.model_rho_element_ =
            allocateVector<VECTOR_REAL_VIEW>(n_elem, "model_rho_elem");

        for (int k = 0; k < ez_; ++k)
        {
          FloatType const centroid_z =
              oz_ + (k + static_cast<FloatType>(0.5)) * hz;
          bool const is_fluid = (centroid_z >= acoustoElasticBoundaryZ_);
          for (int j = 0; j < ey_; ++j)
          {
            for (int i = 0; i < ex_; ++i)
            {
              int const e = i + j * ex_ + k * ex_ * ey_;
              data.model_vp_element_[e] = is_fluid
                                              ? static_cast<FloatType>(1500)
                                              : static_cast<FloatType>(3000);
              data.model_vs_element_[e] =
                  is_fluid ? static_cast<FloatType>(0)
                           : static_cast<FloatType>(1500);
              data.model_rho_element_[e] = is_fluid
                                               ? static_cast<FloatType>(1000)
                                               : static_cast<FloatType>(2000);
            }
          }
        }
      }
    }

    auto model =
        std::make_shared<model::ModelStruct<FloatType, ScalarType, Order>>(
            data);

    model->buildFaceConnectivity();

    return model;
  }

 private:
  FloatType ox_, oy_, oz_;  // Local origin coordinate in 3D
  FloatType global_ox_, global_oy_, global_oz_;  // Global origin
  ScalarType ex_, ey_, ez_;  // Number of elements for each axis (local)
  FloatType lx_, ly_, lz_;   // Domain size (local)
  FloatType global_lx_, global_ly_, global_lz_;  // Domain size (global)
  bool isModelOnNodes_;
  bool isElastic_;
  bool isAcoustoElastic_{false};
  FloatType acoustoElasticBoundaryZ_{static_cast<FloatType>(0)};
};
}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_STRUCT_BUILDER_H_
