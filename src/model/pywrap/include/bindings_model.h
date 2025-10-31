#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <string>

#include "bindings_utils.h"
#include "common_macros.h"
#include "model.h"
#include "model_struct.h"
#include "model_unstruct.h"

namespace py = pybind11;

namespace model
{

// template binder for ModelAPI
template <typename FloatType, typename ScalarType>
void bind_modelapi(py::module_ &m)
{
  using T = model::ModelApi<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("ModelApi");

  py::class_<T, std::shared_ptr<T>>(m, name.c_str())
      .def("node_coord", &T::nodeCoord)
      .def("global_node_index", &T::globalNodeIndex)
      .def("get_model_vp_on_node", &T::getModelVpOnNodes)
      .def("get_model_vp_on_element", &T::getModelVpOnElement)
      .def("get_model_rho_on_node", &T::getModelRhoOnNodes)
      .def("get_model_rho_on_element", &T::getModelRhoOnElement)
      .def("get_model_vs_on_node", &T::getModelVsOnNodes)
      .def("get_model_vs_on_element", &T::getModelVsOnElement)
      .def("get_model_delta_on_node", &T::getModelDeltaOnNodes)
      .def("get_model_delta_on_element", &T::getModelDeltaOnElement)
      .def("get_model_epsilon_on_node", &T::getModelEpsilonOnNodes)
      .def("get_model_epsilon_on_element", &T::getModelEpsilonOnElement)
      .def("get_model_gamma_on_node", &T::getModelGammaOnNodes)
      .def("get_model_gamma_on_element", &T::getModelGammaOnElement)
      .def("get_model_phi_on_node", &T::getModelPhiOnNodes)
      .def("get_model_phi_on_element", &T::getModelPhiOnElement)
      .def("get_model_theta_on_node", &T::getModelThetaOnNodes)
      .def("get_model_theta_on_element", &T::getModelThetaOnElement)
      .def("get_number_of_elements", &T::getNumberOfElements)
      .def("get_number_of_nodes", &T::getNumberOfNodes)
      .def("get_number_of_points_per_element", &T::getNumberOfPointsPerElement)
      .def("get_order", &T::getOrder)
      .def("boundary_type", &T::boundaryType)
      .def("face_normal", &T::faceNormal)
      .def("domain_size", &T::domainSize)
      .def("get_min_spacing", &T::getMinSpacing)
      .def("get_max_speed", &T::getMaxSpeed);
}

// templated binder for one ModelStruct instantiation
template <typename FloatType, typename ScalarType, int Order>
void bind_modelstruct(py::module_ &m)
{
  using Base = model::ModelApi<FloatType, ScalarType>;
  using T = model::ModelStruct<FloatType, ScalarType, Order>;
  using Data = model::ModelStructData<FloatType, ScalarType>;

  std::string name =
      model_class_name<FloatType, ScalarType, Order>("ModelStruct");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str())
      .def(py::init<const Data &>());
}

// templated binder for ModelStructData
template <typename FloatType, typename ScalarType>
void bind_modelstructdata(py::module_ &m)
{
  using Data = model::ModelStructData<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("ModelStructData");

  py::class_<Data>(m, name.c_str())
      .def(py::init<>())
      .def_readwrite("ex", &Data::ex_)
      .def_readwrite("ey", &Data::ey_)
      .def_readwrite("ez", &Data::ez_)
      .def_readwrite("dx", &Data::dx_)
      .def_readwrite("dy", &Data::dy_)
      .def_readwrite("dz", &Data::dz_);
}

// templated binder for ModelUnstruct
template <typename FloatType, typename ScalarType>
void bind_modelunstruct(py::module_ &m)
{
  using Base = model::ModelApi<FloatType, ScalarType>;
  using T = model::ModelUnstruct<FloatType, ScalarType>;
  using Data = model::ModelUnstructData<FloatType, ScalarType>;

  std::string name = model_class_name<FloatType, ScalarType>("ModelUnstruct");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str())
      .def(py::init<const Data &>());
}

// template binder for ModelUnstructData
template <typename FloatType, typename ScalarType>
void bind_modelunstructdata(py::module_ &m)
{
  using Data = model::ModelUnstructData<FloatType, ScalarType>;
  std::string name =
      model_class_name<FloatType, ScalarType>("ModelUnstructData");

  py::class_<Data>(m, name.c_str())
      .def(py::init<
               ScalarType, ScalarType, ScalarType, FloatType, FloatType,
               FloatType, bool, bool,
               Kokkos::Experimental::python_view_type_t<ARRAY_INT_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>,
               Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW>>(),
           py::arg("order"), py::arg("n_element"), py::arg("n_node"),
           py::arg("lx"), py::arg("ly"), py::arg("lz"),
           py::arg("is_model_on_nodes"), py::arg("is_elastic"),
           py::arg("global_node_index"), py::arg("nodes_coords_x"),
           py::arg("nodes_coords_y"), py::arg("nodes_coords_z"),
           py::arg("model_vp_node"), py::arg("model_vp_element"),
           py::arg("model_rho_node"), py::arg("model_rho_element"),
           py::arg("model_vs_node"), py::arg("model_vs_element"),
           py::arg("model_delta_node"), py::arg("model_delta_element"),
           py::arg("model_epsilon_node"), py::arg("model_epsilon_element"),
           py::arg("model_gamma_node"), py::arg("model_gamma_element"),
           py::arg("model_phi_node"), py::arg("model_phi_element"),
           py::arg("model_theta_node"), py::arg("model_theta_element"),
           py::arg("boundaries_t"));
}

}  // namespace model