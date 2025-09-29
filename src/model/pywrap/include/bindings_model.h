#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bindings_utils.h"
#include "model.h"
#include "model_struct.h"
#include "model_unstruct.h"

#include <string>

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
      .def("get_model_vp_on_nodes", &T::getModelVpOnNodes)
      .def("get_model_vp_on_element", &T::getModelVpOnElement)
      .def("get_model_rho_on_nodes", &T::getModelRhoOnNodes)
      .def("get_model_rho_on_element", &T::getModelRhoOnElement)
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

  std::string name = model_class_name<FloatType, ScalarType, Order>("ModelStruct");

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

//template binder for ModelUnstructData
template <typename FloatType, typename ScalarType>
void bind_modelunstructdata(py::module_ &m)
{
  using Data = model::ModelUnstructData<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("ModelUnstructData");

  py::class_<Data>(m, name.c_str())
      .def(py::init<>())
      .def_readwrite("order", &Data::order_)
      .def_readwrite("n_element", &Data::n_element_)
      .def_readwrite("n_node", &Data::n_node_)
      .def_readwrite("lx", &Data::lx_)
      .def_readwrite("ly", &Data::ly_)
      .def_readwrite("lz", &Data::lz_)
      .def_readwrite("global_node_index", &Data::global_node_index_)
      .def_readwrite("nodes_coords_x", &Data::nodes_coords_x_)
      .def_readwrite("nodes_coords_y", &Data::nodes_coords_y_)
      .def_readwrite("nodes_coords_z", &Data::nodes_coords_z_)
      .def_readwrite("model_vp_node", &Data::model_vp_node_)
      .def_readwrite("model_vp_element", &Data::model_vp_element_)
      .def_readwrite("model_rho_node", &Data::model_rho_node_)
      .def_readwrite("model_rho_element", &Data::model_rho_element_)
      .def_readwrite("boundaries_t", &Data::boundaries_t_);
}

}