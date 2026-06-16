#ifndef FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_BUILDER_H_
#define FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_BUILDER_H_

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "bindings_utils.h"
#include "builder.h"
#include "cartesian_struct_builder.h"
#include "cartesian_unstruct_builder.h"

namespace py = pybind11;

namespace model {

// template binder for ModelBuilderBase
template <typename FloatType, typename ScalarType>
void bind_modelbuilderbase(py::module_ &m) {
  using T = model::ModelBuilderBase<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("ModelBuilderBase");

  py::class_<T, std::shared_ptr<T>>(m, name.c_str())
      .def_static("max_order", []() { return T::MAX_ORDER; })
      .def("get_model", &T::getModel, py::arg("free_surface_on_top") = true);
}

// template binder for CartesianStructBuilder
template <typename FloatType, typename ScalarType, int Order>
void bind_cartesian_struct_builder(py::module_ &m) {
  using Base = model::ModelBuilderBase<FloatType, ScalarType>;
  using T = model::CartesianStructBuilder<FloatType, ScalarType, Order>;
  std::string name = model_class_name<FloatType, ScalarType, Order>("CartesianStructBuilder");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str())
      // Standard constructor (no acousto-elastic).
      .def(py::init<ScalarType, FloatType, ScalarType, FloatType, ScalarType, FloatType, bool, bool, FloatType,
                    FloatType, FloatType>(),
           py::arg("ex"), py::arg("hx"), py::arg("ey"), py::arg("hy"), py::arg("ez"), py::arg("hz"),
           py::arg("is_model_on_nodes"), py::arg("is_elastic"), py::arg("ox") = static_cast<FloatType>(0),
           py::arg("oy") = static_cast<FloatType>(0), py::arg("oz") = static_cast<FloatType>(0))
      // Acousto-elastic constructor: exposes the two AE parameters while
      // leaving global domain metrics at their defaults (local = global),
      // which is correct for single-rank / non-MPI usage.
      // For multi-rank, a third overload exposing global_lx/ly/lz would be
      // needed.
      .def(py::init([](ScalarType ex, FloatType hx, ScalarType ey, FloatType hy, ScalarType ez, FloatType hz,
                       bool is_model_on_nodes, bool is_elastic, FloatType ox, FloatType oy, FloatType oz,
                       bool is_acousto_elastic, FloatType acousto_elastic_boundary_z, const std::string &model_file) {
             return T(ex, hx, ey, hy, ez, hz, is_model_on_nodes, is_elastic, ox, oy, oz, static_cast<FloatType>(-1),
                      static_cast<FloatType>(-1), static_cast<FloatType>(-1), static_cast<FloatType>(0),
                      static_cast<FloatType>(0), static_cast<FloatType>(0), is_acousto_elastic,
                      acousto_elastic_boundary_z, static_cast<FloatType>(0), model_file);
           }),
           py::arg("ex"), py::arg("hx"), py::arg("ey"), py::arg("hy"), py::arg("ez"), py::arg("hz"),
           py::arg("is_model_on_nodes"), py::arg("is_elastic"), py::arg("ox") = static_cast<FloatType>(0),
           py::arg("oy") = static_cast<FloatType>(0), py::arg("oz") = static_cast<FloatType>(0),
           py::arg("is_acousto_elastic") = false, py::arg("acousto_elastic_boundary_z") = static_cast<FloatType>(0),
           py::arg("model_file") = std::string{});
}

// template binder for CartesianParams
template <typename FloatType, typename ScalarType>
void bind_cartesian_unstruct_params(py::module_ &m) {
  using Params = model::CartesianParams<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("CartesianParams");

  py::class_<Params, std::shared_ptr<Params>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int, ScalarType, ScalarType, ScalarType, FloatType, FloatType, FloatType, bool, bool>(),
           py::arg("order"), py::arg("ex"), py::arg("ey"), py::arg("ez"), py::arg("lx"), py::arg("ly"), py::arg("lz"),
           py::arg("is_model_on_nodes"), py::arg("is_elastic"))
      .def_readwrite("order", &Params::order)
      .def_readwrite("ex", &Params::ex)
      .def_readwrite("ey", &Params::ey)
      .def_readwrite("ez", &Params::ez)
      .def_readwrite("lx", &Params::lx)
      .def_readwrite("ly", &Params::ly)
      .def_readwrite("lz", &Params::lz)
      .def_readwrite("origin_x", &Params::origin_x)
      .def_readwrite("origin_y", &Params::origin_y)
      .def_readwrite("origin_z", &Params::origin_z)
      .def_readwrite("is_model_on_nodes", &Params::isModelOnNodes)
      .def_readwrite("is_elastic", &Params::isElastic)
      .def_readwrite("is_acousto_elastic", &Params::isAcoustoElastic)
      .def_readwrite("acousto_elastic_boundary_z", &Params::acoustoElasticBoundaryZ);
}

// template binder for CartesianUnstructBuilder
template <typename FloatType, typename ScalarType>
void bind_cartesian_unstruct_builder(py::module_ &m) {
  using Base = model::ModelBuilderBase<FloatType, ScalarType>;
  using T = model::CartesianUnstructBuilder<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("CartesianUnstructBuilder");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<const model::CartesianParams<FloatType, ScalarType> &>());
}

}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_BUILDER_H_
