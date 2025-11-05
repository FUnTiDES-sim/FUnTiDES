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

namespace model
{

// template binder for ModelBuilderBase
template <typename FloatType, typename ScalarType>
void bind_modelbuilderbase(py::module_ &m)
{
  using T = model::ModelBuilderBase<FloatType, ScalarType>;
  std::string name =
      model_class_name<FloatType, ScalarType>("ModelBuilderBase");

  py::class_<T, std::shared_ptr<T>>(m, name.c_str())
      .def_static("max_order", []() { return T::MAX_ORDER; })
      .def("get_model", &T::getModel);
}

// template binder for CartesianStructBuilder
template <typename FloatType, typename ScalarType, int Order>
void bind_cartesian_struct_builder(py::module_ &m)
{
  using Base = model::ModelBuilderBase<FloatType, ScalarType>;
  using T = model::CartesianStructBuilder<FloatType, ScalarType, Order>;
  std::string name =
      model_class_name<FloatType, ScalarType, Order>("CartesianStructBuilder");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str())
      .def(py::init<ScalarType, FloatType, ScalarType, FloatType, ScalarType,
                    FloatType, bool, bool>(),
           py::arg("ex"), py::arg("hx"), py::arg("ey"), py::arg("hy"),
           py::arg("ez"), py::arg("hz"), py::arg("is_model_on_nodes"),
           py::arg("is_elastic"));
}

// template binder for CartesianParams
template <typename FloatType, typename ScalarType>
void bind_cartesian_unstruct_params(py::module_ &m)
{
  using Params = model::CartesianParams<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("CartesianParams");

  py::class_<Params, std::shared_ptr<Params>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<int, ScalarType, ScalarType, ScalarType, FloatType,
                    FloatType, FloatType, bool, bool>(),
           py::arg("order"), py::arg("ex"), py::arg("ey"), py::arg("ez"),
           py::arg("lx"), py::arg("ly"), py::arg("lz"),
           py::arg("is_model_on_nodes"), py::arg("is_elastic"))
      .def_readwrite("order", &Params::order)
      .def_readwrite("ex", &Params::ex)
      .def_readwrite("ey", &Params::ey)
      .def_readwrite("ez", &Params::ez)
      .def_readwrite("lx", &Params::lx)
      .def_readwrite("ly", &Params::ly)
      .def_readwrite("lz", &Params::lz)
      .def_readwrite("is_model_on_nodes", &Params::isModelOnNodes)
      .def_readwrite("is_elastic", &Params::isElastic);
}

// template binder for CartesianUnstructBuilder
template <typename FloatType, typename ScalarType>
void bind_cartesian_unstruct_builder(py::module_ &m)
{
  using Base = model::ModelBuilderBase<FloatType, ScalarType>;
  using T = model::CartesianUnstructBuilder<FloatType, ScalarType>;
  std::string name =
      model_class_name<FloatType, ScalarType>("CartesianUnstructBuilder");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str())
      .def(py::init<>())
      .def(py::init<const model::CartesianParams<FloatType, ScalarType> &>());
}

}  // namespace model
