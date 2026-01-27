#ifndef SRC_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_PARTIONNER_H_
#define SRC_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_PARTIONNER_H_

#include <pybind11/pybind11.h>

#include <string>

#include "bindings_utils.h"
#include "cartesian_partitioner.h"

namespace py = pybind11;

namespace model
{

// template binder for CartesianXPartitioner
template <typename FloatType, typename ScalarType>
void bind_cartesian_partitioner(py::module_ &m)
{
  // Define alias for the class we are binding
  using Partitioner = model::CartesianXPartitioner<FloatType, ScalarType>;

  // Generate the Python class name (e.g., CartesianXPartitioner_f32_i32)
  std::string name =
      model_class_name<FloatType, ScalarType>("CartesianXPartitioner");

  // Bind the class
  py::class_<Partitioner>(m, name.c_str())
      .def(py::init<>())
      .def("partition", &Partitioner::partition, py::arg("global_params"),
           py::arg("rank"), py::arg("size"));
}

}  // namespace model

#endif  // SRC_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_PARTIONNER_H_
