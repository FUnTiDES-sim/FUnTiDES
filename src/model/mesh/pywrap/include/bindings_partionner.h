#ifndef SRC_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_PARTIONNER_H_
#define SRC_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_PARTIONNER_H_

#include <pybind11/pybind11.h>

#include <string>

#include "bindings_utils.h"
#include "cartesian_partitioner.h"

namespace py = pybind11;

namespace model {

/**
 * @brief Registers CartesianXPartitioner<FloatType, ScalarType> in a Python module.
 *
 * The Python class name is built by model_class_name from the base name
 * "CartesianXPartitioner" and the two template arguments. The class exposes a
 * default constructor and partition(global_params, rank, size).
 *
 * @tparam FloatType Floating-point type of the bound partitioner.
 * @tparam ScalarType Integer type of the bound partitioner.
 * @param[in,out] m Python module receiving the class.
 */
template <typename FloatType, typename ScalarType>
void bind_cartesian_partitioner(py::module_ &m) {
  using Partitioner = model::CartesianXPartitioner<FloatType, ScalarType>;

  std::string name = model_class_name<FloatType, ScalarType>("CartesianXPartitioner");

  py::class_<Partitioner>(m, name.c_str())
      .def(py::init<>())
      .def("partition", &Partitioner::partition, py::arg("global_params"), py::arg("rank"), py::arg("size"));
}

}  // namespace model

#endif  // SRC_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_PARTIONNER_H_
