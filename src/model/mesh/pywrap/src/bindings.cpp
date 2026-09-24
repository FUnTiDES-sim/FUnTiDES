#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>

#include "bindings_builder.h"
#include "bindings_face_connectivity.h"
#include "bindings_model.h"
#include "bindings_partionner.h"

namespace py = pybind11;

/**
 * @brief Python module exposing the mesh and model classes.
 *
 * Every class template is registered for the four combinations of
 * FloatType (float, double) and ScalarType (int, long); the order-templated
 * classes (ModelStruct, CartesianStructBuilder) are also registered for each
 * order 1 to 9. Base classes are registered before their derived classes.
 */
PYBIND11_MODULE(model, m) {
  m.attr("__name__") = "pyfuntides.model";

  bindings::bindFaceConnectivityUnstruct<float, int>(m);
  bindings::bindFaceConnectivityUnstruct<double, int>(m);
  bindings::bindFaceConnectivityUnstruct<float, long>(m);
  bindings::bindFaceConnectivityUnstruct<double, long>(m);

  model::bind_anisotropy_type(m);

  model::bind_boundary_flag(m);

  model::bind_modelapi<float, int>(m);
  model::bind_modelapi<double, int>(m);
  model::bind_modelapi<float, long>(m);
  model::bind_modelapi<double, long>(m);

  for (int order = 1; order <= 9; ++order) {
    model::orderDispatch(order, [&](auto order_tag) {
      constexpr int ord = decltype(order_tag)::value;
      model::bind_modelstruct<float, int, ord>(m);
      model::bind_modelstruct<double, int, ord>(m);
      model::bind_modelstruct<float, long, ord>(m);
      model::bind_modelstruct<double, long, ord>(m);
      return nullptr;
    });
  }

  model::bind_modelstructdata<float, int>(m);
  model::bind_modelstructdata<double, int>(m);
  model::bind_modelstructdata<float, long>(m);
  model::bind_modelstructdata<double, long>(m);

  model::bind_modelunstruct<float, int>(m);
  model::bind_modelunstruct<double, int>(m);
  model::bind_modelunstruct<float, long>(m);
  model::bind_modelunstruct<double, long>(m);

  model::bind_modelunstructdata<float, int>(m);
  model::bind_modelunstructdata<double, int>(m);
  model::bind_modelunstructdata<float, long>(m);
  model::bind_modelunstructdata<double, long>(m);

  model::bind_modelbuilderbase<float, int>(m);
  model::bind_modelbuilderbase<double, int>(m);
  model::bind_modelbuilderbase<float, long>(m);
  model::bind_modelbuilderbase<double, long>(m);

  for (int order = 1; order <= 9; ++order) {
    model::orderDispatch(order, [&](auto order_tag) {
      constexpr int ord = decltype(order_tag)::value;
      model::bind_cartesian_struct_builder<float, int, ord>(m);
      model::bind_cartesian_struct_builder<double, int, ord>(m);
      model::bind_cartesian_struct_builder<float, long, ord>(m);
      model::bind_cartesian_struct_builder<double, long, ord>(m);
      return nullptr;
    });
  }

  model::bind_cartesian_unstruct_params<float, int>(m);
  model::bind_cartesian_unstruct_params<double, int>(m);
  model::bind_cartesian_unstruct_params<float, long>(m);
  model::bind_cartesian_unstruct_params<double, long>(m);

  model::bind_cartesian_unstruct_builder<float, int>(m);
  model::bind_cartesian_unstruct_builder<double, int>(m);
  model::bind_cartesian_unstruct_builder<float, long>(m);
  model::bind_cartesian_unstruct_builder<double, long>(m);

  model::bind_cartesian_partitioner<float, int>(m);
  model::bind_cartesian_partitioner<double, int>(m);
  model::bind_cartesian_partitioner<float, long>(m);
  model::bind_cartesian_partitioner<double, long>(m);
}
