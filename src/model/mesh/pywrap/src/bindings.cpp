#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>

#include "bindings_builder.h"
#include "bindings_model.h"

namespace py = pybind11;

PYBIND11_MODULE(model, m)
{
  // Create submodule 'model'
  m.attr("__name__") = "pyfuntides.model";

  // Bind ModelApi
  model::bind_modelapi<float, int>(m);
  model::bind_modelapi<double, int>(m);
  model::bind_modelapi<float, long>(m);
  model::bind_modelapi<double, long>(m);

  // Bind ModelStruct
  for (int order = 1; order <= 3; ++order)
  {
    model::orderDispatch(order, [&](auto order_tag) {
      constexpr int ord = decltype(order_tag)::value;
      model::bind_modelstruct<float, int, ord>(m);
      model::bind_modelstruct<double, int, ord>(m);
      model::bind_modelstruct<float, long, ord>(m);
      model::bind_modelstruct<double, long, ord>(m);
      return nullptr;
    });
  }

  // Bind ModelStructData
  model::bind_modelstructdata<float, int>(m);
  model::bind_modelstructdata<double, int>(m);
  model::bind_modelstructdata<float, long>(m);
  model::bind_modelstructdata<double, long>(m);

  // Bind ModelUnstruct
  model::bind_modelunstruct<float, int>(m);
  model::bind_modelunstruct<double, int>(m);
  model::bind_modelunstruct<float, long>(m);
  model::bind_modelunstruct<double, long>(m);

  // Bind ModelUnstructData
  model::bind_modelunstructdata<float, int>(m);
  model::bind_modelunstructdata<double, int>(m);
  model::bind_modelunstructdata<float, long>(m);
  model::bind_modelunstructdata<double, long>(m);

  // Bind ModelBuilderBase
  model::bind_modelbuilderbase<float, int>(m);
  model::bind_modelbuilderbase<double, int>(m);
  model::bind_modelbuilderbase<float, long>(m);
  model::bind_modelbuilderbase<double, long>(m);

  // Bind CartesianStructBuilder
  for (int order = 1; order <= 3; ++order)
  {
    model::orderDispatch(order, [&](auto order_tag) {
      constexpr int ord = decltype(order_tag)::value;
      model::bind_cartesian_struct_builder<float, int, ord>(m);
      model::bind_cartesian_struct_builder<double, int, ord>(m);
      model::bind_cartesian_struct_builder<float, long, ord>(m);
      model::bind_cartesian_struct_builder<double, long, ord>(m);
      return nullptr;
    });
  }

  // Bind CartesianParams
  model::bind_cartesian_unstruct_params<float, int>(m);
  model::bind_cartesian_unstruct_params<double, int>(m);
  model::bind_cartesian_unstruct_params<float, long>(m);
  model::bind_cartesian_unstruct_params<double, long>(m);

  // Bind CartesianUnstructBuilder
  model::bind_cartesian_unstruct_builder<float, int>(m);
  model::bind_cartesian_unstruct_builder<double, int>(m);
  model::bind_cartesian_unstruct_builder<float, long>(m);
  model::bind_cartesian_unstruct_builder<double, long>(m);
}
