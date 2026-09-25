#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bindings_differentiator.h"
#include "bindings_gradient.h"
#include "bindings_wavefield_view.h"

namespace py = pybind11;

/**
 * @brief Python module entry point, exposed as `pyfuntides.gradient`.
 *
 * Registration order matters: base classes and enums must be bound before the
 * classes that derive from them or use them.
 */
PYBIND11_MODULE(gradient, m) {
  m.attr("__name__") = "pyfuntides.gradient";

  gradient::bind_wavefield_view_base(m);
  gradient::bind_wavefield_view_forward_acoustic(m);
  gradient::bind_wavefield_view_backward_acoustic(m);
  gradient::bind_wavefield_view_forward_elastic(m);
  gradient::bind_wavefield_view_backward_elastic(m);

  gradient::bind_gradient_base(m);
  gradient::bind_gradient_acoustic(m);
  gradient::bind_gradient_elastic(m);

  gradient::bind_data_struct(m);
  gradient::bind_gradient_data_acoustic(m);
  gradient::bind_gradient_data_elastic(m);

  gradient::bind_differentiator_base(m);
  gradient::bind_differentiator_factory(m);
}
