#ifndef FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_VIEW_H_
#define FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_VIEW_H_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include "data_type.h"
#include "wavefield_view.h"
#include "wavefield_view_backward_acoustic.h"
#include "wavefield_view_backward_elastic.h"
#include "wavefield_view_forward_acoustic.h"
#include "wavefield_view_forward_elastic.h"

namespace py = pybind11;

namespace gradient {

/**
 * @brief Registers the Python class @c WavefieldView, exposing only @c print.
 * @param[in,out] m Python module receiving the class.
 * @note Must be called before the registration of the derived view classes.
 */
void bind_wavefield_view_base(py::module_& m) {
  py::class_<WavefieldView, std::shared_ptr<WavefieldView>>(m, "WavefieldView").def("print", &WavefieldView::print);
}

/**
 * @brief Registers the Python class @c WavefieldViewForwardAcoustic.
 *
 * The Python constructor takes one array, @c pn, converted to a Kokkos view.
 * @param[in,out] m Python module receiving the class.
 */
void bind_wavefield_view_forward_acoustic(py::module_& m) {
  py::class_<WavefieldViewForwardAcoustic, WavefieldView, std::shared_ptr<WavefieldViewForwardAcoustic>>(
      m, "WavefieldViewForwardAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<vectorReal>>(), py::arg("pn"))
      .def("print", &WavefieldViewForwardAcoustic::print);
}

/**
 * @brief Registers the Python class @c WavefieldViewBackwardAcoustic.
 *
 * The Python constructor takes three arrays: @c qn, @c qn_prev and @c qn_prev_prev.
 * @param[in,out] m Python module receiving the class.
 * @todo VERIFY: do qn_prev and qn_prev_prev denote the adjoint field at the two previous time steps?
 */
void bind_wavefield_view_backward_acoustic(py::module_& m) {
  py::class_<WavefieldViewBackwardAcoustic, WavefieldView, std::shared_ptr<WavefieldViewBackwardAcoustic>>(
      m, "WavefieldViewBackwardAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>>(),
           py::arg("qn"), py::arg("qn_prev"), py::arg("qn_prev_prev"))
      .def("print", &WavefieldViewBackwardAcoustic::print);
}

/**
 * @brief Registers the Python class @c WavefieldViewForwardElastic.
 *
 * The Python constructor takes three arrays: @c ux_n, @c uy_n and @c uz_n.
 * @param[in,out] m Python module receiving the class.
 */
void bind_wavefield_view_forward_elastic(py::module_& m) {
  py::class_<WavefieldViewForwardElastic, WavefieldView, std::shared_ptr<WavefieldViewForwardElastic>>(
      m, "WavefieldViewForwardElastic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>>(),
           py::arg("ux_n"), py::arg("uy_n"), py::arg("uz_n"))
      .def("print", &WavefieldViewForwardElastic::print);
}

/**
 * @brief Registers the Python class @c WavefieldViewBackwardElastic.
 *
 * The Python constructor takes six arrays: the displacement components @c ux_n, @c uy_n,
 * @c uz_n, then the second time derivative components @c ux_dt2, @c uy_dt2, @c uz_dt2.
 * @param[in,out] m Python module receiving the class.
 */
void bind_wavefield_view_backward_elastic(py::module_& m) {
  py::class_<WavefieldViewBackwardElastic, WavefieldView, std::shared_ptr<WavefieldViewBackwardElastic>>(
      m, "WavefieldViewBackwardElastic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>,
                    Kokkos::Experimental::python_view_type_t<vectorReal>>(),
           py::arg("ux_n"), py::arg("uy_n"), py::arg("uz_n"), py::arg("ux_dt2"), py::arg("uy_dt2"), py::arg("uz_dt2"))
      .def("print", &WavefieldViewBackwardElastic::print);
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_WAVEFIELD_VIEW_H_
