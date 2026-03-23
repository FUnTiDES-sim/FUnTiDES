#ifndef FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_GRADIENTS_H_
#define FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_GRADIENTS_H_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <memory>

#include "common_macros.h"
#include "../../api/include/differentiator.h"
#include "../../api/include/gradient_data.h"
#include "gradients.h"
#include "gradients_acoustic.h"
#include "gradients_elastic.h"

namespace py = pybind11;

namespace gradient
{

/**
 * @brief Bind Gradients base class to Python.
 *
 * Exposes the abstract Gradients interface with virtual methods:
 * - get_num_grads()
 * - get_grads_names()
 * - get_current_grads(index)
 * - print()
 */
void bind_gradients_base(py::module_ &m)
{
  py::class_<Gradients, std::shared_ptr<Gradients>>(m, "Gradients")
      .def("get_num_grads", &Gradients::getNumGrads)
      .def("get_grads_names", &Gradients::getGradsNames)
      .def(
          "get_current_grads",
          [](Gradients &self, int index)
              -> Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW> {
            return self.getCurrentGrads(index);
          },
          py::arg("index"), py::return_value_policy::reference_internal)
      .def("print", &Gradients::print);
}

/**
 * @brief Bind GradientsAcoustic to Python.
 *
 * Exposes acoustic-specific accessors:
 * - get_grad_kappa() -> returns Kokkos View
 * - get_grad_buoyancy() -> returns Kokkos View
 *
 * Views are returned with reference_internal policy to prevent data copies.
 */
void bind_gradients_acoustic(py::module_ &m)
{
  py::class_<GradientsAcoustic, Gradients, std::shared_ptr<GradientsAcoustic>>(
      m, "GradientsAcoustic")
      .def(py::init<const VECTOR_REAL_VIEW &, const VECTOR_REAL_VIEW &>(),
           py::arg("grad_kappa"), py::arg("grad_buoyancy"))
      .def("get_num_grads", &GradientsAcoustic::getNumGrads)
      .def("get_grads_names", &GradientsAcoustic::getGradsNames)
      .def(
          "get_grad_kappa",
          [](GradientsAcoustic &self)
              -> Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW> {
            return self.getCurrentGrads(0);
          },
          py::return_value_policy::reference_internal)
      .def(
          "get_grad_buoyancy",
          [](GradientsAcoustic &self)
              -> Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW> {
            return self.getCurrentGrads(1);
          },
          py::return_value_policy::reference_internal)
      .def("print", &GradientsAcoustic::print);
}

/**
 * @brief Bind GradientsElastic to Python.
 *
 * Exposes elastic-specific accessors:
 * - get_grad_rho() -> returns Kokkos View
 * - get_grad_lambda() -> returns Kokkos View
 * - get_grad_mu() -> returns Kokkos View
 *
 * Views are returned with reference_internal policy to prevent data copies.
 */
void bind_gradients_elastic(py::module_ &m)
{
  py::class_<GradientsElastic, Gradients, std::shared_ptr<GradientsElastic>>(
      m, "GradientsElastic")
      .def(py::init<const VECTOR_REAL_VIEW &, const VECTOR_REAL_VIEW &,
                    const VECTOR_REAL_VIEW &>(),
           py::arg("grad_rho"), py::arg("grad_lambda"), py::arg("grad_mu"))
      .def("get_num_grads", &GradientsElastic::getNumGrads)
      .def("get_grads_names", &GradientsElastic::getGradsNames)
      .def(
          "get_grad_rho",
          [](GradientsElastic &self)
              -> Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW> {
            return self.getCurrentGrads(0);
          },
          py::return_value_policy::reference_internal)
      .def(
          "get_grad_lambda",
          [](GradientsElastic &self)
              -> Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW> {
            return self.getCurrentGrads(1);
          },
          py::return_value_policy::reference_internal)
      .def(
          "get_grad_mu",
          [](GradientsElastic &self)
              -> Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW> {
            return self.getCurrentGrads(2);
          },
          py::return_value_policy::reference_internal)
      .def("print", &GradientsElastic::print);
}

/**
 * @brief Bind GradientData to Python.
 *
 * Exposes the gradient data container:
 * - gradients: access to underlying Gradients object
 * - print()
 */
void bind_gradient_data(py::module_ &m)
{
  py::class_<GradientData, std::shared_ptr<GradientData>>(m, "GradientData")
      .def(py::init<>())
      .def(py::init<const VECTOR_REAL_VIEW &, const VECTOR_REAL_VIEW &>(),
           py::arg("grad_kappa"), py::arg("grad_buoyancy"))
      .def(py::init<const VECTOR_REAL_VIEW &, const VECTOR_REAL_VIEW &,
                    const VECTOR_REAL_VIEW &>(),
           py::arg("grad_rho"), py::arg("grad_lambda"), py::arg("grad_mu"))
      .def_readwrite("gradients", &GradientData::gradients)
      .def("print", &GradientData::print);
}

/**
 * @brief Bind Differentiator to Python.
 *
 * Exposes the abstract gradient computation interface:
 * - compute(forward_field, adjoint_fields, grad_data)
 * - get_order()
 * - is_model_on_nodes()
 * - print()
 *
 * This is the INDEPENDENT gradient computation interface - not part of solver.
 */
void bind_differentiator(py::module_ &m)
{
  py::class_<Differentiator, std::shared_ptr<Differentiator>>(
      m, "Differentiator")
      .def(
          "compute",
          [](Differentiator &self, py::handle forward_py,
             py::list adjoint_py_list, GradientData &grad_data) {
            // Convert numpy arrays to Kokkos views
            auto forward = Kokkos::Experimental::python_view_from_numpy<
                VECTOR_REAL_VIEW>(forward_py);

            // Convert list of adjoint arrays
            std::vector<VECTOR_REAL_VIEW> adjoints;
            for (int i = 0; i < adjoint_py_list.size(); ++i)
            {
              adjoints.push_back(
                  Kokkos::Experimental::python_view_from_numpy<
                      VECTOR_REAL_VIEW>(adjoint_py_list[i]));
            }

            self.compute(forward, adjoints, grad_data);
          },
          py::arg("forward_field"), py::arg("adjoint_fields"),
          py::arg("grad_data"))
      .def("get_order", &Differentiator::getOrder)
      .def("is_model_on_nodes", &Differentiator::isModelOnNodes)
      .def("print", &Differentiator::print);
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_GRADIENTS_H_
