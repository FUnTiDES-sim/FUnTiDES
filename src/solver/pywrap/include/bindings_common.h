#ifndef FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_COMMON_H_
#define FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_COMMON_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <memory>

#include "common_macros.h"
#include "solver_factory.h"
#include "wavefield.h"

namespace py = pybind11;

namespace solver {

/**
 * @brief Binds the Wavefield base interface (api/include/wavefield.h).
 *
 * Method-agnostic: every concrete wavefield of every method derives from it,
 * so this binding belongs to the module's common part, not to any method.
 */
inline void bind_wavefield_base(py::module_& m) {
  // Bind Wavefield (base class)
  py::class_<Wavefield, std::shared_ptr<Wavefield>>(m, "Wavefield")
      .def("swap", &Wavefield::swap)
      .def(
          "get_current_field",
          [](const Wavefield& self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getCurrentField(i));
          },
          py::arg("i"))
      .def(
          "get_previous_field",
          [](const Wavefield& self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPreviousField(i));
          },
          py::arg("i"))
      .def(
          "get_prevprev_field",
          [](const Wavefield& self, int i) {
            return Kokkos::Experimental::python_view_type_t<vectorReal>(self.getPrevPrevField(i));
          },
          py::arg("i"))
      .def("has_prevprev", &Wavefield::hasPrevPrev)
      .def("print", &Wavefield::print);
}

/**
 * @brief Binds solver_factory::createSolver (impl/common).
 *
 * Method-agnostic: the factory dispatches on the method at runtime and every
 * method registers itself through it, so this binding too belongs to the
 * module's common part.
 */
inline void bind_solver_factory(py::module_& m) {
  m.def(
      "create_solver",
      [](utils::enums::methodType method, utils::enums::implemType implem, utils::enums::meshType mesh,
         utils::enums::modelLocationType modelLocation, utils::enums::physicType physic, int order, int order_min) {
        auto solver = solver_factory::createSolver(method, implem, mesh, modelLocation, physic, order, order_min);
        return std::shared_ptr<Solver>(std::move(solver));
      },
      py::arg("method_type"), py::arg("implem_type"), py::arg("mesh_type"), py::arg("model_location"),
      py::arg("physic_type"), py::arg("order"), py::arg("order_min") = 0);
}

}  // namespace solver
#endif  // FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_COMMON_H_
