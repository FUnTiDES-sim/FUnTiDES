#ifndef FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_COMMON_H_
#define FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_COMMON_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

#include "common_macros.h"
#include "solver.h"
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

/**
 * @brief Binds the Solver::DataStruct base, which every method's data class derives from.
 *
 * Must be registered before any method's data class.
 */
inline void bind_data_struct(py::module_& m) {
  py::class_<Solver::DataStruct, std::shared_ptr<Solver::DataStruct>>(m, "DataStruct")
      .def("print", &Solver::DataStruct::print);
}

/**
 * @brief Binds the Solver interface (api/include/solver.h) returned by create_solver.
 */
inline void bind_solver_base(py::module_& m) {
  py::class_<Solver, std::shared_ptr<Solver>>(m, "Solver")
      .def("compute_fe_init", &Solver::computeFEInit, py::arg("model"),
           py::arg("sponge_size") = std::array<float, 3>{0.0f, 0.0f, 0.0f}, py::arg("sponge_surface") = true,
           py::arg("taper_delta") = 0)
      .def("compute_one_step", &Solver::computeOneStep, py::arg("dt"), py::arg("time_sample"), py::arg("data"))
      .def("compute_forces", &Solver::computeForces, py::arg("dt"), py::arg("time_sample"), py::arg("data"))
      .def("update_solution_forward", &Solver::updateSolutionForward, py::arg("dt"), py::arg("data"))
      .def("update_solution_backward", &Solver::updateSolutionBackward, py::arg("dt"), py::arg("data"))
      .def(
          "get_mass_matrix",
          [](Solver& self, bool elastic) -> Kokkos::Experimental::python_view_type_t<vectorReal> {
            return elastic ? self.getMassMatrixElastic() : self.getMassMatrixAcoustic();
          },
          py::arg("elastic") = false, py::return_value_policy::reference_internal)
      .def(
          "get_damping_matrix",
          [](Solver& self, int component) -> Kokkos::Experimental::python_view_type_t<vectorReal> {
            return self.getDampingMatrix(component);
          },
          py::arg("component"), py::return_value_policy::reference_internal)
      .def(
          "get_force_vector",
          [](Solver& self, int component) -> Kokkos::Experimental::python_view_type_t<vectorReal> {
            return self.getForceVector(component);
          },
          py::arg("component"), py::return_value_policy::reference_internal)
      .def(
          "get_interface_coupling_coeff",
          [](Solver& self, int component) -> Kokkos::Experimental::python_view_type_t<vectorReal> {
            return self.getInterfaceCouplingCoeff(component);
          },
          py::arg("component"), py::return_value_policy::reference_internal)
      .def("output_solution_values",
           static_cast<void (Solver::*)(const int&, int&, const vectorReal&, const char*)>(
               &Solver::outputSolutionValues),
           py::arg("t"), py::arg("e"), py::arg("field_global"), py::arg("field_name"))
      .def(
          "set_sls_attenuation",
          [](Solver& self, const std::vector<float>& freqs, const std::vector<float>& coeffs) {
            vectorReal vf;
            if (!freqs.empty()) {
              vf = allocateVector<vectorReal>(freqs.size(), "sls_freqs");
              for (size_t i = 0; i < freqs.size(); ++i) vf[i] = freqs[i];
            }
            vectorReal vc;
            if (!coeffs.empty()) {
              vc = allocateVector<vectorReal>(coeffs.size(), "sls_coeffs");
              for (size_t i = 0; i < coeffs.size(); ++i) vc[i] = coeffs[i];
            }
            self.setSLSAttenuation(vf, vc);
          },
          py::arg("reference_frequencies"), py::arg("anelasticity_coefficients") = std::vector<float>{})
      .def("set_anisotropy_type", &Solver::setAnisotropyType, py::arg("anisotropy_type"))
      .def("set_interface_property_convention", &Solver::setInterfacePropertyConvention, py::arg("convention"))
      .def(
          "set_solid_on_interface_nodes",
          [](Solver& self, const py::array_t<float, py::array::c_style | py::array::forcecast>& vp,
             const py::array_t<float, py::array::c_style | py::array::forcecast>& vs,
             const py::array_t<float, py::array::c_style | py::array::forcecast>& rho) {
            if (vp.size() != vs.size() || vp.size() != rho.size())
              throw std::invalid_argument("set_solid_on_interface_nodes: vp, vs and rho differ in size");
            // One bulk copy per array: converting element by element through std::vector
            // costs about a minute per call on a production-size mesh.
            auto to_view = [](const py::array_t<float, py::array::c_style | py::array::forcecast>& a,
                              const char* name) {
              vectorReal out = allocateVector<vectorReal>(a.size(), name);
              Kokkos::View<const float*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> in(a.data(),
                                                                                                        a.size());
              Kokkos::deep_copy(out, in);
              return out;
            };
            self.setSolidOnInterfaceNodes(to_view(vp, "solid_vp_node"), to_view(vs, "solid_vs_node"),
                                          to_view(rho, "solid_rho_node"));
          },
          py::arg("vp"), py::arg("vs"), py::arg("rho"))
      .def(
          "set_element_tags",
          [](Solver& self, const py::array_t<int, py::array::c_style | py::array::forcecast>& tags) {
            vectorInt vt;
            if (tags.size() > 0) {
              vt = allocateVector<vectorInt>(tags.size(), "external_element_tags");
              Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> in(tags.data(),
                                                                                                      tags.size());
              Kokkos::deep_copy(vt, in);
            }
            self.setElementTags(vt);
          },
          py::arg("tags"));
}

}  // namespace solver
#endif  // FUNTIDES_SOLVER_PYWRAP_INCLUDE_BINDINGS_COMMON_H_
