#ifndef FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_
#define FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <memory>

#include "common_macros.h"
#include "rhs_acoustic.h"
#include "rhs_elastic.h"
#include "sem_solver.h"
#include "solver_factory.h"
#include "wavefield_acoustic.h"
#include "wavefield_elastic.h"

namespace py = pybind11;

namespace solver
{
namespace fe
{

void bind_data_struct(py::module_ &m)
{
  py::class_<Solver::DataStruct, std::shared_ptr<Solver::DataStruct>>(
      m, "DataStruct")
      .def("print", &Solver::DataStruct::print);
}

void bind_acoustic_solver_data(py::module_ &m)
{
  py::class_<SEMsolverDataAcoustic, Solver::DataStruct,
             std::shared_ptr<SEMsolverDataAcoustic>>(m, "SEMsolverDataAcoustic")
      .def(py::init<const WavefieldAcoustic &, const RhsAcoustic &>(),
           py::arg("wavefield"), py::arg("rhs"))
      .def("swap_wavefields", &SEMsolverDataAcoustic::swapWavefields)
      .def("print", &SEMsolverDataAcoustic::print);
}

void bind_elastic_solver_data(py::module_ &m)
{
  py::class_<SEMsolverDataElastic, Solver::DataStruct,
             std::shared_ptr<SEMsolverDataElastic>>(m, "SEMsolverDataElastic")
      .def(py::init<const WavefieldElastic &, const RhsElastic &>(),
           py::arg("wavefield"), py::arg("rhs"))
      .def("swap_wavefields", &SEMsolverDataElastic::swapWavefields)
      .def("print", &SEMsolverDataElastic::print);
}

void bind_sem_solver_base(py::module_ &m)
{
  py::class_<Solver, std::shared_ptr<Solver>>(m, "Solver")
      .def("compute_fe_init", &Solver::computeFEInit, py::arg("model"),
           py::arg("sponge_size") = std::array<float, 3>{0.0f, 0.0f, 0.0f},
           py::arg("sponge_surface") = true, py::arg("taper_delta") = 0)
      .def("compute_one_step", &Solver::computeOneStep, py::arg("dt"),
           py::arg("time_sample"), py::arg("data"))
      .def("compute_forces", &Solver::computeForces, py::arg("dt"),
           py::arg("time_sample"), py::arg("data"))
      .def("update_solution", &Solver::updateSolution, py::arg("dt"),
           py::arg("data"))
      .def(
          "get_mass_matrix",
          [](Solver &self)
              -> Kokkos::Experimental::python_view_type_t<VECTOR_REAL_VIEW> {
            return self.getMassMatrix();
          },
          py::return_value_policy::reference_internal)
      .def("output_solution_values", &Solver::outputSolutionValues,
           py::arg("t"), py::arg("e"), py::arg("field_global"),
           py::arg("field_name"))
      .def(
          "set_sls_attenuation",
          [](Solver &self, const std::vector<float> &freqs,
             const std::vector<float> &coeffs) {
            VECTOR_REAL_VIEW vf;
            if (!freqs.empty())
            {
              vf = allocateVector<VECTOR_REAL_VIEW>(freqs.size(), "sls_freqs");
              for (size_t i = 0; i < freqs.size(); ++i) vf[i] = freqs[i];
            }
            VECTOR_REAL_VIEW vc;
            if (!coeffs.empty())
            {
              vc =
                  allocateVector<VECTOR_REAL_VIEW>(coeffs.size(), "sls_coeffs");
              for (size_t i = 0; i < coeffs.size(); ++i) vc[i] = coeffs[i];
            }
            self.setSLSAttenuation(vf, vc);
          },
          py::arg("reference_frequencies"),
          py::arg("anelasticity_coefficients") = std::vector<float>{});
}

void bind_solver_factory(py::module_ &m)
{
  m.def(
      "create_solver",
      [](utils::enums::methodType method, utils::enums::implemType implem,
         utils::enums::meshType mesh,
         utils::enums::modelLocationType modelLocation,
         utils::enums::physicType physic, int order) {
        auto solver = solver_factory::createSolver(
            method, implem, mesh, modelLocation, physic, order);
        return std::shared_ptr<Solver>(std::move(solver));
      },
      py::arg("method_type"), py::arg("implem_type"), py::arg("mesh_type"),
      py::arg("model_location"), py::arg("physic_type"), py::arg("order"));
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_SOLVER_H_