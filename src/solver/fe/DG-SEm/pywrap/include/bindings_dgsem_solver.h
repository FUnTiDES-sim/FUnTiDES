#ifndef FUNTIDES_SOLVER_FE_DG_SEM_PYWRAP_INCLUDE_BINDINGS_DGSEM_SOLVER_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_PYWRAP_INCLUDE_BINDINGS_DGSEM_SOLVER_H_
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "dg-sem_rhs_acoustic.h"
#include "dg-sem_solver_data.h"
#include "rhs.h"

namespace py = pybind11;

namespace solver {
namespace fe {

/**
 * @brief Registers the coupled DG-SEM acoustic source class in a Python module.
 * @param[in,out] m Python module receiving the class.
 * @note The base class Rhs must already be registered in Python.
 */
void bind_dgsem_rhs_acoustic(py::module_ &m) {
  py::class_<DGSEMRhsAcoustic, Rhs, std::shared_ptr<DGSEMRhsAcoustic>>(m, "DGSEMRhsAcoustic")
      .def(py::init<Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>,
                    Kokkos::Experimental::python_view_type_t<vectorInt>,
                    Kokkos::Experimental::python_view_type_t<arrayReal>>(),
           py::arg("dg_acoustic_term"), py::arg("sem_acoustic_term"), py::arg("element"), py::arg("weights"))
      .def("print", &DGSEMRhsAcoustic::print);
}

/**
 * @brief Registers the coupled DG-SEM acoustic solver data class in a Python module.
 * @param[in,out] m Python module receiving the class.
 * @note The base class Solver::DataStruct must already be registered in Python.
 */
void bind_dgsem_acoustic_data(py::module_ &m) {
  py::class_<DGSEMsolverData, Solver::DataStruct, std::shared_ptr<DGSEMsolverData>>(m, "DGSEMsolverData")
      .def(py::init<const DGSEMWavefieldAcoustic &, const DGSEMRhsAcoustic &>(), py::arg("wavefield"), py::arg("rhs"))
      .def("swap_wavefields", &DGSEMsolverData::swapWavefields)
      .def("print", &DGSEMsolverData::print);
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_DG_SEM_PYWRAP_INCLUDE_BINDINGS_DGSEM_SOLVER_H_
