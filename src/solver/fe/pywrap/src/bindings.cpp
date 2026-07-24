#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bindings_device_utils.h"
#include "bindings_rhs.h"
#include "bindings_sem_enums.h"
#include "bindings_sem_solver.h"
#include "bindings_wavefield.h"

#ifdef COMPILE_DG
#include "bindings_dg_solver.h"
#include "bindings_dg_wavefield.h"
#endif
#ifdef COMPILE_DG_SEM
#include "bindings_dgsem_solver.h"
#include "bindings_dgsem_wavefield.h"
#endif
#ifdef COMPILE_DG_PADAPTIVE
#include "bindings_dgpadaptive_solver.h"
#include "bindings_dgpadaptive_wavefield.h"
#endif

namespace py = pybind11;

PYBIND11_MODULE(solver, m) {
  // Create submodule 'solver'
  m.attr("__name__") = "pyfuntides.solver";

  // Bind Enums
  solver::fe::bind_all_sem_enums(m);

  // Bind device-side utilities (ghost exchange, trace gather, diagnostics)
  solver::fe::bind_device_utils(m);

  // Bind RHS
  solver::fe::bind_rhs_base(m);
  solver::fe::bind_rhs_acoustic(m);
  solver::fe::bind_rhs_elastic(m);
  solver::fe::bind_rhs_acoustoelastic(m);

  // Bind Wavefield
  solver::fe::bind_wavefield_base(m);
  solver::fe::bind_wavefield_acoustic(m);
  solver::fe::bind_wavefield_elastic(m);
  solver::fe::bind_wavefield_acoustoelastic(m);

  // Bind Data Structures
  solver::fe::bind_data_struct(m);
  solver::fe::bind_acoustic_solver_data(m);
  solver::fe::bind_elastic_solver_data(m);
  solver::fe::bind_acoustoelastic_solver_data(m);

  // Bind Solver
  solver::fe::bind_sem_solver_base(m);
  solver::fe::bind_solver_factory(m);

#ifdef COMPILE_DG
  // Bind DG
  solver::fe::bind_dg_wavefield_acoustic(m);
  solver::fe::bind_dg_acoustic_data(m);
#endif

#ifdef COMPILE_DG_SEM
  // Bind DG-SEM
  solver::fe::bind_dgsem_wavefield_acoustic(m);
  solver::fe::bind_dgsem_rhs_acoustic(m);
  solver::fe::bind_dgsem_acoustic_data(m);
#endif

#ifdef COMPILE_DG_PADAPTIVE
  // Bind DG p-adaptive
  solver::fe::bind_dgpadaptive_wavefield_acoustic(m);
  solver::fe::bind_dgpadaptive_rhs_acoustic(m);
  solver::fe::bind_dgpadaptive_acoustic_data(m);
#endif
}
