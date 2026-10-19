#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
  solver::bind_all_sem_enums(m);

  // Bind RHS
  solver::bind_rhs_base(m);
  solver::bind_rhs_acoustic(m);
  solver::bind_rhs_elastic(m);
  solver::bind_rhs_acoustoelastic(m);

  // Bind Wavefield
  solver::bind_wavefield_base(m);
  solver::bind_wavefield_acoustic(m);
  solver::bind_wavefield_elastic(m);
  solver::bind_wavefield_acoustoelastic(m);

  // Bind Data Structures
  solver::bind_data_struct(m);
  solver::bind_acoustic_solver_data(m);
  solver::bind_elastic_solver_data(m);
  solver::bind_acoustoelastic_solver_data(m);

  // Bind Solver
  solver::bind_sem_solver_base(m);
  solver::bind_solver_factory(m);

#ifdef COMPILE_DG
  // Bind DG
  solver::bind_dg_wavefield_acoustic(m);
  solver::bind_dg_acoustic_data(m);
#endif

#ifdef COMPILE_DG_SEM
  // Bind DG-SEM
  solver::bind_dgsem_wavefield_acoustic(m);
  solver::bind_dgsem_rhs_acoustic(m);
  solver::bind_dgsem_acoustic_data(m);
#endif

#ifdef COMPILE_DG_PADAPTIVE
  // Bind DG p-adaptive
  solver::bind_dgpadaptive_wavefield_acoustic(m);
  solver::bind_dgpadaptive_rhs_acoustic(m);
  solver::bind_dgpadaptive_acoustic_data(m);
#endif
}
