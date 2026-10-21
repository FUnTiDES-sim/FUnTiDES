#ifndef FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_

/**
 * @file sem_solver_impl.h
 * @brief Umbrella over the out-of-line definition fragments of SEMsolver.
 *
 * The definitions live in thematic fragments, all included here in one
 * fixed order so that every translation unit including this umbrella
 * compiles the same set of definitions:
 *
 * - sem_solver_step.h     — computeFEInit and the one-step drivers
 *                           (computeForces, updateSolution*, reset*, applyRHSTerm,
 *                           computeElementContributions dispatcher);
 * - sem_solver_kernels.h  — the shared capability probes, z-deformed geometry
 *                           helpers and the attenuation dispatcher;
 * - sem_solver_acoustic_impl.h — the acoustic element kernels;
 * - sem_solver_elastic_impl.h — the elastic element kernels and the TTI
 *                           tensor builders (computeCMatrix, computeTtiCompact);
 * - sem_solver_update.h   — updateFieldsForward / updateFieldsBackward;
 * - sem_solver_init.h     — setup and initialisation: global mass matrix,
 *                           damping matrix, array allocation, sponge values,
 *                           diagnostic output;
 * - sem_solver_element_subset.h — the masked / node-list entry points
 *                           (subset assemblies and field updates).
 *
 * Include the fragment you need when adding new code; keep this umbrella
 * for the explicit-instantiation units.
 */

#include "sem_solver_step.h"

#include "sem_solver_kernels.h"

#include "sem_solver_acoustic_impl.h"

#include "sem_solver_elastic_impl.h"

#include "sem_solver_update.h"

#include "sem_solver_init.h"

#include "sem_solver_element_subset.h"

// Named physics solver types — after all fragments, where SEMsolver's
// definition is complete.
#include "acoustic_solver.h"
#include "elastic_solver.h"

#endif  // FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
