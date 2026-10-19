#ifndef FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_

/**
 * @file sem_solver_impl.h
 * @brief Umbrella over the out-of-line definitions of SEMsolver.
 *
 * The former single 2700-line implementation header is split into thematic
 * fragments, all included here in the original definition order so that every
 * translation unit that includes this umbrella compiles the exact same set of
 * definitions as before the split:
 *
 * - sem_solver_step.h     — computeFEInit and the one-step drivers
 *                           (computeForces, updateSolution*, reset*, applyRHSTerm,
 *                           computeElementContributions dispatcher);
 * - sem_solver_kernels.h  — the shared capability probes, z-deformed geometry
 *                           helpers, the attenuation dispatcher, and the TTI
 *                           tensor builders (computeCMatrix, computeTtiCompact);
 * - sem_solver_acoustic_kernels.h — the acoustic element kernels
 *                           (physics-policy migration, step 1: moved verbatim
 *                           from sem_solver_kernels.h);
 * - sem_solver_elastic_kernels.h — the elastic element kernels
 *                           (physics-policy migration, step 2: moved verbatim
 *                           from sem_solver_kernels.h);
 * - sem_solver_update.h   — updateFieldsForward / updateFieldsBackward;
 * - sem_solver_init.h    — setup and initialisation: global mass matrix,
 *                           damping matrix, array allocation, sponge values,
 *                           diagnostic output;
 * - sem_solver_element_subset.h — the masked / node-list entry points
 *                           (subset assemblies and field updates).
 *
 * The split is purely mechanical: no function body was edited, moved between
 * fragments or reordered (verified by reassembling the fragments into a file
 * byte-identical to the original). Include the fragment you need when adding
 * new code; keep this umbrella for the explicit-instantiation units.
 */

#include "sem_solver_step.h"

#include "sem_solver_kernels.h"

#include "sem_solver_acoustic_kernels.h"

#include "sem_solver_elastic_kernels.h"

#include "sem_solver_update.h"

#include "sem_solver_init.h"

#include "sem_solver_element_subset.h"

// Named physics solver types (physics-policy migration, step 3) — after all
// fragments, where SEMsolver's definition is complete.
#include "acoustic_solver.h"
#include "elastic_solver.h"

#endif  // FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
