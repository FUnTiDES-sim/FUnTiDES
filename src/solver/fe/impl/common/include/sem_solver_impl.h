#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_

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
 * - sem_solver_kernels.h  — the element kernels (acoustic GEMM/flat/team,
 *                           attenuation, elastic Iso/Tti, z-deformed variants)
 *                           and their detail:: capability probes;
 * - sem_solver_update.h   — updateFieldsForward / updateFieldsBackward;
 * - sem_solver_assembly.h — mass matrix, damping, array allocation, sponge
 *                           values, diagnostics output, computeCMatrix and the
 *                           masked / node-list entry points.
 *
 * The split is purely mechanical: no function body was edited, moved between
 * fragments or reordered (verified by reassembling the fragments into a file
 * byte-identical to the original). Include the fragment you need when adding
 * new code; keep this umbrella for the explicit-instantiation units.
 */

#include "sem_solver_step.h"

#include "sem_solver_kernels.h"

#include "sem_solver_update.h"

#include "sem_solver_assembly.h"

#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_IMPL_H_
