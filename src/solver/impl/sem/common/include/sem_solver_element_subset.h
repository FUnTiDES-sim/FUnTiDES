#ifndef FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_ELEMENT_SUBSET_H_
#define FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_ELEMENT_SUBSET_H_

/**
 * @file sem_solver_element_subset.h
 * @brief Entry points operating on a subset of elements or nodes: the masked
 * mass/damping assemblies, the node-list element contributions and the
 * node-list field updates. These serve the p-adaptive and interface-coupling
 * paths, where only a subset of the domain is assembled or advanced.
 */

#include <data_type.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>

#include "Integrals.h"
#include "elastic_flux.h"
#include "mesh_type_traits.h"
#include "sem_solver.h"

namespace solver {
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrixMasked(
    const vectorInt& elem_mask, int active_value) {
  Kokkos::deep_copy(massMatrixGlobal_, 0.0f);
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeGlobalMassMatrix();
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrixMasked(
    const vectorInt& elem_mask, int active_value) {
  for (int f = 0; f < kNumFields; ++f) Kokkos::deep_copy(dampingMatrixGlobal_[f], 0.0f);
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeDampingMatrix();
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsMasked(
    const DataType& data, const vectorInt& elem_mask, int active_value) {
  m_element_mask_ = elem_mask;
  m_mask_active_value_ = active_value;
  m_mask_enabled_ = true;
  computeElementContributions(data);
  m_mask_enabled_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributionsFromList(
    const DataType& data, const vectorInt& elem_list, int n_elems) {
  m_elem_list_ = elem_list;
  m_n_elem_list_ = n_elems;
  m_list_mode_ = true;
  computeElementContributions(data);
  m_list_mode_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListForward(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_node_list_ = node_list;
  m_n_node_list_ = n_nodes;
  m_node_list_mode_ = true;
  updateFieldsForward(dt, data);
  m_node_list_mode_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListForwardAndReset(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_reset_forces_in_update_ = true;
  updateFieldsFromListForward(dt, data, node_list, n_nodes);
  m_reset_forces_in_update_ = false;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::updateFieldsFromListBackward(
    float dt, const DataType& data, const vectorInt& node_list, int n_nodes) {
  m_node_list_ = node_list;
  m_n_node_list_ = n_nodes;
  m_node_list_mode_ = true;
  updateFieldsBackward(dt, data);
  m_node_list_mode_ = false;
}


}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_ELEMENT_SUBSET_H_
