#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_

#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Integrals.h"
#include "data_type.h"
#include "dg-sem_solver.h"
#include "dg_penalty.h"

namespace solver {
namespace fe {

//============================================================================
// computeFEInit
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& sponge_size, const bool surface_sponge,
    const float taper_delta) {
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh_ = *typed;
  } else {
    throw std::runtime_error("DGSEMsolver: incompatible mesh type in computeFEInit");
  }

  m_face_connectivity_.build(m_mesh_);

  // Initialise sub-solvers (mass/damping matrices are overridden below).
  m_SEm_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  // Before its init: BuildDGInteriorFaceList() below emits ids in this numbering.
  m_DG_solver_.setFaceConnectivity(m_face_connectivity_);
  m_DG_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);

  m_penalty_factor_ = m_DG_solver_.getPenaltyFactor();

  allocateFEarrays();

  TagElements();
  std::cout << "DGSEMsolver: " << num_SEm_elements_ << " SEm elements, " << num_DG_elements_ << " DG elements."
            << std::endl;

  // Re-assemble the SEM sub-solver's global mass/damping matrices restricted to the SEM
  // subdomain. The full-mesh assembly done inside m_SEm_solver_.computeFEInit() above also
  // accumulated contributions from DG-tagged elements, inflating the mass of the shared
  // interface nodes (~2x) — the SEM weak form only owns the SEM elements (stiffness runs on
  // SEm_elem_list_) — which acted as a heavy strip along the DG-SEM interface and produced a
  // spurious partial reflection of waves crossing it.
  m_SEm_solver_.computeGlobalMassMatrixMasked(m_element_type_, kElementTypeSEM);
  m_SEm_solver_.computeDampingMatrixMasked(m_element_type_, kElementTypeSEM);

  TagNodes();
  std::cout << "DGSEMsolver: " << num_interface_faces_ << " interface faces." << std::endl;
}

//============================================================================
// allocateFEarrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::allocateFEarrays() {
  int const nElem = m_mesh_.getNumberOfElements();
  m_element_type_ = allocateVector<vectorInt>(nElem, "DGSEMElementType");
}

//============================================================================
// TagElements
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::TagElements() {
  int const nElem = m_mesh_.getNumberOfElements();
  int n_dg = 0;
  int n_sem = 0;

  if (m_external_element_type_.size() == static_cast<size_t>(nElem)) {
    // Caller-provided split (setElementTags()): skip the Z-threshold heuristic entirely, since
    // it only cuts the intended plane while the mesh is flat.
    for (int e = 0; e < nElem; ++e) {
      m_element_type_[e] = m_external_element_type_[e];
      if (m_element_type_[e] == kElementTypeDG)
        ++n_dg;
      else
        ++n_sem;
    }
  } else {
    int const mid = ORDER / 2;
    for (int e = 0; e < nElem; ++e) {
      int const gIdx = m_mesh_.globalNodeIndex(e, mid, mid, mid);
      float const zCoord = m_mesh_.nodeCoord(gIdx, 2);
      if (zCoord < DG_SEM_interface_z_) {
        m_element_type_[e] = kElementTypeDG;
        ++n_dg;
      } else {
        m_element_type_[e] = kElementTypeSEM;
        ++n_sem;
      }
    }
  }

  num_DG_elements_ = n_dg;
  num_SEm_elements_ = n_sem;

  DG_elem_list_ = allocateVector<vectorInt>(num_DG_elements_, "DGElemList");
  SEm_elem_list_ = allocateVector<vectorInt>(num_SEm_elements_, "SEmElemList");
  int idg = 0;
  int isem = 0;
  for (int e = 0; e < nElem; ++e) {
    if (m_element_type_[e] == kElementTypeDG)
      DG_elem_list_[idg++] = e;
    else
      SEm_elem_list_[isem++] = e;
  }
}

//============================================================================
// TagNodes
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::TagNodes() {
  int const nNode = m_mesh_.getNumberOfNodes();
  int const nElem = m_mesh_.getNumberOfElements();
  int const dim = ORDER + 1;

  vectorInt dg_count = allocateVector<vectorInt>(nNode, "dgCount");
  vectorInt sem_count = allocateVector<vectorInt>(nNode, "semCount");

  Kokkos::parallel_for(
      "TagNodes_initCount", nNode, KOKKOS_LAMBDA(const int i) {
        dg_count[i] = 0;
        sem_count[i] = 0;
      });
  FENCE

  auto elem_type = m_element_type_;
  auto mesh_local = m_mesh_;

  Kokkos::parallel_for(
      "TagNodes_mainLoop", nElem, KOKKOS_LAMBDA(const int e) {
        if (e >= nElem) return;

        int const etype = elem_type[e];
        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const gIdx = mesh_local.globalNodeIndex(e, i, j, k);
              if (etype == kElementTypeDG) {
                ATOMICADD(dg_count[gIdx], 1);
              } else {
                ATOMICADD(sem_count[gIdx], 1);
              }
            }
      });
  FENCE

  // Host mirrors: dg_count/sem_count are filled by device kernel; need host access for face loop.
  auto h_dg_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, dg_count);
  auto h_sem_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, sem_count);

  // Iterate over m_face_connectivity_ faces — same index space as coupling kernels.
  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());
  int n_interface = 0;
  for (int f = 0; f < num_faces_fc; ++f) {
    if (m_face_connectivity_.isBoundaryFace(f)) continue;
    bool face_on_interface = true;
    for (int j = 0; j < knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_dg_count(gn) == 0 || h_sem_count(gn) == 0) {
        face_on_interface = false;
        break;
      }
    }
    if (face_on_interface) ++n_interface;
  }

  num_interface_faces_ = n_interface;
  m_interface_face_indices_ = allocateVector<vectorInt>(num_interface_faces_, "interfaceFaceIndices");

  int idx = 0;
  for (int f = 0; f < num_faces_fc; ++f) {
    if (m_face_connectivity_.isBoundaryFace(f)) continue;
    bool face_on_interface = true;
    for (int j = 0; j < knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_dg_count(gn) == 0 || h_sem_count(gn) == 0) {
        face_on_interface = false;
        break;
      }
    }
    if (face_on_interface) m_interface_face_indices_[idx++] = f;
  }

  // Build compact SEM node list.
  {
    int n_sem = 0;
    for (int n = 0; n < nNode; ++n)
      if (h_sem_count(n) > 0) ++n_sem;

    num_SEm_nodes_ = n_sem;
    SEm_node_list_ = allocateVector<vectorInt>(n_sem, "SEmNodeList");

    int isem = 0;
    for (int n = 0; n < nNode; ++n)
      if (h_sem_count(n) > 0) SEm_node_list_[isem++] = n;
  }

  BuildDGInteriorFaceList();
  std::cout << "DGSEMsolver: " << m_n_DG_interior_faces_ << " DG interior faces." << std::endl;
}

//============================================================================
// BuildDGInteriorFaceList
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::BuildDGInteriorFaceList() {
  auto h_elem_type = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_element_type_);
  auto h_iface = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_interface_face_indices_);

  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());

  std::vector<bool> is_iface(num_faces_fc, false);
  for (int i = 0; i < num_interface_faces_; ++i) is_iface[h_iface(i)] = true;

  // Collect faces adjacent to at least one DG element and not on the DG-SEM interface.
  std::vector<int> result;
  result.reserve(num_faces_fc / 2);
  // Collect ALL faces adjacent to at least one DG element (interior + interface).
  // Used by the fused face kernel (DG-DG flux + DG-SEM coupling in one launch).
  std::vector<int> result_all;
  result_all.reserve(num_faces_fc / 2);
  for (int f = 0; f < num_faces_fc; ++f) {
    if (is_iface[f]) continue;
    int const oe = m_face_connectivity_.elemOwner(f);
    bool dg_adj = (h_elem_type(oe) == kElementTypeDG);
    if (dg_adj) result.push_back(f);
  }
  for (int f = 0; f < num_faces_fc; ++f) {
    int const oe = m_face_connectivity_.elemOwner(f);
    bool dg_adj = (h_elem_type(oe) == kElementTypeDG);
    if (dg_adj) result_all.push_back(f);
  }

  m_n_DG_interior_faces_ = static_cast<int>(result.size());
  m_DG_interior_face_list_ = allocateVector<vectorInt>(m_n_DG_interior_faces_, "DGInteriorFaceList");
  auto h_list = Kokkos::create_mirror_view(m_DG_interior_face_list_);
  for (int i = 0; i < m_n_DG_interior_faces_; ++i) h_list(i) = result[i];
  Kokkos::deep_copy(m_DG_interior_face_list_, h_list);

  m_n_DG_all_faces_ = static_cast<int>(result_all.size());
  m_DG_all_face_list_ = allocateVector<vectorInt>(m_n_DG_all_faces_, "DGAllFaceList");
  auto h_all = Kokkos::create_mirror_view(m_DG_all_face_list_);
  for (int i = 0; i < m_n_DG_all_faces_; ++i) h_all(i) = result_all[i];
  Kokkos::deep_copy(m_DG_all_face_list_, h_all);
}

//============================================================================
// applyVerletFused - single-launch Verlet update for both sub-domains
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::applyVerletFused(const float& dt,
                                                                                                const DataType& data) {
  float const dt_local = dt;
  float const dt2_local = dt * dt;

  // DG side: element-wise Verlet over the DG element list.
  auto mesh_local = m_mesh_;
  auto const p_DG_cur = data.m_wavefield.m_DGacoustic.getCurrentField(0);
  auto const p_DG_prev = data.m_wavefield.m_DGacoustic.getPreviousField(0);
  arrayReal mass_dg = m_DG_solver_.m_mass_local_;
  arrayReal stiff_dg = m_DG_solver_.m_stiff_local_;
  arrayReal damp_dg = m_DG_solver_.m_damp_local_;
  arrayReal rhs_dg = m_DG_solver_.m_rhs_elem_;
  auto dg_elem_list = DG_elem_list_;
  int const n_dg = num_DG_elements_;
  constexpr int kPPE = dgSolver::kPointsPerElement;

  // SEM side: node-wise Verlet over the SEM node list.
  auto const p_SEM_cur = data.m_wavefield.m_SEMacoustic.getCurrentField(0);
  auto const p_SEM_prev = data.m_wavefield.m_SEMacoustic.getPreviousField(0);
  vectorReal mass_sem = m_SEm_solver_.getMassMatrixAcoustic();
  vectorReal damp_sem = m_SEm_solver_.getDampingMatrix(0);
  vectorReal work_sem = m_SEm_solver_.getForceVector(0);
  vectorReal taper_sem = m_SEm_solver_.getSpongeTaperCoeff();
  auto sem_node_list = SEm_node_list_;
  int const n_sem = num_SEm_nodes_;

  int const n_total = n_dg + n_sem;
  Kokkos::parallel_for(
      "DG-SEM Fused Verlet", n_total, KOKKOS_LAMBDA(const int idx) {
        if (idx < n_dg) {
          int const e = dg_elem_list[idx];
          for (int i = 0; i < kPPE; ++i) {
            float const M = mass_dg(e, i);
            float const K = stiff_dg(e, i) + rhs_dg(e, i);
            float const D = damp_dg(e, i);
            p_DG_prev(e, i) =
                (2.0f * M * p_DG_cur(e, i) - dt2_local * K - (M - 0.5f * dt_local * D) * p_DG_prev(e, i)) /
                (M + 0.5f * dt_local * D);
          }
        } else {
          int const I = sem_node_list[idx - n_dg];
          if (mass_sem(I) <= 0.0f) return;
          if (mesh_local.isFreeSurface(I)) {
            p_SEM_cur(I) = 0.0f;
            p_SEM_prev(I) = 0.0f;
          } else {
            float next_val = (2.0f * mass_sem(I) * p_SEM_cur(I) -
                              (mass_sem(I) - 0.5f * dt_local * damp_sem(I)) * p_SEM_prev(I) -
                              dt2_local * work_sem(I));
            p_SEM_prev(I) = next_val / (mass_sem(I) + 0.5f * dt_local * damp_sem(I));
            p_SEM_prev(I) *= taper_sem(I);
            p_SEM_cur(I) *= taper_sem(I);
          }
        }
      });
}

//============================================================================
// computeOneStep  (staggered DG-SEM coupling scheme)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeOneStep(const float& dt,
                                                                                              const int& timeSample,
                                                                                              DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  int const nNode = m_mesh_.getNumberOfNodes();

  if (myData.isDistributed) {
    throw std::runtime_error(
        "computeOneStep called in distributed mode. Use computeForces() -> "
        "synchronize() -> updateSolutionForward().");
  }

  // Sub-solver data views are constructed once and reused throughout the step.
  DGsolverDataAcoustic DG_data(myData.m_wavefield.m_DGacoustic, myData.m_rhs.m_rhs_DGacoustic);

  SEMsolverData<utils::enums::physicType::kAcoustic> SEm_data(myData.m_wavefield.m_SEMacoustic,
                                                              myData.m_rhs.m_rhs_SEMacoustic);

  // =========================================================================
  // DG: volume + DG-DG interior flux + DG-SEM interface coupling (fused face kernel)
  // =========================================================================

  m_DG_solver_.m_list_mode_ = true;
  m_DG_solver_.m_elem_list_ = DG_elem_list_;
  m_DG_solver_.m_n_elem_list_ = num_DG_elements_;
  m_DG_solver_.m_face_list_ = m_DG_all_face_list_;
  m_DG_solver_.m_n_face_list_ = m_n_DG_all_faces_;

  m_DG_solver_.computeVolumeAndBoundary(num_DG_elements_, DG_data.getCurrentField(0), timeSample, DG_data);

  // =========================================================================
  // SEM: source + stiffness (Neumann = 0 at interface until coupling kernel).
  // Must run BEFORE the fused face kernel: the coupling branch writes work_sem,
  // and resetGlobalVectors zeroes it — so the reset must precede the coupling.
  m_SEm_solver_.resetGlobalVectors(nNode);
  // Fold the SEM source term into the element-contribution kernel (removes the
  // separate applyRHSTerm launch). Disabled again after the step.
  m_SEm_solver_.setRhsTimeSample(timeSample);
  m_SEm_solver_.computeElementContributionsFromList(SEm_data, SEm_elem_list_, num_SEm_elements_);
  m_SEm_solver_.setRhsTimeSample(-1);

  // Fused face kernel: DG-DG interior flux + DG-SEM interface coupling in one launch.
  // Interface faces (owner/neighbor in different domains) take the SIPG coupling branch,
  // writing both stiff_dg and work_sem (after the SEM reset, so the coupling survives).
  {
    auto const p_sem = myData.m_wavefield.m_SEMacoustic.getCurrentField(0);
    auto work_sem = m_SEm_solver_.getForceVector(0);
    m_DG_solver_.computeBoundaryDampingAndInterfaceFlux(m_n_DG_all_faces_, DG_data.getCurrentField(0), &p_sem,
                                                        &work_sem, &m_element_type_);
  }

  // =========================================================================
  // Both Verlots — fused into a single launch (DG element-wise + SEM node-wise).
  // =========================================================================

  m_DG_solver_.m_list_mode_ = false;
  applyVerletFused(dt, myData);
  // Final fence: synchronizes the GPU with the host so the caller (unit test,
  // benchmark loop) reads a completed step. Removing it would make time_s
  // measure kernel-launch time instead of GPU execution (illusory speedup).
  FENCE
}

//============================================================================
// outputSolutionValues (SEM solution output: p[nNodes])
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const vectorReal& field, const char* fieldName) {
  m_SEm_solver_.outputSolutionValues(t, e, field, fieldName);
}

//============================================================================
// outputSolutionValues (DG solution output: p[nElem][nDof])
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const arrayReal& field, const char* fieldName) {
  m_DG_solver_.outputSolutionValues(t, e, field, fieldName);
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_
