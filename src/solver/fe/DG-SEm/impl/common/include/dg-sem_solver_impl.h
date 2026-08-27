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
  // Hand our connectivity to the DG sub-solver before its own init: the DG
  // interior face list built below (BuildDGInteriorFaceList) holds face ids in
  // this numbering and is fed to the DG face kernels through m_face_list_.
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
  for (int f = 0; f < num_faces_fc; ++f) {
    if (is_iface[f]) continue;
    int const oe = m_face_connectivity_.elemOwner(f);
    bool dg_adj = (h_elem_type(oe) == kElementTypeDG);
    if (dg_adj) result.push_back(f);
  }

  m_n_DG_interior_faces_ = static_cast<int>(result.size());
  m_DG_interior_face_list_ = allocateVector<vectorInt>(m_n_DG_interior_faces_, "DGInteriorFaceList");
  auto h_list = Kokkos::create_mirror_view(m_DG_interior_face_list_);
  for (int i = 0; i < m_n_DG_interior_faces_; ++i) h_list(i) = result[i];
  Kokkos::deep_copy(m_DG_interior_face_list_, h_list);
}

//============================================================================
// ApplyCoupling — SIPG flux: SEM pressure → DG stiff_local_ and DG stiff_local_ → SEM pressure
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES,
          utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::ApplyCoupling(const DataType& data) {
  auto mesh_local = m_mesh_;
  auto face_connectivity_local = m_face_connectivity_;
  auto const p_DG = data.m_wavefield.m_DGacoustic.getCurrentField(0);
  auto const p_SEM = data.m_wavefield.m_SEMacoustic.getCurrentField(0);

  auto iface_list = m_interface_face_indices_;
  int const n_iface = num_interface_faces_;
  auto element_type = m_element_type_;
  vectorReal work_sem = m_SEm_solver_.getForceVector(0);
  arrayReal stiff_dg = m_DG_solver_.m_stiff_local_;
  auto const face_to_elem_dof = dgSolver::kFaceToElemDof;
  real_t const penalty_local = m_penalty_factor_;

  Kokkos::parallel_for(
      "ApplyCouplingSEMToDG", n_iface, KOKKOS_LAMBDA(const int _loop_idx) {
        int const f = iface_list(_loop_idx);

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        bool const owner_is_dg = (element_type(owner_e) == kElementTypeDG);
        int const dg_e = owner_is_dg ? owner_e : neighbor_e;
        int const sem_e = owner_is_dg ? neighbor_e : owner_e;
        int const fid_dg = owner_is_dg ? fid_o : fid_n;
        int const fid_sem = owner_is_dg ? fid_n : fid_o;

        auto dg_to_sem = [&](int i) {
          return owner_is_dg ? face_connectivity_local.getNeighborFaceDof(f, i)
                             : face_connectivity_local.getOwnerFaceDof(f, i);
        };
        auto sem_to_dg = [&](int i) {
          return owner_is_dg ? face_connectivity_local.getOwnerFaceDof(f, i)
                             : face_connectivity_local.getNeighborFaceDof(f, i);
        };

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni = face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float dg_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(dg_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        dg_coords[iv + 2 * jv + 4 * kv]);
        }

        float sem_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(sem_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        sem_coords[iv + 2 * jv + 4 * kv]);
        }

        float normal_dg[3];
        mesh_local.faceNormal(dg_e, static_cast<model::CubicFace>(fid_dg), normal_dg);

        real_t const inv_rho_dg = 1.0f / mesh_local.getModelRhoOnElement(dg_e);
        real_t const gamma_dg = computeSIPGPenalty<ORDER>(faceCoords, dg_coords, penalty_local);
        real_t const inv_rho_sem = 1.0f / mesh_local.getModelRhoOnElement(sem_e);
        real_t const gamma_sem = computeSIPGPenalty<ORDER>(faceCoords, sem_coords, penalty_local);

        // Face-sized accumulator indexed by DG-side face dof: the coupling flux only touches
        // the shared face's (ORDER+1)^2 dofs, so an element-sized array forced a 7x larger
        // per-thread local-memory footprint and an all-element atomic flush of mostly zeros.
        float stiff_dg_local[knumNodesPerFace] = {0};

        INTEGRAL_TYPE::computeInterfaceFluxTerm(
            faceCoords, dg_coords, fid_dg, [&](const int i, const int j, const int k, const real_t val) {
              int const ei = face_to_elem_dof[fid_dg][i];
              int const ej = face_to_elem_dof[fid_dg][j];
              int const gn_j = face_connectivity_local.getGlobalNodeFromFace(f, dg_to_sem(j));
              float const nk = normal_dg[k];
              stiff_dg_local[i] += inv_rho_dg * nk * (-0.5f * val * p_DG(dg_e, ej) + 0.5f * val * p_SEM(gn_j));
              stiff_dg_local[j] += inv_rho_dg * nk * (-0.5f * val * p_DG(dg_e, ei));
              ATOMICADD(work_sem(gn_j), inv_rho_dg * nk * (0.5f * val * p_DG(dg_e, ei)));
            });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const ei = face_to_elem_dof[fid_dg][i];
          int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, dg_to_sem(i));
          stiff_dg_local[i] +=
              gamma_dg * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) * (p_DG(dg_e, ei) - p_SEM(gn_i));
        }

        INTEGRAL_TYPE::computeInterfaceFluxTerm(
            faceCoords, sem_coords, fid_sem, [&](const int i, const int j, const int k, const real_t val) {
              int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, i);
              int const gn_j = face_connectivity_local.getGlobalNodeFromFace(f, j);
              int const sd_j = sem_to_dg(j);
              int const ej_perm = face_to_elem_dof[fid_dg][sd_j];
              float const nk = -normal_dg[k];  // SEM outward = -DG outward
              stiff_dg_local[sd_j] += inv_rho_sem * nk * (0.5f * val * p_SEM(gn_i));
              ATOMICADD(work_sem(gn_i),
                        inv_rho_sem * nk * (-0.5f * val * p_SEM(gn_j) + 0.5f * val * p_DG(dg_e, ej_perm)));
              ATOMICADD(work_sem(gn_j), inv_rho_sem * nk * (-0.5f * val * p_SEM(gn_i)));
            });

        for (int i = 0; i < knumNodesPerFace; ++i) {
          int const gn_i = face_connectivity_local.getGlobalNodeFromFace(f, i);
          int const ei_perm = face_to_elem_dof[fid_dg][sem_to_dg(i)];
          ATOMICADD(work_sem(gn_i), inv_rho_sem * gamma_sem * INTEGRAL_TYPE::computeDampingTerm(i, faceCoords) *
                                        (p_SEM(gn_i) - p_DG(dg_e, ei_perm)));
        }

        for (int i = 0; i < knumNodesPerFace; ++i)
          ATOMICADD(stiff_dg(dg_e, face_to_elem_dof[fid_dg][i]), stiff_dg_local[i]);
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
  // DG: volume + DG-DG interior flux (interface faces excluded from face list)
  // =========================================================================

  m_DG_solver_.m_list_mode_ = true;
  m_DG_solver_.m_elem_list_ = DG_elem_list_;
  m_DG_solver_.m_n_elem_list_ = num_DG_elements_;
  m_DG_solver_.m_face_list_ = m_DG_interior_face_list_;
  m_DG_solver_.m_n_face_list_ = m_n_DG_interior_faces_;

  m_DG_solver_.applyRHSTerm(timeSample, dt, DG_data);
  FENCE
  m_DG_solver_.computeVolumeAndBoundary(num_DG_elements_, DG_data.getCurrentField(0));
  FENCE
  m_DG_solver_.computeBoundaryDampingAndInterfaceFlux(m_n_DG_interior_faces_, DG_data.getCurrentField(0));
  FENCE

  // =========================================================================
  // SEM: source + stiffness (Neumann = 0 at interface until coupling kernel)
  // =========================================================================

  m_SEm_solver_.resetGlobalVectors(nNode);
  FENCE
  m_SEm_solver_.applyRHSTerm(timeSample, dt, SEm_data);
  FENCE
  m_SEm_solver_.computeElementContributionsFromList(SEm_data, SEm_elem_list_, num_SEm_elements_);
  FENCE

  // =========================================================================
  // Symmetric SIPG interface coupling: both sides read p^n (no temporal lag).
  // =========================================================================

  ApplyCoupling(myData);
  FENCE

  // =========================================================================
  // Both Verlots
  // =========================================================================

  m_DG_solver_.applyVerlet(num_DG_elements_, dt, DG_data.getCurrentField(0), DG_data.getPreviousField(0));
  m_DG_solver_.m_list_mode_ = false;
  FENCE

  m_SEm_solver_.updateFieldsFromListForward(dt, SEm_data, SEm_node_list_, num_SEm_nodes_);
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
