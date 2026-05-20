#ifndef FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "Integrals.h"
#include "data_type.h"
#include "dg-sem_solver.h"
#include "dg-sem_solver_data.h"

namespace solver {
namespace fe {

//============================================================================
// computeFEInit
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeFEInit(
    model::ModelApi<float, int>& mesh_in, const std::array<float, 3>& sponge_size, const bool surface_sponge,
    const float taper_delta) {
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh_ = *typed;
  } else {
    throw std::runtime_error("DGSEMsolver: incompatible mesh type in computeFEInit");
  }

  // Initialise sub-solvers (mass/damping matrices are overridden below).
  m_SEm_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  m_DG_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);

  allocateFEarrays();

  TagElements();
  std::cout << "DGSEMsolver: " << num_SEm_elements_ << " SEm elements, " << num_DG_elements_
            << " DG elements." << std::endl;

  TagNodes();
  std::cout << "DGSEMsolver: " << num_interface_faces_ << " interface faces." << std::endl;
}

//============================================================================
// allocateFEarrays
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::allocateFEarrays() {
  int const nElem = m_mesh_.getNumberOfElements();
  int const nNode = m_mesh_.getNumberOfNodes();

  m_element_type_ = allocateVector<vectorInt>(nElem, "DGSEMElementType");
}

//============================================================================
// TagElements
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::TagElements() {
  int const nElem = m_mesh_.getNumberOfElements();
  int n_dg = 0;
  int n_sem = 0;
  
  int const mid = ORDER / 2;
  for (int e = 0; e < nElem; ++e) {
    int const gIdx = m_mesh_.globalNodeIndex(e, mid, mid, mid);
    int const zCoord = m_mesh_.nodeCoord(gIdx, 2);
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

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
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

  int const num_faces_mesh = mesh_local.getNumberOfFaces();
  int n_interface = 0;
  for (int f = 0; f < num_faces_mesh; ++f) {
    std::vector<int> face_nodes(knumNodesPerFace);
    for (int j = 0; j < knumNodesPerFace; ++j) {
      face_nodes[j] = mesh_local.getGlobalNodeFromFace(f, j);
    }
    const bool face_on_interface = std::all_of(face_nodes.begin(), face_nodes.end(),
        [&](int n) { return dg_count[n] > 0 && sem_count[n] > 0; }
    );
    if (face_on_interface) ++n_interface;
  }
  
  num_interface_faces_ = n_interface;
  m_interface_face_indices_ = allocateVector<vectorInt>(num_interface_faces_, "interfaceFaceIndices");

  int idx = 0;
  for (int f = 0; f < num_faces_mesh; ++f) {
    std::vector<int> face_nodes(knumNodesPerFace);
    for (int j = 0; j < knumNodesPerFace; ++j) {
      face_nodes[j] = mesh_local.getGlobalNodeFromFace(f, j);
    }
    const bool face_on_interface = std::all_of(face_nodes.begin(), face_nodes.end(),
        [&](int n) { return dg_count[n] > 0 && sem_count[n] > 0; }
    );
    if (face_on_interface) m_interface_face_indices_[++idx] = f;
  }

  // Build compact SEM node lists for updateFieldsFromList.
  {
    int n_sem = 0;
    for (int n = 0; n < nNode; ++n) {
      if (sem_count[n] > 0) ++n_sem;
    }

    num_SEm_nodes_ = n_sem;
    SEm_node_list_ = allocateVector<vectorInt>(n_sem, "SEmNodeList");

    int isem = 0;
    for (int n = 0; n < nNode; ++n) {
      if (sem_count[n] > 0) SEm_node_list_[isem++] = n;
    }
  }

}

//============================================================================
// ApplyCouplingDGToSEM
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::ApplyCouplingDGToSEM(
    float dt, const DataType& data) {
  float const dt2 = dt * dt;
  auto mesh_local = m_mesh_;
  auto face_connectivity_local = m_face_connectivity_;

  auto p_SEM_next = data.m_wavefield.m_SEMacoustic.getPreviousField(0);
  auto p_DG_next = data.m_wavefield.m_DGacoustic.getPreviousField(0);
  auto iface_list = m_interface_face_indices_;
  int const n_iface = num_interface_faces_;
  auto element_type = m_element_type_;

  Kokkos::parallel_for(
      "ApplyCouplingDGToSEM_Loop", n_iface, KOKKOS_LAMBDA(const int _loop_idx) {
        int const f = iface_list(_loop_idx);

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        bool const owner_is_sem = element_type(owner_e) == kElementTypeSEM;
        int const sem_e = owner_is_sem ? owner_e : neighbor_e;
        int const dg_e = owner_is_sem ? neighbor_e : owner_e;
        int const fid_sem = owner_is_sem ? fid_o : fid_n;
        int const fid_dg = owner_is_sem ? fid_n : fid_o;

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni = face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float sem_coords[8][3];
        auto const eIdx_o = mesh_local.elementIndex(sem_e);
        for (int kv = 0; kv < 2; ++kv)
          for (int jv = 0; jv < 2; ++jv)
            for (int iv = 0; iv < 2; ++iv)
              mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx_o, iv, jv, kv),
                                      sem_coords[iv + 2 * jv + 4 * kv]);

        float normal[3];
        mesh_local.faceNormal(sem_e, static_cast<model::CubicFace>(fid_sem), normal);

        real_t const inv_rho = 1.0f / mesh_local.getModelRhoOnElement(sem_e);

        INTEGRAL_TYPE::computeInterfaceFluxTerm(
            faceCoords, sem_coords, fid_sem, [&](const int i, const int j, const int k, const real_t val) {
              int const global_node = face_connectivity_local.getGlobalNodeFromFace(f, j);
              int const ei_perm = faceLocalToElemLocal(static_cast<model::CubicFace>(fid_dg), face_connectivity_local.getNeighborFaceDof(f, i), ORDER);
              float const nk = normal[k];
              ATOMICADD(p_SEM_next(global_node), dt2 * inv_rho * 0.5f * val * nk * p_DG_next(dg_e, ei_perm));
            });
      });
}

//============================================================================
// ApplyCouplingSEMToDG
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::ApplyCouplingSEMToDG(
    const DataType& data) {
  auto face_connectivity_local = m_face_connectivity_;
  auto p_DG_next = data.m_wavefield.m_DGacoustic.getPreviousField(0);
  auto p_SEM_curr = data.m_wavefield.m_SEMacoustic.getCurrentField(0);
  
  auto iface_list = m_interface_face_indices_;
  int const n_iface = num_interface_faces_;
  auto element_type = m_element_type_;

  Kokkos::parallel_for(
      "ApplyCouplingSEMToDG_Loop", n_iface, KOKKOS_LAMBDA(const int _loop_idx) {
        int const f = iface_list(_loop_idx);

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        bool const owner_is_sem = element_type(owner_e) == kElementTypeSEM;
        int const sem_e = owner_is_sem ? owner_e : neighbor_e;
        int const fid_sem = owner_is_sem ? fid_o : fid_n;
        
        for (int j = 0; j < knumNodesPerFace; ++j) {
          int const global_node = face_connectivity_local.getGlobalNodeFromFace(f, j);
          int const elem_dof = faceLocalToElemLocal(static_cast<model::CubicFace>(fid_sem), j, ORDER);
          p_DG_next(sem_e, elem_dof) = p_SEM_curr(global_node);
        }
      });
}

//============================================================================
// computeOneStep  (staggered DG-SEM coupling scheme)
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeOneStep(const float& dt,
                                                                                                 const int& timeSample,
                                                                                                 DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  int const nNode = m_mesh_.getNumberOfNodes();

  // Sub-solver data views are constructed once and reused throughout the step.
  DGsolverDataAcoustic DG_data(myData.m_wavefield.m_DGacoustic, myData.m_rhs.m_rhs_DGacoustic);
  
  SEMsolverData<utils::enums::physicType::kAcoustic> SEm_data(myData.m_wavefield.m_SEMacoustic,
                                                                   myData.m_rhs.m_rhs_SEMacoustic);

  // =========================================================================
  // DG STEP
  // =========================================================================

  // 1. Apply DG source term.
  m_DG_solver_.computeForces(dt, timeSample, DG_data);
  FENCE

  // 2. DG Verlet: u^{n+1} written into DG_data.getPreviousField().
  m_DG_solver_.updateFieldsFromList(dt, DG_data, DG_elem_list_, num_DG_elements_);
  FENCE

  // 3. SEM→DG coupling.
  ApplyCouplingSEMToDG(myData);
  FENCE

  // =========================================================================
  // SEM STEP
  // =========================================================================

  // 4. Reset SEM work vector.
  m_SEm_solver_.resetGlobalVectors(nNode);
  FENCE

  // 5. Apply SEM source term.
  m_SEm_solver_.applyRHSTerm(timeSample, dt, SEm_data);
  FENCE

  // 6. Compute SEM stiffness on the right list of elements.
  m_SEm_solver_.computeElementContributionsFromList(SEm_data, SEm_elem_list_, num_SEm_elements_);
  FENCE

  // 7. SEM Verlet: p^{n+1} written into SEm_data.getPreviousField().
  m_SEm_solver_.updateFieldsFromList(dt, SEm_data, SEm_node_list_, num_SEm_nodes_);
  FENCE

  // 8. DG→SEM coupling post-Verlet.
  ApplyCouplingDGToSEM(dt, myData);
  FENCE
}

//============================================================================
// outputSolutionValues (SEM solution output: p[nNodes])
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const vectorReal& field, const char* fieldName) {
  m_SEm_solver_.outputSolutionValues(t, e, field, fieldName);
}

//============================================================================
// outputSolutionValues (DG solution output: p[nElem][nDof])
//============================================================================

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGSEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const arrayReal& field, const char* fieldName) {
  m_DG_solver_.outputSolutionValues(t, e, field, fieldName);
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_SEM_IMPL_COMMON_INCLUDE_DG_SEM_SOLVER_IMPL_H_
