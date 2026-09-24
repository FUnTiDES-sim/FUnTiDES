#ifndef FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_IMPL_H_
#define FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_IMPL_H_

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Integrals.h"
#include "dg_padaptive_solver.h"
#include "dg_padaptive_solver_data.h"
#include "dg_solver.h"
#include "dg_solver_impl.h"

namespace solver {
namespace fe {

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::computeFEInit(model::ModelApi<float, int>& mesh_in,
                                               const std::array<float, 3>& sponge_size, const bool surface_sponge,
                                               const float taper_delta) {
  if (auto* typed = dynamic_cast<MESH_TYPE*>(&mesh_in)) {
    m_mesh_ = *typed;
  } else {
    throw std::runtime_error("DGPAdaptiveSolver: incompatible mesh type in computeFEInit");
  }

  m_face_connectivity_.build(m_mesh_);

  m_pMin_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);
  m_pMax_solver_.computeFEInit(mesh_in, sponge_size, surface_sponge, taper_delta);

  // The pMin sub-solver built its face connectivity at ORDER_MIN, but the shared mesh has order
  // ORDER_MAX. Without this rebuild, its interior "Plus" faces (kXPlus/kYPlus/kZPlus) fail the
  // neighbor node-id matching, because the face-normal coordinate stays at ORDER_MIN instead of
  // the far edge of the mesh. No-op for pMax, whose order already equals the mesh order.
  if constexpr (ORDER_MIN != ORDER_MAX) {
    m_pMin_solver_.rebuildFaceConnectivityGeometry(ORDER_MAX);
  }

  m_penalty_factor_ = m_pMax_solver_.getPenaltyFactor();  // both solvers have the same penalty factor

  allocateFEarrays();

  ComputeMortarProjection();

  std::cout << "DGPAdaptiveSolver: ORDER_MIN=" << ORDER_MIN << ", ORDER_MAX=" << ORDER_MAX << std::endl;

  TagElements();
  std::cout << "DGPAdaptiveSolver: " << num_pMin_elements_ << " pMin elements, " << num_pMax_elements_
            << " pMax elements." << std::endl;

  TagNodes();
  std::cout << "DGPAdaptiveSolver: " << num_interface_faces_ << " interface faces." << std::endl;
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::allocateFEarrays() {
  int const nElem = m_mesh_.getNumberOfElements();
  m_element_type_ = allocateVector<vectorInt>(nElem, "pMinpMaxElementType");
  m_p1d_projection_ = allocateArray2D<arrayReal>(ORDER_MAX + 1, ORDER_MIN + 1, "p1dProjectionMatrix");
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::TagElements() {
  int const nElem = m_mesh_.getNumberOfElements();
  int n_pMin = 0;
  int n_pMax = 0;

  if (m_external_element_type_.size() == static_cast<size_t>(nElem)) {
    // Caller-provided split (setElementTags()): the Z-threshold heuristic is skipped, since
    // it only cuts the intended plane while the mesh is flat.
    for (int e = 0; e < nElem; ++e) {
      m_element_type_[e] = m_external_element_type_[e];
      if (m_element_type_[e] == kElementTypePMin)
        ++n_pMin;
      else
        ++n_pMax;
    }
  } else {
    for (int e = 0; e < nElem; ++e) {
      int const gIdx = m_mesh_.globalNodeIndex(e, 0, 0, 0);
      float const zCoord = m_mesh_.nodeCoord(gIdx, 2);
      if (zCoord < pAdaptive_interface_z_) {
        m_element_type_[e] = kElementTypePMin;
        ++n_pMin;
      } else {
        m_element_type_[e] = kElementTypePMax;
        ++n_pMax;
      }
    }
  }

  num_pMin_elements_ = n_pMin;
  num_pMax_elements_ = n_pMax;

  pMin_elem_list_ = allocateVector<vectorInt>(num_pMin_elements_, "pMinElemList");
  pMax_elem_list_ = allocateVector<vectorInt>(num_pMax_elements_, "pMaxElemList");
  int ipMin = 0;
  int ipMax = 0;
  for (int e = 0; e < nElem; ++e) {
    if (m_element_type_[e] == kElementTypePMin)
      pMin_elem_list_[ipMin++] = e;
    else
      pMax_elem_list_[ipMax++] = e;
  }
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::TagNodes() {
  int const nNode = m_mesh_.getNumberOfNodes();
  int const nElem = m_mesh_.getNumberOfElements();
  int const dim = ORDER_MAX + 1;

  vectorInt pMin_count = allocateVector<vectorInt>(nNode, "pMinCount");
  vectorInt pMax_count = allocateVector<vectorInt>(nNode, "pMaxCount");

  Kokkos::parallel_for(
      "TagNodes_initCount", nNode, KOKKOS_LAMBDA(const int i) {
        pMin_count[i] = 0;
        pMax_count[i] = 0;
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
              if (etype == kElementTypePMin) {
                ATOMICADD(pMin_count[gIdx], 1);
              } else {
                ATOMICADD(pMax_count[gIdx], 1);
              }
            }
      });
  FENCE

  // The counts are filled on the device; the face loops below run on the host.
  auto h_pMin_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, pMin_count);
  auto h_pMax_count = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, pMax_count);

  // Iterate over the faces of m_face_connectivity_, the same index space as the coupling kernels.
  // A face is on the interface when all its nodes belong to both a pMin and a pMax element.
  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());
  int n_interface = 0;
  for (int f = 0; f < num_faces_fc; ++f) {
    if (m_face_connectivity_.isBoundaryFace(f)) continue;
    bool face_on_interface = true;
    for (int j = 0; j < pMaxSolver::knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_pMin_count(gn) == 0 || h_pMax_count(gn) == 0) {
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
    for (int j = 0; j < pMaxSolver::knumNodesPerFace; ++j) {
      int const gn = m_face_connectivity_.getGlobalNodeFromFace(f, j);
      if (h_pMin_count(gn) == 0 || h_pMax_count(gn) == 0) {
        face_on_interface = false;
        break;
      }
    }
    if (face_on_interface) m_interface_face_indices_[idx++] = f;
  }

  BuildInteriorFaceLists();
  std::cout << "DGPAdaptiveSolver: " << m_n_pMin_interior_faces_ << " pMin interior faces." << std::endl;
  std::cout << "DGPAdaptiveSolver: " << m_n_pMax_interior_faces_ << " pMax interior faces." << std::endl;

  BuildInterfaceElementList();
  std::cout << "DGPAdaptiveSolver: " << m_n_iface_pMin_elements_ << " interface-adjacent pMin elements." << std::endl;
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::BuildInteriorFaceLists() {
  auto h_elem_type = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_element_type_);
  auto h_iface = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_interface_face_indices_);

  int const num_faces_fc = static_cast<int>(m_face_connectivity_.getNumberOfFaces());

  std::vector<bool> is_iface(num_faces_fc, false);
  for (int i = 0; i < num_interface_faces_; ++i) is_iface[h_iface(i)] = true;

  // Faces whose owner is a pMin element and which are not on the pMin-pMax interface.
  std::vector<int> result_pMin;
  result_pMin.reserve(num_faces_fc / 2);
  for (int f = 0; f < num_faces_fc; ++f) {
    if (is_iface[f]) continue;
    int const oe = m_face_connectivity_.elemOwner(f);
    bool pMin_adj = (h_elem_type(oe) == kElementTypePMin);
    if (pMin_adj) result_pMin.push_back(f);
  }

  m_n_pMin_interior_faces_ = static_cast<int>(result_pMin.size());
  m_pMin_interior_face_list_ = allocateVector<vectorInt>(m_n_pMin_interior_faces_, "pMinInteriorFaceList");
  auto h_pMin_list = Kokkos::create_mirror_view(m_pMin_interior_face_list_);
  for (int i = 0; i < m_n_pMin_interior_faces_; ++i) h_pMin_list(i) = result_pMin[i];
  Kokkos::deep_copy(m_pMin_interior_face_list_, h_pMin_list);

  // Faces whose owner is a pMax element and which are not on the pMin-pMax interface.
  std::vector<int> result_pMax;
  result_pMax.reserve(num_faces_fc / 2);
  for (int f = 0; f < num_faces_fc; ++f) {
    if (is_iface[f]) continue;
    int const oe = m_face_connectivity_.elemOwner(f);
    bool pMax_adj = (h_elem_type(oe) == kElementTypePMax);
    if (pMax_adj) result_pMax.push_back(f);
  }

  m_n_pMax_interior_faces_ = static_cast<int>(result_pMax.size());
  m_pMax_interior_face_list_ = allocateVector<vectorInt>(m_n_pMax_interior_faces_, "pMaxInteriorFaceList");
  auto h_pMax_list = Kokkos::create_mirror_view(m_pMax_interior_face_list_);
  for (int i = 0; i < m_n_pMax_interior_faces_; ++i) h_pMax_list(i) = result_pMax[i];
  Kokkos::deep_copy(m_pMax_interior_face_list_, h_pMax_list);
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::ComputeMortarProjection() {
  // Only the 1D factor is stored, of size (ORDER_MAX+1) x (ORDER_MIN+1): the projection applied by
  // the coupling is the threefold tensor product of this matrix, and
  // ProlongPMinField()/RestrictPMinStiff() apply it one direction at a time. Storing the 2D or 3D
  // form would cost (ORDER_MAX+1)^d (ORDER_MIN+1)^d entries and turn each application into a dense
  // product for no gain.
  for (int k = 0; k < ORDER_MAX + 1; ++k)
    for (int m = 0; m < ORDER_MIN + 1; ++m)
      m_p1d_projection_(k, m) =
          INTEGRAL_TYPE_MIN::BasisType::value(m, INTEGRAL_TYPE_MAX::BasisType::parentSupportCoord(k));
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::BuildInterfaceElementList() {
  int const nElem = m_mesh_.getNumberOfElements();
  auto h_elem_type = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_element_type_);
  auto h_iface = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, m_interface_face_indices_);

  std::vector<int> slot(nElem, -1);
  std::vector<int> elems;
  for (int i = 0; i < num_interface_faces_; ++i) {
    int const f = h_iface(i);
    int const sides[2] = {m_face_connectivity_.elemOwner(f), m_face_connectivity_.elemNeighbor(f)};
    for (int s = 0; s < 2; ++s) {
      int const e = sides[s];
      // An element facing the interface through several faces still needs raising only once.
      if (h_elem_type(e) != kElementTypePMin || slot[e] >= 0) continue;
      slot[e] = static_cast<int>(elems.size());
      elems.push_back(e);
    }
  }

  m_n_iface_pMin_elements_ = static_cast<int>(elems.size());
  m_iface_pMin_elem_list_ = allocateVector<vectorInt>(m_n_iface_pMin_elements_, "ifacePMinElemList");
  m_pMin_elem_to_slot_ = allocateVector<vectorInt>(nElem, "pMinElemToSlot");

  auto h_list = Kokkos::create_mirror_view(m_iface_pMin_elem_list_);
  auto h_slot = Kokkos::create_mirror_view(m_pMin_elem_to_slot_);
  for (int i = 0; i < m_n_iface_pMin_elements_; ++i) h_list(i) = elems[i];
  for (int e = 0; e < nElem; ++e) h_slot(e) = slot[e];
  Kokkos::deep_copy(m_iface_pMin_elem_list_, h_list);
  Kokkos::deep_copy(m_pMin_elem_to_slot_, h_slot);

  m_pMin_prolonged_field_ =
      allocateArray2D<arrayReal>(m_n_iface_pMin_elements_, pMaxSolver::kPointsPerElement, "pMinProlongedField");
  m_pMin_prolonged_stiff_ =
      allocateArray2D<arrayReal>(m_n_iface_pMin_elements_, pMaxSolver::kPointsPerElement, "pMinProlongedStiff");
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::ProlongPMinField(const DataType& data) {
  auto const pField_pMin = data.m_wavefield.m_pMinAcoustic.getCurrentField(0);
  auto elem_list = m_iface_pMin_elem_list_;
  arrayReal prolonged_field = m_pMin_prolonged_field_;
  arrayReal prolonged_stiff = m_pMin_prolonged_stiff_;
  arrayReal p1d = m_p1d_projection_;

  Kokkos::parallel_for(
      "ProlongPMinField", m_n_iface_pMin_elements_, KOKKOS_LAMBDA(const int slot) {
        constexpr int kNMin = kNumDofs1dMin;
        constexpr int kNMax = kNumDofs1dMax;
        int const e = elem_list(slot);

        // One direction at a time: kNMax*kNMin^3 + kNMax^2*kNMin^2 + kNMax^3*kNMin multiply-adds
        // instead of the (kNMax*kNMin)^3 a dense 3D matrix would cost.
        float tmp_i[kNMax * kNMin * kNMin];
        for (int c = 0; c < kNMin; ++c)
          for (int b = 0; b < kNMin; ++b)
            for (int a = 0; a < kNMax; ++a) {
              float acc = 0.0f;
              for (int i = 0; i < kNMin; ++i) acc += p1d(a, i) * pField_pMin(e, i + kNMin * (b + kNMin * c));
              tmp_i[a + kNMax * (b + kNMin * c)] = acc;
            }

        float tmp_ij[kNMax * kNMax * kNMin];
        for (int c = 0; c < kNMin; ++c)
          for (int b = 0; b < kNMax; ++b)
            for (int a = 0; a < kNMax; ++a) {
              float acc = 0.0f;
              for (int j = 0; j < kNMin; ++j) acc += p1d(b, j) * tmp_i[a + kNMax * (j + kNMin * c)];
              tmp_ij[a + kNMax * (b + kNMax * c)] = acc;
            }

        for (int c = 0; c < kNMax; ++c)
          for (int b = 0; b < kNMax; ++b)
            for (int a = 0; a < kNMax; ++a) {
              float acc = 0.0f;
              for (int k = 0; k < kNMin; ++k) acc += p1d(c, k) * tmp_ij[a + kNMax * (b + kNMax * k)];
              int const d = a + kNMax * (b + kNMax * c);
              prolonged_field(slot, d) = acc;
              prolonged_stiff(slot, d) = 0.0f;
            }
      });
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::RestrictPMinStiff() {
  auto elem_list = m_iface_pMin_elem_list_;
  arrayReal prolonged_stiff = m_pMin_prolonged_stiff_;
  arrayReal stiff_pMin = m_pMin_solver_.m_stiff_local_;
  arrayReal p1d = m_p1d_projection_;

  Kokkos::parallel_for(
      "RestrictPMinStiff", m_n_iface_pMin_elements_, KOKKOS_LAMBDA(const int slot) {
        constexpr int kNMin = kNumDofs1dMin;
        constexpr int kNMax = kNumDofs1dMax;
        int const e = elem_list(slot);

        // Transpose of ProlongPMinField(), contracted in the reverse order for the same reason.
        float tmp_i[kNMin * kNMax * kNMax];
        for (int c = 0; c < kNMax; ++c)
          for (int b = 0; b < kNMax; ++b)
            for (int i = 0; i < kNMin; ++i) {
              float acc = 0.0f;
              for (int a = 0; a < kNMax; ++a) acc += p1d(a, i) * prolonged_stiff(slot, a + kNMax * (b + kNMax * c));
              tmp_i[i + kNMin * (b + kNMax * c)] = acc;
            }

        float tmp_ij[kNMin * kNMin * kNMax];
        for (int c = 0; c < kNMax; ++c)
          for (int j = 0; j < kNMin; ++j)
            for (int i = 0; i < kNMin; ++i) {
              float acc = 0.0f;
              for (int b = 0; b < kNMax; ++b) acc += p1d(b, j) * tmp_i[i + kNMin * (b + kNMax * c)];
              tmp_ij[i + kNMin * (j + kNMin * c)] = acc;
            }

        // No atomic: one thread per pMin element, and the sub-solver kernels are fenced, so this
        // is the only writer of this element's row at this point.
        for (int k = 0; k < kNMin; ++k)
          for (int j = 0; j < kNMin; ++j)
            for (int i = 0; i < kNMin; ++i) {
              float acc = 0.0f;
              for (int c = 0; c < kNMax; ++c) acc += p1d(c, k) * tmp_ij[i + kNMin * (j + kNMin * c)];
              stiff_pMin(e, i + kNMin * (j + kNMin * k)) += acc;
            }
      });
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::ApplyCoupling(const DataType& data) {
  auto mesh_local = m_mesh_;
  auto face_connectivity_local = m_face_connectivity_;
  auto const pField_pMax = data.m_wavefield.m_pMaxAcoustic.getCurrentField(0);
  // pMin is read and written at ORDER_MAX here: ProlongPMinField() raised it beforehand and
  // RestrictPMinStiff() maps the result back. Both sides then run the same code path.
  arrayReal pField_pMin_up = m_pMin_prolonged_field_;
  arrayReal stiff_pMin_up = m_pMin_prolonged_stiff_;
  auto elem_to_slot = m_pMin_elem_to_slot_;

  auto iface_list = m_interface_face_indices_;
  int const n_iface = num_interface_faces_;
  auto element_type = m_element_type_;
  arrayReal stiff_pMax = m_pMax_solver_.m_stiff_local_;
  auto const face_to_elem_dof = pMaxSolver::kFaceToElemDof;
  auto const face_to_elem_dof_depth = pMaxSolver::kFaceToElemDofAtDepth;
  real_t const penalty_local = m_penalty_factor_;

  Kokkos::parallel_for(
      "ApplyCoupling", n_iface, KOKKOS_LAMBDA(const int _loop_idx) {
        int const f = iface_list(_loop_idx);

        int const owner_e = face_connectivity_local.elemOwner(f);
        int const neighbor_e = face_connectivity_local.elemNeighbor(f);
        int const fid_o = face_connectivity_local.localFaceOwner(f);
        int const fid_n = face_connectivity_local.localFaceNeighbor(f);

        bool const owner_is_pMin = (element_type(owner_e) == kElementTypePMin);
        int const pMin_e = (owner_is_pMin) ? owner_e : neighbor_e;
        int const pMax_e = (owner_is_pMin) ? neighbor_e : owner_e;
        int const fid_pMin = (owner_is_pMin) ? fid_o : fid_n;
        int const fid_pMax = (owner_is_pMin) ? fid_n : fid_o;

        auto pMin_to_pMax = [&](int i) {
          return owner_is_pMin ? face_connectivity_local.getNeighborFaceDof(f, i)
                               : face_connectivity_local.getOwnerFaceDof(f, i);
        };
        auto pMax_to_pMin = [&](int i) {
          return owner_is_pMin ? face_connectivity_local.getOwnerFaceDof(f, i)
                               : face_connectivity_local.getNeighborFaceDof(f, i);
        };

        float faceCoords[4][3];
        for (int j = 0; j < 4; ++j) {
          int const gni =
              face_connectivity_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE_MAX::meshIndexToLinearIndex2D(j));
          for (int d = 0; d < 3; ++d) faceCoords[j][d] = mesh_local.nodeCoord(gni, d);
        }

        float pMin_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(pMin_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        pMin_coords[iv + 2 * jv + 4 * kv]);
        }
        float pMax_coords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(pMax_e);
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv),
                                        pMax_coords[iv + 2 * jv + 4 * kv]);
        }

        real_t const inv_rho_min = 1.0f / mesh_local.getModelRhoOnElement(pMin_e);
        real_t const inv_rho_max = 1.0f / mesh_local.getModelRhoOnElement(pMax_e);

        float normal_pMin[3];
        mesh_local.faceNormal(pMin_e, static_cast<model::CubicFace>(fid_pMin), normal_pMin);

        real_t const gamma_min = computeSIPGPenalty<ORDER_MIN>(faceCoords, pMin_coords, penalty_local);
        real_t const gamma_max = computeSIPGPenalty<ORDER_MAX>(faceCoords, pMax_coords, penalty_local);
        // Symmetric penalty: the same gamma on both sides of the interface, otherwise the SIPG
        // bilinear form loses symmetry across the hp-nonconforming face (spurious reflection).
        real_t const gamma_iface = (gamma_min > gamma_max) ? gamma_min : gamma_max;

        int const slot = elem_to_slot(pMin_e);

        // Face-sized accumulators, one per side, each in its own side's face numbering.
        // Element-sized rows would put 2*(ORDER_MAX+1)^3 floats per thread in local memory and
        // route every quadrature update through L1.
        float stiff_min[pMaxSolver::knumNodesPerFace] = {0};
        float stiff_max[pMaxSolver::knumNodesPerFace] = {0};

        float const neg_normal_pMin[3] = {-normal_pMin[0], -normal_pMin[1], -normal_pMin[2]};
        real_t const half_min = 0.5f * inv_rho_min;
        real_t const half_max = 0.5f * inv_rho_max;

        // One quadrature point at a time, both sides fused, at ORDER_MAX resolution. The contracted
        // callbacks fold sum_k C_ijk * n_k, so each contribution fires once instead of once per
        // physical direction. They always fire with j == q, so everything derived from j is
        // hoisted here.
        for (int q = 0; q < pMaxSolver::knumNodesPerFace; ++q) {
          // Face-normal accumulators, one per side. They carry the SIPG consistency channel
          // sum_k C_ijk n_k restricted to the depth direction. On an axis-aligned face this is
          // the only non-zero channel (invJ is diagonal, so the two tangential factors vanish);
          // dropping it would leave the side with nothing but its penalty.
          float norm_min[ORDER_MAX + 1] = {0};
          float norm_max[ORDER_MAX + 1] = {0};

          int const q_to_max = pMin_to_pMax(q);
          int const q_to_min = pMax_to_pMin(q);
          int const ej_min = face_to_elem_dof[fid_pMin][q];
          int const ej_min_perm = face_to_elem_dof[fid_pMax][q_to_max];
          int const ej_max = face_to_elem_dof[fid_pMax][q];
          int const ej_max_perm = face_to_elem_dof[fid_pMin][q_to_min];

          real_t const dp_min = half_min * (pField_pMax(pMax_e, ej_min_perm) - pField_pMin_up(slot, ej_min));
          real_t const dp_max = half_max * (pField_pMin_up(slot, ej_max_perm) - pField_pMax(pMax_e, ej_max));

          // The two other contributions of each side land on fixed slots with equal magnitude and
          // opposite sign, so one register carries both and is flushed once below.
          float acc_min = 0.0f;
          float acc_max = 0.0f;

          // --- pMin side (outward normal = normal_pMin[]) ---
          INTEGRAL_TYPE_MAX::computeInterfaceFluxTermAt(
              q, faceCoords, pMin_coords, fid_pMin, normal_pMin,
              [&](const int i, const int, const real_t val) {
                stiff_min[i] += val * dp_min;
                acc_min -= half_min * val * pField_pMin_up(slot, face_to_elem_dof[fid_pMin][i]);
              },
              [&](const int m, const int, const real_t val) {
                norm_min[m] += val * dp_min;
                acc_min -= half_min * val * pField_pMin_up(slot, face_to_elem_dof_depth[fid_pMin][q][m]);
              });

          // --- pMax side (outward normal = -normal_pMin[]) ---
          INTEGRAL_TYPE_MAX::computeInterfaceFluxTermAt(
              q, faceCoords, pMax_coords, fid_pMax, neg_normal_pMin,
              [&](const int i, const int, const real_t val) {
                stiff_max[i] += val * dp_max;
                acc_max -= half_max * val * pField_pMax(pMax_e, face_to_elem_dof[fid_pMax][i]);
              },
              [&](const int m, const int, const real_t val) {
                norm_max[m] += val * dp_max;
                acc_max -= half_max * val * pField_pMax(pMax_e, face_to_elem_dof_depth[fid_pMax][q][m]);
              });

          // Negation is exact in IEEE-754, so mirroring each register onto the opposite side gives
          // the value the callbacks would have accumulated there.
          stiff_min[q] += acc_min;
          stiff_max[q_to_max] -= acc_min;
          stiff_max[q] += acc_max;
          stiff_min[q_to_min] -= acc_max;

          // Off-face dofs, so they bypass the face-sized flush. Atomic because several faces of the
          // same element write them.
          for (int m = 0; m <= ORDER_MAX; ++m) {
            ATOMICADD(stiff_pMin_up(slot, face_to_elem_dof_depth[fid_pMin][q][m]), norm_min[m]);
            ATOMICADD(stiff_pMax(pMax_e, face_to_elem_dof_depth[fid_pMax][q][m]), norm_max[m]);
          }
        }

        // SIPG penalty and atomic write-back, fused: both sides use the same damping term at face
        // dof i and the same symmetric gamma, so the weight is computed once instead of once per
        // side.
        for (int i = 0; i < pMaxSolver::knumNodesPerFace; ++i) {
          real_t const pen_i = gamma_iface * INTEGRAL_TYPE_MAX::computeDampingTerm(i, faceCoords);

          int const ei_min = face_to_elem_dof[fid_pMin][i];
          int const ei_min_perm = face_to_elem_dof[fid_pMax][pMin_to_pMax(i)];
          stiff_min[i] += pen_i * (pField_pMin_up(slot, ei_min) - pField_pMax(pMax_e, ei_min_perm));

          int const ei_max = face_to_elem_dof[fid_pMax][i];
          int const ei_max_perm = face_to_elem_dof[fid_pMin][pMax_to_pMin(i)];
          stiff_max[i] += pen_i * (pField_pMax(pMax_e, ei_max) - pField_pMin_up(slot, ei_max_perm));

          ATOMICADD(stiff_pMin_up(slot, ei_min), stiff_min[i]);
          ATOMICADD(stiff_pMax(pMax_e, ei_max), stiff_max[i]);
        }
      });
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::computeOneStep(const float& dt, const int& timeSample, DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);

  if (myData.isDistributed) {
    throw std::runtime_error(
        "computeOneStep called in distributed mode. Use computeForces() -> "
        "synchronize() -> updateSolutionForward().");
  }

  // Sub-solver data views are constructed once and reused throughout the step.
  DGsolverDataAcoustic pMin_data(myData.m_wavefield.m_pMinAcoustic, myData.m_rhs.m_rhs_pMinAcoustic);
  DGsolverDataAcoustic pMax_data(myData.m_wavefield.m_pMaxAcoustic, myData.m_rhs.m_rhs_pMaxAcoustic);

  // pMin DG: volume terms and pMin-pMin interior fluxes (interface faces are excluded from the
  // face list).
  m_pMin_solver_.m_list_mode_ = true;
  m_pMin_solver_.m_elem_list_ = pMin_elem_list_;
  m_pMin_solver_.m_n_elem_list_ = num_pMin_elements_;
  m_pMin_solver_.m_face_list_ = m_pMin_interior_face_list_;
  m_pMin_solver_.m_n_face_list_ = m_n_pMin_interior_faces_;

  m_pMin_solver_.applyRHSTerm(timeSample, dt, pMin_data);
  FENCE
  m_pMin_solver_.computeVolumeAndBoundary(num_pMin_elements_, pMin_data.getCurrentField(0));
  FENCE
  m_pMin_solver_.computeBoundaryDampingAndInterfaceFlux(m_n_pMin_interior_faces_, pMin_data.getCurrentField(0));
  FENCE

  // pMax DG: volume terms and pMax-pMax interior fluxes (interface faces are excluded from the
  // face list).
  m_pMax_solver_.m_list_mode_ = true;
  m_pMax_solver_.m_elem_list_ = pMax_elem_list_;
  m_pMax_solver_.m_n_elem_list_ = num_pMax_elements_;
  m_pMax_solver_.m_face_list_ = m_pMax_interior_face_list_;
  m_pMax_solver_.m_n_face_list_ = m_n_pMax_interior_faces_;

  m_pMax_solver_.applyRHSTerm(timeSample, dt, pMax_data);
  FENCE
  m_pMax_solver_.computeVolumeAndBoundary(num_pMax_elements_, pMax_data.getCurrentField(0));
  FENCE
  m_pMax_solver_.computeBoundaryDampingAndInterfaceFlux(m_n_pMax_interior_faces_, pMax_data.getCurrentField(0));
  FENCE

  // Symmetric SIPG interface coupling: both sides read p^n (no temporal lag). Both
  // m_stiff_local_ arrays are complete at this point (FENCE above). The pMin side is raised to
  // ORDER_MAX before the coupling and restricted after, so the interface kernel sees two
  // same-order elements and evaluates the full SIPG form on both.
  ProlongPMinField(myData);
  FENCE

  ApplyCoupling(myData);
  FENCE

  RestrictPMinStiff();
  FENCE

  // Time update of both sub-solvers.
  m_pMin_solver_.applyVerlet(num_pMin_elements_, dt, pMin_data.getCurrentField(0), pMin_data.getPreviousField(0));
  m_pMin_solver_.m_list_mode_ = false;
  FENCE

  m_pMax_solver_.applyVerlet(num_pMax_elements_, dt, pMax_data.getCurrentField(0), pMax_data.getPreviousField(0));
  m_pMax_solver_.m_list_mode_ = false;
  FENCE
}

template <int ORDER_MIN, int ORDER_MAX, template <int, int> class INTEGRAL_SELECTOR, int IMPL_TAG, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES, utils::enums::physicType PHYSICS>
void DGPAdaptiveSolver<ORDER_MIN, ORDER_MAX, INTEGRAL_SELECTOR, IMPL_TAG, MESH_TYPE, IS_MODEL_ON_NODES,
                       PHYSICS>::outputSolutionValues(const int& t, int& e, const arrayReal& field,
                                                      const char* fieldName) {
  cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
       << " after computeOneStep = " << field(e, 0) << endl;
}

}  // namespace fe
}  // namespace solver

#endif  // FUNTIDES_SOLVER_FE_DG_PADAPTIVE_IMPL_COMMON_INCLUDE_DG_PADAPTIVE_SOLVER_IMPL_H_
