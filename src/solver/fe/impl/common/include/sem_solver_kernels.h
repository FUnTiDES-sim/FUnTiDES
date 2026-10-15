#ifndef FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_
#define FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_

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
namespace fe {
namespace detail {
/// @brief True when T declares a nested type `TeamGemm`, i.e. the integral back-end provides the team GEMM kernels.
template <typename, typename = void>
struct has_team_gemm : std::false_type {};
template <typename T>
struct has_team_gemm<T, std::void_t<typename T::TeamGemm>> : std::true_type {};

/// @brief True when T declares a nested type `ZDeformedSumFact`, i.e. the integral back-end provides
/// the z-deformed sum-factorization kernels.
template <typename, typename = void>
struct has_z_deformed_sumfact : std::false_type {};
template <typename T>
struct has_z_deformed_sumfact<T, std::void_t<typename T::ZDeformedSumFact>> : std::true_type {};
}  // namespace detail

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_Gemm(
    const DataType& data) {
  if constexpr (detail::has_team_gemm<INTEGRAL_TYPE>::value) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    constexpr int kStride = INTEGRAL_TYPE::numNodes * 6;

    auto mesh_local = m_mesh;
    int const dim = mesh_local.getOrder() + 1;
    int const pointsPerElem = dim * dim * dim;
    int const nElemsFull = mesh_local.getNumberOfElements();

    // W = w * alpha * B is precomputed once for every element, kStride values per element.
    // It runs on the first call (the benchmark warmup), so it stays out of the timed loop.
    if (!gemmMetricsReady_) {
      gemmMetrics_ = allocateVector<vectorReal>(static_cast<size_t>(nElemsFull) * kStride, "gemmMetrics");
      auto W_global = gemmMetrics_;
      auto mesh_pc = m_mesh;
      Kokkos::parallel_for(
          "Gemm Precompute Metrics", Kokkos::RangePolicy<ExecSpace>(0, nElemsFull), KOKKOS_LAMBDA(const int e) {
            float cornerCoords[8][3];
            auto const eIdx = mesh_pc.elementIndex(e);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_pc.vertexCoords(mesh_pc.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);

            real_t inv_density = 0.0f;
            if constexpr (!IS_MODEL_ON_NODES) {
              inv_density = 1.0f / mesh_pc.getModelRhoOnElement(e);
            }

            INTEGRAL_TYPE::computeElementMetrics(
                cornerCoords,
                [&](const int qa, const int qb, const int qc) -> real_t {
                  if constexpr (IS_MODEL_ON_NODES) {
                    int const g = mesh_pc.globalNodeIndex(e, qa, qb, qc);
                    return 1.0f / mesh_pc.getModelRhoOnNodes(g);
                  } else {
                    return inv_density;
                  }
                },
                W_global.data() + static_cast<size_t>(e) * kStride);
          });
      Kokkos::fence();
      gemmMetricsReady_ = true;
    }

    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : nElemsFull;

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) local_workVectorsGlobal[f] = workVectorsGlobal_[f];

    auto W_global = gemmMetrics_;

    // 1D derivative operator, built on host, captured by value (no device View).
    Kokkos::Array<real_t, INTEGRAL_TYPE::num1dNodes * INTEGRAL_TYPE::num1dNodes> D_arr;
    INTEGRAL_TYPE::fillDerivativeMatrix(D_arr.data());

    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView2D = Kokkos::View<float**, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    TeamPolicyType policy(n_iter, Kokkos::AUTO);
    // Scratch holds the gather/scatter buffers and the 12 [N][N^2] GEMM tensors.
    // W lives in global memory, so it needs no scratch: less shared memory than
    // the streaming variant, which also relaxes the occupancy limit slightly.
    size_t bytes_fields = ScratchView2D::shmem_size(kNumFields, pointsPerElem) * 2;
    size_t bytes_gemm = ScratchView2D::shmem_size(INTEGRAL_TYPE::num1dNodes, INTEGRAL_TYPE::numNodesPerFace) * 12;
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_gemm));

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic Gemm", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          ScratchView2D localFields(team.team_scratch(0), kNumFields, pointsPerElem);
          ScratchView2D localWork(team.team_scratch(0), kNumFields, pointsPerElem);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pointsPerElem), [&](const int localIdx) {
            int i = localIdx % dim;
            int j = (localIdx / dim) % dim;
            int k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            for (int f = 0; f < kNumFields; ++f) {
              localFields(f, localIdx) = data.getCurrentField(f)(globalIdx);
              localWork(f, localIdx) = 0.0f;
            }
          });
          team.team_barrier();

          real_t const* W_ptr = W_global.data() + static_cast<size_t>(elementNumber) * kStride;
          INTEGRAL_TYPE::computeStiffnessOperatorTeamVector(team, &localFields(0, 0), &localWork(0, 0), W_ptr,
                                                            D_arr.data());

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, pointsPerElem), [&](const int localIdx) {
            int i = localIdx % dim;
            int j = (localIdx / dim) % dim;
            int k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork(f, localIdx));
            }
          });
        });
  } else {
    throw std::runtime_error(
        "computeElementContributions_Acoustic_Gemm: INTEGRAL_TYPE has no GEMM path "
        "(use Qk_Hexahedron_Tensorial).");
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic(
    const DataType& data) {
  if constexpr (detail::has_team_gemm<INTEGRAL_TYPE>::value) {
    computeElementContributions_Acoustic_Gemm(data);
  } else {
    if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
      prepareZDeformedGeometry();
      if (zDeformedEnabled_) {
        computeElementContributions_Acoustic_TeamZ(data);
        return;
      }
    }
    computeElementContributions_Acoustic_Flat(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::prepareZDeformedGeometry() {
  if (zDeformedReady_) return;
  zDeformedReady_ = true;
  zDeformedEnabled_ = false;
  zDeformedGeom_ = vectorReal();
  zDeformedNodes_ = vectorInt();
  zDeformedNodesList_ = vectorInt();
  zDeformedInvRhoNodes_ = vectorReal();
  zDeformedBasisTab_ = vectorReal();

  bool const is_acoustic = (PHYSICS == utils::enums::physicType::kAcoustic);

  auto mesh_pc = m_mesh;
  int const nElems = mesh_pc.getNumberOfElements();
  vectorReal geom = allocateVector<vectorReal>(nElems * kZGeomStride, "zDeformedGeom");
  float* geomPtr = geom.data();

  int nonConforming = 0;
  Kokkos::parallel_reduce(
      "ZDeformed Geometry", Kokkos::RangePolicy<>(0, nElems),
      KOKKOS_LAMBDA(const int e, int& bad) {
        float X[8][3];
        auto const eIdx = mesh_pc.elementIndex(e);
        int I = 0;
        for (int kv = 0; kv < 2; ++kv)
          for (int jv = 0; jv < 2; ++jv)
            for (int iv = 0; iv < 2; ++iv) mesh_pc.vertexCoords(mesh_pc.globalVertexIndex(eIdx, iv, jv, kv), X[I++]);

        // Vertex k has indices (k % 2, (k % 4) / 2, k / 4): x may follow the first only, y the second only.
        float const hx = X[1][0] - X[0][0];
        float const hy = X[2][1] - X[0][1];
        float const tolX = 1e-4f * Kokkos::fabs(hx);
        float const tolY = 1e-4f * Kokkos::fabs(hy);
        for (int k = 0; k < 8; ++k) {
          if (Kokkos::fabs(X[k][0] - X[k % 2][0]) > tolX || Kokkos::fabs(X[k][1] - X[2 * ((k % 4) / 2)][1]) > tolY) {
            ++bad;
            break;
          }
        }

        geomPtr[e] = 0.5f * hx;
        geomPtr[nElems + e] = 0.5f * hy;
        for (int k = 0; k < 8; ++k) geomPtr[(2 + k) * nElems + e] = X[k][2];

        // With the model on nodes, 1/rho comes from zDeformedInvRhoNodeTable().
        if (is_acoustic && !IS_MODEL_ON_NODES)
          geomPtr[(kZGeomStride - 1) * nElems + e] = 1.0f / mesh_pc.getModelRhoOnElement(e);
        else
          geomPtr[(kZGeomStride - 1) * nElems + e] = 0.0f;
      },
      nonConforming);

  if (nonConforming > 0) return;

  zDeformedGeom_ = geom;
  zDeformedEnabled_ = true;
  if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    vectorReal basisTab = allocateVector<vectorReal>(INTEGRAL_TYPE::kBasisTableSize, "zDeformedBasisTab");
    auto h_basisTab = Kokkos::create_mirror_view(basisTab);
    for (int idx = 0; idx < INTEGRAL_TYPE::kBasisTableSize; ++idx)
      h_basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx);
    Kokkos::deep_copy(basisTab, h_basisTab);
    zDeformedBasisTab_ = basisTab;
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
int const* SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::zDeformedNodeTable(
    bool element_major) {
  bool const list_on = m_list_mode_;
  int const n_visited = list_on ? m_n_elem_list_ : m_mesh.getNumberOfElements();
  int const* list_key = list_on ? m_elem_list_.data() : nullptr;
  if (zDeformedNodes_.extent(0) == static_cast<size_t>(n_visited) * kPointsPerElement &&
      zDeformedNodesList_.data() == list_key && zDeformedNodesElementMajor_ == element_major)
    return zDeformedNodes_.data();

  auto mesh_pc = m_mesh;
  auto list_local = m_elem_list_;
  // Filled on the device, so that the kernels never read the mesh connectivity, whose pages the
  // host may pull back under unified memory.
  vectorInt nodes = allocateVector<vectorInt>(n_visited * kPointsPerElement, "zDeformedNodes");
  int* nodesPtr = nodes.data();
  Kokkos::parallel_for(
      "ZDeformed Node Table", Kokkos::RangePolicy<>(0, n_visited), KOKKOS_LAMBDA(const int i) {
        constexpr int dim = ORDER + 1;
        int const e = list_on ? list_local[i] : i;
        // The team kernels read the indices element by element, the flat kernel node by node.
        for (int q = 0; q < kPointsPerElement; ++q)
          nodesPtr[element_major ? i * kPointsPerElement + q : q * n_visited + i] =
              mesh_pc.globalNodeIndex(e, q % dim, (q / dim) % dim, q / (dim * dim));
      });

  zDeformedNodes_ = nodes;
  zDeformedNodesList_ = list_on ? m_elem_list_ : vectorInt();
  zDeformedNodesElementMajor_ = element_major;
  return zDeformedNodes_.data();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
float const* SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::zDeformedInvRhoNodeTable() {
  if constexpr (IS_MODEL_ON_NODES) {
    int const n_nodes = m_mesh.getNumberOfNodes();
    if (zDeformedInvRhoNodes_.extent(0) == static_cast<size_t>(n_nodes)) return zDeformedInvRhoNodes_.data();
    auto mesh_pc = m_mesh;
    vectorReal table = allocateVector<vectorReal>(n_nodes, "zDeformedInvRhoNodes");
    float* tablePtr = table.data();
    Kokkos::parallel_for(
        "ZDeformed InvRho Node Table", Kokkos::RangePolicy<>(0, n_nodes),
        KOKKOS_LAMBDA(const int g) { tablePtr[g] = 1.0f / mesh_pc.getModelRhoOnNodes(g); });
    zDeformedInvRhoNodes_ = table;
  }
  return zDeformedInvRhoNodes_.data();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_FlatZ(
    const DataType& data) {
  if constexpr (detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    auto mesh_local = m_mesh;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = mesh_local.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto geom = zDeformedGeom_;
    float const* invRhoNodes = zDeformedInvRhoNodeTable();
    int const* elemNodes = zDeformedNodeTable();
    auto force = workVectorsGlobal_[0];

    // Input and output arrays of one element live in registers: 3 blocks per SM leave room for both.
    using Policy = Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, 3>>;

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic FlatZ", detail::lightWeight(Policy(0, n_iter)),
        KOKKOS_LAMBDA(const int _loop_idx) {
          int const e = list_on ? list_local[_loop_idx] : _loop_idx;
          constexpr int dim = ORDER + 1;

          float u[kPointsPerElement];
          for (int k = 0; k < dim; ++k)
            for (int j = 0; j < dim; ++j)
              for (int i = 0; i < dim; ++i)
                u[i + j * dim + k * dim * dim] =
                    data.getCurrentField(0)(elemNodes[(i + j * dim + k * dim * dim) * n_iter + _loop_idx]);

          real_t const J00 = geom(e);
          real_t const J11 = geom(nElems + e);
          real_t Z[8];
          for (int k = 0; k < 8; ++k) Z[k] = geom((2 + k) * nElems + e);
          real_t const invRhoElem = geom((kZGeomStride - 1) * nElems + e);

          auto get_alpha = [&](const int q) -> real_t {
            if (IS_MODEL_ON_NODES) return invRhoNodes[elemNodes[q * n_iter + _loop_idx]];
            return invRhoElem;
          };

          float v[kPointsPerElement];
          INTEGRAL_TYPE::computeStiffnessTermSumFactZDeformedScatter(J00, J11, Z, u, v, get_alpha);

          for (int k = 0; k < dim; ++k)
            for (int j = 0; j < dim; ++j)
              for (int i = 0; i < dim; ++i)
                ATOMICADD(force[elemNodes[(i + j * dim + k * dim * dim) * n_iter + _loop_idx]],
                          v[i + j * dim + k * dim * dim]);
        });
  } else {
    throw std::runtime_error(
        "computeElementContributions_Acoustic_FlatZ: INTEGRAL_TYPE has no z-deformed sum-factorization path.");
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_TeamZ(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kAcoustic ||
                !detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    throw std::runtime_error("computeElementContributions_Acoustic_TeamZ: needs an acoustic z-deformed setup.");
  } else {
    constexpr int n = ORDER + 1;
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    constexpr int kMinTeamsPerSM = std::clamp(1024 / kPreferredTeamSize, 1, 16);

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<kPreferredTeamSize, kMinTeamsPerSM>>;
    using TeamMember = typename TeamPolicyType::member_type;

    // One thread per node and one per geometry entry are required. The last geometry entry, 1/rho of
    // the element, is only read when the model lives on elements.
    constexpr int kGeomLoaded = IS_MODEL_ON_NODES ? kZGeomStride - 1 : kZGeomStride;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    if (team_size < std::max(kPointsPerElement, kGeomLoaded)) {
      computeElementContributions_Acoustic_FlatZ(data);
      return;
    }

    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = m_mesh.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto zgeom = zDeformedGeom_;
    int const* elemNodes = zDeformedNodeTable(true);
    float const* invRhoNodes = zDeformedInvRhoNodeTable();
    float const* basisTabGlobal = zDeformedBasisTab_.data();
    auto force = workVectorsGlobal_[0];

    // Same choice as the TTI team kernel: up to 8 elements per team, at least four waves of teams.
    constexpr int kMaxElemsPerTeam = 8;
    constexpr int kMinWaves = 4;
    int const resident_teams = std::max<int>(1, ExecSpace::concurrency() / kPreferredTeamSize);
    int const elems_per_team = std::clamp(n_iter / (kMinWaves * resident_teams), 1, kMaxElemsPerTeam);
    TeamPolicyType policy((n_iter + elems_per_team - 1) / elems_per_team, team_size);

    // Shared scratch, every block on a 16-byte boundary: rows of D padded to a multiple of 4, element
    // geometry, nodal values, then the three fluxes.
    constexpr int kRowChunks = (n + 3) / 4;
    constexpr int kRowStride = 4 * kRowChunks;
    constexpr int kGeomPadded = ((kZGeomStride + 3) / 4) * 4;
    static_assert(kGeomPadded == 12, "the geometry is read as three 16-byte blocks");
    constexpr int kGeomOffset = n * kRowStride;
    constexpr int kValueOffset = kGeomOffset + kGeomPadded;
    constexpr int kFluxOffset = kValueOffset + kPointsPerElement;
    constexpr int kScratchFloats = kFluxOffset + 3 * kPointsPerElement;
    constexpr size_t kScratchBytes = kScratchFloats * sizeof(float);
    policy.set_scratch_size(0, Kokkos::PerTeam(kScratchBytes + alignof(detail::Float4)));

    Kokkos::parallel_for(
        "Solver Element Contribution Acoustic TeamZ", detail::lightWeight(policy),
        KOKKOS_LAMBDA(const TeamMember& team) {
          float* scratch =
              static_cast<float*>(team.team_scratch(0).get_shmem_aligned(kScratchBytes, alignof(detail::Float4)));
          // Captured here rather than first inside an if constexpr branch, which nvcc rejects.
          float const* inv_rho_nodes = invRhoNodes;
          auto const geom = zgeom;
          float* rows = scratch;
          float* geom_sh = scratch + kGeomOffset;
          float* u = scratch + kValueOffset;
          float* G = scratch + kFluxOffset;
          detail::Float4 const* rows4 = reinterpret_cast<detail::Float4 const*>(rows);
          detail::Float4 const* geom4 = reinterpret_cast<detail::Float4 const*>(geom_sh);

          // Thread q owns point q, which is also node q.
          int const q = team.team_rank();
          bool const owns = q < kPointsPerElement;
          int const qa = q % n;
          int const qb = (q / n) % n;
          int const qc = q / (n * n);

          // Filled before the first barrier, which covers it.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n * kRowStride), [&](const int idx) {
            int const r = idx / kRowStride;
            int const p = idx % kRowStride;
            rows[idx] = p < n ? basisTabGlobal[r * n + p] : 0.0f;
          });

          // The same for every element of the team: basis derivatives along the three lines through q,
          // the quadrature weight and the interpolation coordinates of q.
          float d_xi[n];
          float d_eta[n];
          float d_zeta[n];
          float weight = 0.0f;
          float alpha[3] = {0.0f, 0.0f, 0.0f};
          if (owns) {
            for (int i = 0; i < n; ++i) {
              d_xi[i] = basisTabGlobal[i * n + qa];
              d_eta[i] = basisTabGlobal[i * n + qb];
              d_zeta[i] = basisTabGlobal[i * n + qc];
            }
            weight = basisTabGlobal[n * n + qa] * basisTabGlobal[n * n + qb] * basisTabGlobal[n * n + qc];
            alpha[0] = basisTabGlobal[n * n + n + qa];
            alpha[1] = basisTabGlobal[n * n + n + qb];
            alpha[2] = basisTabGlobal[n * n + n + qc];
          }

          int const first_rank = team.league_rank() * elems_per_team;
          int const n_here = Kokkos::min(elems_per_team, n_iter - first_rank);

          // Two-stage pipeline held in registers: node index and element number two elements ahead, values
          // and one geometry entry one element ahead, so that a load never waits on its index.
          int ahead_node = 0;
          int ahead_e = 0;
          int next_node = 0;
          float next_u = 0.0f;
          float next_inv_rho = 0.0f;
          float next_geom = 0.0f;
          auto const load_index = [&](int rank) {
            if (owns) ahead_node = elemNodes[rank * kPointsPerElement + q];
            if (q < kGeomLoaded) ahead_e = list_on ? list_local[rank] : rank;
          };
          auto const load_values = [&]() {
            if (owns) {
              next_node = ahead_node;
              next_u = data.getCurrentField(0)(next_node);
              if constexpr (IS_MODEL_ON_NODES) next_inv_rho = inv_rho_nodes[next_node];
            }
            if (q < kGeomLoaded) next_geom = geom(q * nElems + ahead_e);
          };
          load_index(first_rank);
          load_values();
          if (n_here > 1) load_index(first_rank + 1);

          for (int k = 0; k < n_here; ++k) {
            int const node = next_node;
            float inv_rho = next_inv_rho;

            if (owns) u[q] = next_u;
            if (q < kGeomLoaded) geom_sh[q] = next_geom;
            team.team_barrier();

            if (k + 1 < n_here) load_values();
            if (k + 2 < n_here) load_index(first_rank + k + 2);

            // Entries J00, J11, the eight vertex z and 1/rho of the element, in that order.
            detail::Float4 const g0 = geom4[0];
            detail::Float4 const g1 = geom4[1];
            detail::Float4 const g2 = geom4[2];
            float const J00 = g0.v[0];
            float const J11 = g0.v[1];
            float const Z[8] = {g0.v[2], g0.v[3], g1.v[0], g1.v[1], g1.v[2], g1.v[3], g2.v[0], g2.v[1]};
            if constexpr (!IS_MODEL_ON_NODES) inv_rho = g2.v[2];

            if (owns) {
              float du[3] = {0.0f, 0.0f, 0.0f};
              for (int i = 0; i < n; ++i) {
                du[0] += d_xi[i] * u[i + qb * n + qc * n * n];
                du[1] += d_eta[i] * u[qa + i * n + qc * n * n];
                du[2] += d_zeta[i] * u[qa + qb * n + i * n * n];
              }

              float dZ[3][4];
              INTEGRAL_TYPE::zDeformedHalfDifferences(Z, dZ);
              float invJ[3][3];
              float const detJ = INTEGRAL_TYPE::invJacobianZDeformedBilinear(alpha, J00, J11, dZ, invJ);
              float const a = invJ[0][0];
              float const b = invJ[1][1];
              float const c = invJ[2][0];
              float const d = invJ[2][1];
              float const ez = invJ[2][2];
              // B = det(J) J^-1 J^-T, whose xy entry is zero.
              float const scale = weight * inv_rho * detJ;
              float const B0 = a * a;
              float const B1 = b * b;
              float const B2 = c * c + d * d + ez * ez;
              float const B3 = b * d;
              float const B4 = a * c;
              G[q] = scale * (B0 * du[0] + B4 * du[2]);
              G[kPointsPerElement + q] = scale * (B1 * du[1] + B3 * du[2]);
              G[2 * kPointsPerElement + q] = scale * (B4 * du[0] + B3 * du[1] + B2 * du[2]);
            }
            team.team_barrier();

            // Transpose: node (qa, qb, qc) gathers D(qa, p) G_xi(p, qb, qc), and so on.
            if (owns) {
              float v = 0.0f;
              for (int c = 0; c < kRowChunks; ++c) {
                detail::Float4 const ra = rows4[qa * kRowChunks + c];
                detail::Float4 const rb = rows4[qb * kRowChunks + c];
                detail::Float4 const rc = rows4[qc * kRowChunks + c];
                for (int j = 0; j < 4; ++j) {
                  int const p = 4 * c + j;
                  if (p < n) {
                    v += ra.v[j] * G[p + qb * n + qc * n * n];
                    v += rb.v[j] * G[kPointsPerElement + qa + p * n + qc * n * n];
                    v += rc.v[j] * G[2 * kPointsPerElement + qa + qb * n + p * n * n];
                  }
                }
              }
              ATOMICADD(force[node], v);
            }
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Acoustic_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  using Policy = Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>;

  Kokkos::parallel_for(
      "Solver Element Contribution Acoustic", Policy(0, n_iter), KOKKOS_LAMBDA(const int _loop_idx) {
        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        constexpr int dim = ORDER + 1;

        float localFields[kNumFields][kPointsPerElement];
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int k = 0; k < dim; ++k) {
          for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        real_t inv_density = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          inv_density = 1.0f / mesh_local.getModelRhoOnElement(elementNumber);
        }

        INTEGRAL_TYPE::computeStiffnessTermSumFact(
            cornerCoords, localFields[0], localWork[0], [&](const int qa, const int qb, const int qc) -> real_t {
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                return 1.0f / mesh_local.getModelRhoOnNodes(gIndex);
              } else {
                return inv_density;
              }
            });

        for (int k = 0; k < dim; ++k) {
          for (int j = 0; j < dim; ++j) {
            for (int i = 0; i < dim; ++i) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) {
                ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributions(
    const DataType& data) {
  if (!attenuationEnabled_ || nSls_ <= 0) return;

  if constexpr (PHYSICS == utils::enums::physicType::kAcoustic)
    computeAttenuationContributionsAcoustic(data);
  else
    computeAttenuationContributionsElastic(data);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsAcoustic(
    const DataType& data) {
  auto mesh_local = m_mesh;
  Kokkos::Array<vectorReal, kNumFields> local_attenuationWorkVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];

  Kokkos::parallel_for(
      "Solver Attenuation Contributions Acoustic",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWorkA[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
            }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        real_t inv_density_q = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          inv_density_q =
              1.0f / (mesh_local.getModelRhoOnElement(elementNumber) * mesh_local.getModelQpOnElement(elementNumber));
        }

        INTEGRAL_TYPE::computeStiffnessTerm(
            cornerCoords,
            [&](const int qa, const int qb, const int qc) {
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                inv_density_q = 1.0f / (mesh_local.getModelRhoOnNodes(gIndex) * mesh_local.getModelQpOnNodes(gIndex));
              }
            },
            [&](const int i, const int j, const real_t val) {
              localWorkA[0][i] += inv_density_q * val * localFields[0][j];
            });

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f)
                ATOMICADD(local_attenuationWorkVectorsGlobal[f][globalIdx], localWorkA[f][localIdx]);
            }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeAttenuationContributionsElastic(
    const DataType& data) {
  if (anisotropyType_ != model::AnisotropyType::kIso) return;

  auto mesh_local = m_mesh;

  Kokkos::Array<vectorReal, kNumFields> local_attenuationWorkVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) local_attenuationWorkVectorsGlobal[f] = attenuationWorkVectorsGlobal_[f];

  Kokkos::parallel_for(
      "Solver Attenuation Contributions Elastic",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        constexpr int dim = ORDER + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWorkA[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f) localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
            }

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        // Per (p, r) pair of reference directions: the six distinct coefficients of the
        // isotropic stiffness contracted with the inverse Jacobian.
        // a0, a1, a2 multiply the diagonal displacement terms; a3, b0, b1 the xy, xz, yz couplings.
        struct CJPacked {
          float a0, a1, a2, a3;
          float b0, b1;
        };
        CJPacked CJflat[3 * 3];

        INTEGRAL_TYPE::computeStiffNessTermwithJac(
            cornerCoords,
            [&](int qa, int qb, int qc, float const(&J)[3][3]) {
              float vp, vs, rho, qp, qs;
              if constexpr (IS_MODEL_ON_NODES) {
                int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                vp = mesh_local.getModelVpOnNodes(gIndex);
                vs = mesh_local.getModelVsOnNodes(gIndex);
                rho = mesh_local.getModelRhoOnNodes(gIndex);
                qp = mesh_local.getModelQpOnNodes(gIndex);
                qs = mesh_local.getModelQsOnNodes(gIndex);
              } else {
                vp = mesh_local.getModelVpOnElement(elementNumber);
                vs = mesh_local.getModelVsOnElement(elementNumber);
                rho = mesh_local.getModelRhoOnElement(elementNumber);
                qp = mesh_local.getModelQpOnElement(elementNumber);
                qs = mesh_local.getModelQsOnElement(elementNumber);
              }

              float mu = rho * vs * vs;
              float lambda = rho * (vp * vp - 2.0f * vs * vs);
              float lambdap2mua = (lambda + 2.0f * mu) / qp;
              mu = mu / qs;
              lambda = lambdap2mua - 2.0f * mu;
              float lambda_plus_2mu = lambda + 2.0f * mu;

              for (int p = 0; p < 3; ++p) {
                float const Jp0 = J[p][0], Jp1 = J[p][1], Jp2 = J[p][2];
                for (int r = 0; r < 3; ++r) {
                  float const Jr0 = J[r][0], Jr1 = J[r][1], Jr2 = J[r][2];
                  int const idx = p * 3 + r;
                  CJflat[idx].a0 = lambda_plus_2mu * Jp0 * Jr0 + mu * (Jp1 * Jr1 + Jp2 * Jr2);
                  CJflat[idx].a1 = mu * Jp0 * Jr0 + lambda_plus_2mu * Jp1 * Jr1 + mu * Jp2 * Jr2;
                  CJflat[idx].a2 = mu * (Jp0 * Jr0 + Jp1 * Jr1) + lambda_plus_2mu * Jp2 * Jr2;
                  CJflat[idx].a3 = lambda * Jp0 * Jr1 + mu * Jp1 * Jr0;
                  CJflat[idx].b0 = lambda * Jp0 * Jr2 + mu * Jp2 * Jr0;
                  CJflat[idx].b1 = lambda * Jp1 * Jr2 + mu * Jp2 * Jr1;
                }
              }
            },
            [&](int i, int j, float val, const int p, const int r) {
              int const idx = p * 3 + r;
              float const uxj = localFields[0][j];
              float const uyj = localFields[1][j];
              float const uzj = localFields[2][j];
              localWorkA[0][i] += val * (CJflat[idx].a0 * uxj + CJflat[idx].a3 * uyj + CJflat[idx].b0 * uzj);
              localWorkA[1][i] += val * (CJflat[idx].a3 * uxj + CJflat[idx].a1 * uyj + CJflat[idx].b1 * uzj);
              localWorkA[2][i] += val * (CJflat[idx].b0 * uxj + CJflat[idx].b1 * uyj + CJflat[idx].a2 * uzj);
            });

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;
              for (int f = 0; f < kNumFields; ++f)
                ATOMICADD(local_attenuationWorkVectorsGlobal[f][globalIdx], localWorkA[f][localIdx]);
            }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Element Contribution Iso Flat",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;

        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float mu_e = 0.0f, lambda_e = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            mu_e = rho_e * vs_e * vs_e;
            lambda_e = rho_e * (vp_e * vp_e - 2.0f * vs_e * vs_e);
          }

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork,
                                                        [&](int qa, int qb, int qc, float const(&J_inv)[3][3],
                                                            float const(&grad_u_ref)[3][3], float(&flux)[3][3]) {
                                                          float mu, lambda;
                                                          if constexpr (IS_MODEL_ON_NODES) {
                                                            int const gIndex =
                                                                mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                                                            float const vp = mesh_local.getModelVpOnNodes(gIndex);
                                                            float const vs = mesh_local.getModelVsOnNodes(gIndex);
                                                            float const rho = mesh_local.getModelRhoOnNodes(gIndex);
                                                            mu = rho * vs * vs;
                                                            lambda = rho * (vp * vp - 2.0f * vs * vs);
                                                          } else {
                                                            mu = mu_e;
                                                            lambda = lambda_e;
                                                          }

                                                          flux::elasticFluxIso(J_inv, mu, lambda, grad_u_ref, flux);
                                                        });

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Iso_Flat(data);
  } else {
    computeElementContributions_Iso_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Iso_Team(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  // avoid extended __host__ __device__ lambda cannot first-capture variable in
  // constexpr-if context
  (void)mesh_local;
  (void)list_on;
  (void)list_local;
  (void)n_iter;
  (void)local_workVectorsGlobal;

  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;

    // A team owns one element: gather and scatter are spread over its threads,
    // which makes them coalesced, and the per-element buffers move from
    // per-thread local memory to shared memory.
    //
    // Every team range covers kPointsPerElement, so a larger team only adds
    // threads that idle at the barriers. Round up to a whole warp so no partial
    // warp idles where kPointsPerElement is not a multiple of 32, then clamp to
    // what the backend can actually provide: a host backend caps the team size
    // at its thread count, and asking for more aborts at launch. A
    // TeamThreadRange is correct at any team size.
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    // Constant-Jacobian meshes carry the element geometry (9 inverse-Jacobian
    // entries + the determinant) in a small per-team scratch buffer.
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;
    TeamPolicyType policy(n_iter, team_size);
    // One buffer carries the displacements in and the forces out: the stiffness
    // kernel overwrites it after the barrier that ends its read phase. Halving
    // this buffer is what lifts the occupancy off its shared-memory limit.
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Iso Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float mu_e = 0.0f, lambda_e = 0.0f;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            mu_e = rho_e * vs_e * vs_e;
            lambda_e = rho_e * (vp_e * vp_e - 2.0f * vs_e * vs_e);
          }

          auto const iso_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            float mu, lambda;
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float const vp = mesh_local.getModelVpOnNodes(gIndex);
              float const vs = mesh_local.getModelVsOnNodes(gIndex);
              float const rho = mesh_local.getModelRhoOnNodes(gIndex);
              mu = rho * vs * vs;
              lambda = rho * (vp * vp - 2.0f * vs * vs);
            } else {
              mu = mu_e;
              lambda = lambda_e;
            }
            flux::elasticFluxIso(J_inv, mu, lambda, grad_u_ref, flux);
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              // invJacobianTransformation accumulates into J, so it must start at zero.
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), iso_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), iso_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti_Flat(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Element Contribution Vti Flat",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
      KOKKOS_LAMBDA(const int _loop_idx) {
        // avoid extended __host__ __device__ lambda cannot first-capture
        // variable in constexpr-if context
        (void)local_workVectorsGlobal;
        int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

        int const dim = mesh_local.getOrder() + 1;
        float localFields[kNumFields][kPointsPerElement] = {{0}};
        float localWork[kNumFields][kPointsPerElement] = {{0}};

        for (int i = 0; i < dim; ++i) {
          for (int j = 0; j < dim; ++j) {
            for (int k = 0; k < dim; ++k) {
              int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
              int const localIdx = i + j * dim + k * dim * dim;

              for (int f = 0; f < kNumFields; ++f) {
                localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
              }
            }
          }
        }

        if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float c11_e = 0, c12_e = 0, c13_e = 0, c33_e = 0, c44_e = 0, c66_e = 0;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            float const delta_e = mesh_local.getModelDeltaOnElement(elementNumber);
            float const epsilon_e = mesh_local.getModelEpsilonOnElement(elementNumber);
            float const gamma_e = mesh_local.getModelGammaOnElement(elementNumber);
            float const rho_vp2 = rho_e * vp_e * vp_e;
            float const rho_vs2 = rho_e * vs_e * vs_e;
            c33_e = rho_vp2;
            c44_e = rho_vs2;
            c11_e = rho_vp2 * (1.0f + 2.0f * epsilon_e);
            c66_e = rho_vs2 * (1.0f + 2.0f * gamma_e);
            float const vp2_vs2 = vp_e * vp_e - vs_e * vs_e;
            c13_e = rho_e * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp_e * vp_e * delta_e * vp2_vs2) - rho_vs2;
            c12_e = c11_e - 2.0f * c66_e;
          }

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(
              cornerCoords, localFields, localWork,
              [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                  float(&flux)[3][3]) {
                float c11, c12, c13, c33, c44, c66;
                if constexpr (IS_MODEL_ON_NODES) {
                  int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
                  float const vp = mesh_local.getModelVpOnNodes(gIndex);
                  float const vs = mesh_local.getModelVsOnNodes(gIndex);
                  float const rho = mesh_local.getModelRhoOnNodes(gIndex);
                  float const delta = mesh_local.getModelDeltaOnNodes(gIndex);
                  float const epsilon = mesh_local.getModelEpsilonOnNodes(gIndex);
                  float const gamma = mesh_local.getModelGammaOnNodes(gIndex);
                  float const rho_vp2 = rho * vp * vp;
                  float const rho_vs2 = rho * vs * vs;
                  c33 = rho_vp2;
                  c44 = rho_vs2;
                  c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
                  c66 = rho_vs2 * (1.0f + 2.0f * gamma);
                  float const vp2_vs2 = vp * vp - vs * vs;
                  c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;
                  c12 = c11 - 2.0f * c66;
                } else {
                  c11 = c11_e;
                  c12 = c12_e;
                  c13 = c13_e;
                  c33 = c33_e;
                  c44 = c44_e;
                  c66 = c66_e;
                }

                flux::elasticFluxVti(J_inv, c11, c12, c13, c33, c44, c66, grad_u_ref, flux);
              });

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Vti_Flat(data);
  } else {
    computeElementContributions_Vti_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Vti_Team(
    const DataType& data) {
  auto mesh_local = m_mesh;
  bool const list_on = m_list_mode_;
  auto list_local = m_elem_list_;
  int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

  std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_workVectorsGlobal[f] = workVectorsGlobal_[f];
  }

  // avoid extended __host__ __device__ lambda cannot first-capture variable in
  // constexpr-if context
  (void)mesh_local;
  (void)list_on;
  (void)list_local;
  (void)local_workVectorsGlobal;

  if constexpr (PHYSICS == utils::enums::physicType::kElastic) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;
    // Warp-aligned preferred size, clamped to what the backend can provide (a
    // host backend caps it at its thread count; asking for more aborts).
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;

    TeamPolicyType policy(n_iter, team_size);
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Vti Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          float c11_e = 0, c12_e = 0, c13_e = 0, c33_e = 0, c44_e = 0, c66_e = 0;
          if constexpr (!IS_MODEL_ON_NODES) {
            float const vp_e = mesh_local.getModelVpOnElement(elementNumber);
            float const vs_e = mesh_local.getModelVsOnElement(elementNumber);
            float const rho_e = mesh_local.getModelRhoOnElement(elementNumber);
            float const delta_e = mesh_local.getModelDeltaOnElement(elementNumber);
            float const epsilon_e = mesh_local.getModelEpsilonOnElement(elementNumber);
            float const gamma_e = mesh_local.getModelGammaOnElement(elementNumber);
            float const rho_vp2 = rho_e * vp_e * vp_e;
            float const rho_vs2 = rho_e * vs_e * vs_e;
            c33_e = rho_vp2;
            c44_e = rho_vs2;
            c11_e = rho_vp2 * (1.0f + 2.0f * epsilon_e);
            c66_e = rho_vs2 * (1.0f + 2.0f * gamma_e);
            float const vp2_vs2 = vp_e * vp_e - vs_e * vs_e;
            c13_e = rho_e * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp_e * vp_e * delta_e * vp2_vs2) - rho_vs2;
            c12_e = c11_e - 2.0f * c66_e;
          }

          auto const vti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            float c11, c12, c13, c33, c44, c66;
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float const vp = mesh_local.getModelVpOnNodes(gIndex);
              float const vs = mesh_local.getModelVsOnNodes(gIndex);
              float const rho = mesh_local.getModelRhoOnNodes(gIndex);
              float const delta = mesh_local.getModelDeltaOnNodes(gIndex);
              float const epsilon = mesh_local.getModelEpsilonOnNodes(gIndex);
              float const gamma = mesh_local.getModelGammaOnNodes(gIndex);
              float const rho_vp2 = rho * vp * vp;
              float const rho_vs2 = rho * vs * vs;
              c33 = rho_vp2;
              c44 = rho_vs2;
              c11 = rho_vp2 * (1.0f + 2.0f * epsilon);
              c66 = rho_vs2 * (1.0f + 2.0f * gamma);
              float const vp2_vs2 = vp * vp - vs * vs;
              c13 = rho * sqrtf(vp2_vs2 * vp2_vs2 + 2.0f * vp * vp * delta * vp2_vs2) - rho_vs2;
              c12 = c11 - 2.0f * c66;
            } else {
              c11 = c11_e;
              c12 = c12_e;
              c13 = c13_e;
              c33 = c33_e;
              c44 = c44_e;
              c66 = c66_e;
            }
            flux::elasticFluxVti(J_inv, c11, c12, c13, c33, c44, c66, grad_u_ref, flux);
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), vti_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), vti_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

// Stores the compact TTI description of every node (see computeTtiCompact()) into cttiNodes_.
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::precomputeTtiTensorsOnNodes() {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic || !IS_MODEL_ON_NODES) {
    return;
  } else {
    if (cttiNodesReady_) {
      return;
    }

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    int const numNodes = m_mesh.getNumberOfNodes();
    cttiNodes_ = allocateArray2D<CttiView>(numNodes, flux::kTtiCompactSize, "cttiNodes");

    auto ctti = cttiNodes_;
    auto mesh_pc = m_mesh;
    Kokkos::parallel_for(
        "Tti Precompute Compact On Nodes", Kokkos::RangePolicy<ExecSpace>(0, numNodes), KOKKOS_LAMBDA(const int g) {
          float p[flux::kTtiCompactSize];
          computeTtiCompact(mesh_pc.getModelVpOnNodes(g), mesh_pc.getModelVsOnNodes(g), mesh_pc.getModelRhoOnNodes(g),
                            mesh_pc.getModelDeltaOnNodes(g), mesh_pc.getModelEpsilonOnNodes(g),
                            mesh_pc.getModelGammaOnNodes(g), mesh_pc.getModelPhiOnNodes(g),
                            mesh_pc.getModelThetaOnNodes(g), p);
          for (int k = 0; k < flux::kTtiCompactSize; ++k) ctti(g, k) = p[k];
        });
    Kokkos::fence();

    cttiNodesReady_ = true;
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_Flat(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic) {
  } else {
    precomputeTtiTensorsOnNodes();

    auto mesh_local = m_mesh;
    [[maybe_unused]] auto ctti_local = cttiNodes_;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    Kokkos::parallel_for(
        "Solver Element Contribution Tti Flat",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_iter),
        KOKKOS_LAMBDA(const int _loop_idx) {
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          int const dim = mesh_local.getOrder() + 1;
          float localFields[kNumFields][kPointsPerElement] = {{0}};
          float localWork[kNumFields][kPointsPerElement] = {{0}};

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  localFields[f][localIdx] = data.getCurrentField(f)(globalIdx);
                }
              }
            }
          }

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          // Captured here, not at first use: nvcc rejects an extended lambda that first-captures
          // a variable inside an if-constexpr branch.
          [[maybe_unused]] auto const ctti = ctti_local;

          // Model on elements only: one tensor for the whole element.
          [[maybe_unused]] float CTTI[6][6] = {};
          if constexpr (!IS_MODEL_ON_NODES) {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

          auto ttiFlux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                             float(&flux)[3][3]) {
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = mesh_local.globalNodeIndex(elementNumber, qa, qb, qc);
              float p[flux::kTtiCompactSize];
              for (int k = 0; k < flux::kTtiCompactSize; ++k) p[k] = ctti(gIndex, k);
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            } else {
              flux::elasticFluxTti(J_inv, CTTI, grad_u_ref, flux);
            }
          };

          INTEGRAL_TYPE::computeElasticStiffnessSumFact(cornerCoords, localFields, localWork, ttiFlux);

          for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
              for (int k = 0; k < dim; ++k) {
                int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
                int const localIdx = i + j * dim + k * dim * dim;

                for (int f = 0; f < kNumFields; ++f) {
                  ATOMICADD(local_workVectorsGlobal[f][globalIdx], localWork[f][localIdx]);
                }
              }
            }
          }
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti(
    const DataType& data) {
  if constexpr (ORDER <= kMaxOrderForFlatElastic) {
    computeElementContributions_Tti_Flat(data);
  } else {
    computeElementContributions_Tti_Team(data);
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_Team(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic) {
  } else {
    if constexpr (!HasConstantJacobian<MESH_TYPE>::value && IS_MODEL_ON_NODES &&
                  detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
      prepareZDeformedGeometry();
      if (zDeformedEnabled_) {
        precomputeTtiTensorsOnNodes();
        computeElementContributions_Tti_TeamZ(data);
        return;
      }
    }
    precomputeTtiTensorsOnNodes();

    auto mesh_local = m_mesh;
    [[maybe_unused]] auto ctti_local = cttiNodes_;
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const n_iter = list_on ? m_n_elem_list_ : mesh_local.getNumberOfElements();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchView1D = Kokkos::View<float*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    constexpr int dim = ORDER + 1;
    // Warp-aligned preferred size, clamped to what the backend can provide (a
    // host backend caps it at its thread count; asking for more aborts).
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    constexpr bool kConstJac = HasConstantJacobian<MESH_TYPE>::value;

    TeamPolicyType policy(n_iter, team_size);
    size_t const bytes_fields = ScratchView1D::shmem_size(kNumFields * kPointsPerElement);
    size_t const bytes_flux = ScratchView1D::shmem_size(9 * kPointsPerElement);
    size_t const bytes_geom = kConstJac ? ScratchView1D::shmem_size(10) : 0;
    size_t const bytes_tab = ScratchView1D::shmem_size(INTEGRAL_TYPE::kBasisTableSize);
    size_t const bytes_idx = ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_fields + bytes_flux + bytes_geom + bytes_tab + bytes_idx));

    Kokkos::parallel_for(
        "Solver Element Contribution Tti Team", policy, KOKKOS_LAMBDA(const TeamMember& team) {
          int const _loop_idx = team.league_rank();
          int const elementNumber = list_on ? list_local[_loop_idx] : _loop_idx;

          // Displacements in, forces out, same storage.
          ScratchView1D localFields(team.team_scratch(0), kNumFields * kPointsPerElement);
          ScratchView1D fluxScratch(team.team_scratch(0), 9 * kPointsPerElement);
          ScratchView1D basisTab(team.team_scratch(0), INTEGRAL_TYPE::kBasisTableSize);
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Filled before the gather barrier, so it needs no barrier of its own.
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                               [&](const int idx) { basisTab(idx) = INTEGRAL_TYPE::basisTableEntry(idx); });

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const i = localIdx % dim;
            int const j = (localIdx / dim) % dim;
            int const k = localIdx / (dim * dim);
            int const globalIdx = mesh_local.globalNodeIndex(elementNumber, i, j, k);
            nodeIdx(localIdx) = globalIdx;
            for (int f = 0; f < kNumFields; ++f)
              localFields(f * kPointsPerElement + localIdx) = data.getCurrentField(f)(globalIdx);
          });
          team.team_barrier();

          float cornerCoords[8][3];
          {
            auto const eIdx = mesh_local.elementIndex(elementNumber);
            int I = 0;
            for (int kv = 0; kv < 2; ++kv)
              for (int jv = 0; jv < 2; ++jv)
                for (int iv = 0; iv < 2; ++iv)
                  mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
          }

          // Captured here, not at first use: nvcc rejects an extended lambda that first-captures
          // a variable inside an if-constexpr branch.
          [[maybe_unused]] auto const ctti = ctti_local;

          // Model on elements only: one tensor for the whole element.
          [[maybe_unused]] float CTTI[6][6] = {};
          if constexpr (!IS_MODEL_ON_NODES) {
            mesh_local.getCTensorOnElement(elementNumber, CTTI);
          }

          auto const tti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                    float(&flux)[3][3]) {
            if constexpr (IS_MODEL_ON_NODES) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float p[flux::kTtiCompactSize];
              for (int k = 0; k < flux::kTtiCompactSize; ++k) p[k] = ctti(gIndex, k);
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            } else {
              flux::elasticFluxTti(J_inv, CTTI, grad_u_ref, flux);
            }
          };

          if constexpr (HasConstantJacobian<MESH_TYPE>::value) {
            ScratchView1D geom(team.team_scratch(0), 10);
            if (team.team_rank() == 0) {
              float J_inv[3][3] = {{0}};
              float const detJ = INTEGRAL_TYPE::invJacobianTransformation(0, 0, 0, cornerCoords, J_inv);
              for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) geom(a * 3 + b) = J_inv[a][b];
              geom(9) = detJ;
            }
            team.team_barrier();
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, &geom(0), &localFields(0), &localFields(0),
                                                              &fluxScratch(0), tti_flux, &basisTab(0));
          } else {
            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeam(team, cornerCoords, &localFields(0), &localFields(0),
                                                              &fluxScratch(0), tti_flux, &basisTab(0));
          }
          team.team_barrier();

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
            int const globalIdx = nodeIdx(localIdx);
            for (int f = 0; f < kNumFields; ++f) {
              ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields(f * kPointsPerElement + localIdx));
            }
          });
        });
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeElementContributions_Tti_TeamZ(
    const DataType& data) {
  if constexpr (PHYSICS != utils::enums::physicType::kElastic || !IS_MODEL_ON_NODES ||
                !detail::has_z_deformed_sumfact<INTEGRAL_TYPE>::value) {
    throw std::runtime_error("computeElementContributions_Tti_TeamZ: needs an elastic on-nodes z-deformed setup.");
  } else {
    auto ctti_local = cttiNodes_;
    // Two 16-byte loads per node instead of eight scalar ones: each scalar load of a warp
    // touches the same cache lines again.
    static_assert(flux::kTtiCompactSize == 8, "two 16-byte loads per node");
    if (ctti_local.stride(0) != flux::kTtiCompactSize ||
        reinterpret_cast<std::uintptr_t>(ctti_local.data()) % alignof(detail::Float4) != 0)
      throw std::runtime_error("computeElementContributions_Tti_TeamZ: cttiNodes_ not packed for 16-byte loads.");
    detail::Float4 const* ctti_quads = reinterpret_cast<detail::Float4 const*>(ctti_local.data());
    bool const list_on = m_list_mode_;
    auto list_local = m_elem_list_;
    int const nElems = m_mesh.getNumberOfElements();
    int const n_iter = list_on ? m_n_elem_list_ : nElems;
    auto zgeom = zDeformedGeom_;
    int const* elemNodes = zDeformedNodeTable();
    float const* basisTabGlobal = zDeformedBasisTab_.data();

    std::array<std::remove_reference_t<decltype(workVectorsGlobal_[0])>, kNumFields> local_workVectorsGlobal;
    for (int f = 0; f < kNumFields; ++f) {
      local_workVectorsGlobal[f] = workVectorsGlobal_[f];
    }

    constexpr int dim = ORDER + 1;
    constexpr int kPreferredTeamSize = ((kPointsPerElement + 31) / 32) * 32;
    // Caps registers at 64 per thread; at most 16 teams, the lowest resident-block limit per SM.
    constexpr int kMinTeamsPerSM = std::clamp(1024 / kPreferredTeamSize, 1, 16);

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using TeamPolicyType = Kokkos::TeamPolicy<ExecSpace, Kokkos::LaunchBounds<kPreferredTeamSize, kMinTeamsPerSM>>;
    using TeamMember = typename TeamPolicyType::member_type;
    using ScratchViewInt = Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace::scratch_memory_space,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    int const team_size = std::min<int>(kPreferredTeamSize, ExecSpace::concurrency());
    // Several consecutive elements per team, fewer on small meshes so that the league still fills the device.
    constexpr int kMaxElemsPerTeam = 8;
    constexpr int kMinWaves = 4;
    int const resident_teams = std::max<int>(1, ExecSpace::concurrency() / kPreferredTeamSize);
    // J00, J11 and the eight vertex z; the last entry of the geometry table is unused in elastic.
    constexpr int kGeomEntries = kZGeomStride - 1;
    bool const pipelined = team_size >= std::max(kPointsPerElement, kGeomEntries);
    int const elems_per_team = pipelined ? std::clamp(n_iter / (kMinWaves * resident_teams), 1, kMaxElemsPerTeam) : 1;
    TeamPolicyType policy((n_iter + elems_per_team - 1) / elems_per_team, team_size);
    // One float block: fields (then forces), fluxes, basis table and element geometry. The geometry is
    // read by the whole team as three 16-byte blocks, so it starts on a 16-byte boundary.
    constexpr int kFluxOffset = kNumFields * kPointsPerElement;
    constexpr int kTabOffset = kFluxOffset + 9 * kPointsPerElement;
    constexpr int kGeomOffset = ((kTabOffset + INTEGRAL_TYPE::kBasisTableSize + 3) / 4) * 4;
    constexpr int kGeomPadded = ((kGeomEntries + 3) / 4) * 4;
    static_assert(kGeomPadded == 12, "the geometry is read as three 16-byte blocks");
    constexpr int kScratchFloats = kGeomOffset + kGeomPadded;
    constexpr size_t kScratchBytes = kScratchFloats * sizeof(float);
    size_t const bytes = kScratchBytes + alignof(detail::Float4) + ScratchViewInt::shmem_size(kPointsPerElement);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes));

    Kokkos::parallel_for(
        "Solver Element Contribution Tti TeamZ", detail::lightWeight(policy), KOKKOS_LAMBDA(const TeamMember& team) {
          float* scratch =
              static_cast<float*>(team.team_scratch(0).get_shmem_aligned(kScratchBytes, alignof(detail::Float4)));
          // Global index of each node, read once here and reused by the callback and the scatter.
          ScratchViewInt nodeIdx(team.team_scratch(0), kPointsPerElement);
          // Displacements in, forces out, same storage.
          float* localFields = scratch;
          float* fluxScratch = scratch + kFluxOffset;
          float* basisTab = scratch + kTabOffset;
          float* geomSh = scratch + kGeomOffset;
          detail::Float4 const* geom4 = reinterpret_cast<detail::Float4 const*>(geomSh);

          int const first_rank = team.league_rank() * elems_per_team;
          int const n_here = Kokkos::min(elems_per_team, n_iter - first_rank);
          int const point = team.team_rank();
          bool const owns_point = pipelined && point < kPointsPerElement;
          bool const carries_geom = pipelined && point < kGeomEntries;

          // Same two-stage pipeline as computeElementContributions_Acoustic_TeamZ().
          int ahead_node = 0;
          int ahead_e = 0;
          int next_node = 0;
          float next_fields[kNumFields] = {};
          float next_geom = 0.0f;
          auto const load_index = [&](int rank) {
            // Indexed by the league position, so the gather does not wait for the element list.
            if (owns_point) ahead_node = elemNodes[rank * kPointsPerElement + point];
            if (carries_geom) ahead_e = list_on ? list_local[rank] : rank;
          };
          auto const load_values = [&]() {
            if (owns_point) {
              next_node = ahead_node;
              for (int f = 0; f < kNumFields; ++f) next_fields[f] = data.getCurrentField(f)(next_node);
            }
            if (carries_geom) next_geom = zgeom(point * nElems + ahead_e);
          };
          if (pipelined) {
            load_index(first_rank);
            load_values();
            if (n_here > 1) load_index(first_rank + 1);
          }

          for (int k = 0; k < n_here; ++k) {
            int const rank = first_rank + k;

            if (k == 0) {
              // Filled before the gather barrier, so it needs no barrier of its own.
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, INTEGRAL_TYPE::kBasisTableSize),
                                   [&](const int idx) { basisTab[idx] = basisTabGlobal[idx]; });
            }

            if (pipelined) {
              if (owns_point) {
                nodeIdx(point) = next_node;
                for (int f = 0; f < kNumFields; ++f) localFields[f * kPointsPerElement + point] = next_fields[f];
              }
              if (carries_geom) geomSh[point] = next_geom;
            } else {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
                int const globalIdx = elemNodes[rank * kPointsPerElement + localIdx];
                nodeIdx(localIdx) = globalIdx;
                for (int f = 0; f < kNumFields; ++f)
                  localFields[f * kPointsPerElement + localIdx] = data.getCurrentField(f)(globalIdx);
              });
              int const elementNumber = list_on ? list_local[rank] : rank;
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kGeomEntries),
                                   [&](const int j) { geomSh[j] = zgeom(j * nElems + elementNumber); });
            }
            team.team_barrier();

            if (k + 1 < n_here) load_values();
            if (k + 2 < n_here) load_index(rank + 2);

            // Entries J00, J11 and the eight vertex z, in that order.
            detail::Float4 const g0 = geom4[0];
            detail::Float4 const g1 = geom4[1];
            detail::Float4 const g2 = geom4[2];
            float const Z[8] = {g0.v[2], g0.v[3], g1.v[0], g1.v[1], g1.v[2], g1.v[3], g2.v[0], g2.v[1]};

            auto const tti_flux = [&](int qa, int qb, int qc, float const(&J_inv)[3][3], float const(&grad_u_ref)[3][3],
                                      float(&flux)[3][3]) {
              int const gIndex = nodeIdx(qa + qb * dim + qc * dim * dim);
              float p[flux::kTtiCompactSize];
              detail::Float4 const lo = ctti_quads[2 * static_cast<size_t>(gIndex)];
              detail::Float4 const hi = ctti_quads[2 * static_cast<size_t>(gIndex) + 1];
              for (int j = 0; j < 4; ++j) {
                p[j] = lo.v[j];
                p[j + 4] = hi.v[j];
              }
              flux::elasticFluxTtiCompact(J_inv, p, grad_u_ref, flux);
            };

            INTEGRAL_TYPE::computeElasticStiffnessSumFactTeamZDeformed(team, g0.v[0], g0.v[1], Z, localFields,
                                                                       localFields, fluxScratch, tti_flux, basisTab);
            team.team_barrier();

            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, kPointsPerElement), [&](const int localIdx) {
              int const globalIdx = nodeIdx(localIdx);
              for (int f = 0; f < kNumFields; ++f) {
                ATOMICADD(local_workVectorsGlobal[f][globalIdx], localFields[f * kPointsPerElement + localIdx]);
              }
            });
          }
        });
  }
}


}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_IMPL_COMMON_INCLUDE_SEM_SOLVER_KERNELS_H_
