#ifndef FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_INIT_H_
#define FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_INIT_H_

/**
 * @file sem_solver_init.h
 * @brief Setup and initialisation of SEMsolver: global mass and damping
 * matrices, array allocation, sponge-layer values, diagnostic output.
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
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeGlobalMassMatrix() {
  auto mesh_local = m_mesh;
  auto local_massMatrixGlobal = massMatrixGlobal_;
  bool const mask_enabled = m_mask_enabled_;
  auto element_mask = m_element_mask_;
  int const mask_active_value = m_mask_active_value_;

  Kokkos::parallel_for(
      "Solver Compute GMatrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (mask_enabled && element_mask[elementNumber] != mask_active_value) return;
        float massMatrixLocal[kPointsPerElement] = {0};
        int const dim = mesh_local.getOrder() + 1;

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh_local.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh_local.vertexCoords(mesh_local.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        INTEGRAL_TYPE::computeMassTerm(cornerCoords, [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

        // Acoustic: 1 / (vp^2 rho). Elastic: rho.
        real_t model_factor = 0.0f;
        if constexpr (!IS_MODEL_ON_NODES) {
          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            model_factor =
                1.0f / (mesh_local.getModelVpOnElement(elementNumber) * mesh_local.getModelVpOnElement(elementNumber) *
                        mesh_local.getModelRhoOnElement(elementNumber));
          } else {
            model_factor = mesh_local.getModelRhoOnElement(elementNumber);
          }
        }

        for (int i = 0; i < mesh_local.getNumberOfPointsPerElement(); ++i) {
          int x = i % dim;
          int z = (i / dim) % dim;
          int y = i / (dim * dim);
          int const gIndex = mesh_local.globalNodeIndex(elementNumber, x, y, z);

          if constexpr (IS_MODEL_ON_NODES) {
            if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
              model_factor = 1.0f / (mesh_local.getModelVpOnNodes(gIndex) * mesh_local.getModelVpOnNodes(gIndex) *
                                     mesh_local.getModelRhoOnNodes(gIndex));
            } else {
              model_factor = mesh_local.getModelRhoOnNodes(gIndex);
            }
          }

          massMatrixLocal[i] *= model_factor;
          ATOMICADD(local_massMatrixGlobal[gIndex], massMatrixLocal[i]);
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::computeDampingMatrix() {
  auto mesh_local = m_mesh;
  bool const mask_enabled = m_mask_enabled_;
  auto element_mask = m_element_mask_;
  int const mask_active_value = m_mask_active_value_;

  std::array<std::remove_reference_t<decltype(dampingMatrixGlobal_[0])>, kNumFields> local_dampingMatrixGlobal;
  for (int f = 0; f < kNumFields; ++f) {
    local_dampingMatrixGlobal[f] = dampingMatrixGlobal_[f];
  }

  Kokkos::parallel_for(
      "Solver Compute Damping Matrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh_local.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        (void)local_dampingMatrixGlobal;
        if (mask_enabled && element_mask[elementNumber] != mask_active_value) return;
        for (int i = 0; i < 6; ++i) {
          int f = mesh_local.getGlobalFace(elementNumber, static_cast<model::CubicFace>(i));

          // Only boundary faces carry an absorbing term.
          if (!mesh_local.isBoundaryFace(f)) continue;

          float coords[4][3];
          for (int j = 0; j < 4; ++j) {
            int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, INTEGRAL_TYPE::meshIndexToLinearIndex2D(j));
            for (int d = 0; d < 3; ++d) {
              coords[j][d] = mesh_local.nodeCoord(globalNodeIndex, d);
            }
          }

          if constexpr (PHYSICS == utils::enums::physicType::kAcoustic) {
            real_t model_rho = 0.0f;
            real_t model_vp = 0.0f;
            real_t alpha = 0.0f;

            if constexpr (!IS_MODEL_ON_NODES) {
              model_rho = mesh_local.getModelRhoOnElement(elementNumber);
              model_vp = mesh_local.getModelVpOnElement(elementNumber);
              alpha = 1.0 / (model_rho * model_vp);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Free-surface nodes are not damped.
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                model_rho = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                model_vp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                alpha = 1.0 / (model_rho * model_vp);
              }

              real_t localIncrement = alpha * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrement);
            }
          } else {
            float normal[3];
            mesh_local.faceNormal(elementNumber, static_cast<model::CubicFace>(i), normal);
            real_t nx = normal[0], ny = normal[1], nz = normal[2];

            real_t density, velocityVp, velocityVs;

            if constexpr (!IS_MODEL_ON_NODES) {
              density = mesh_local.getModelRhoOnElement(elementNumber);
              velocityVp = mesh_local.getModelVpOnElement(elementNumber);
              velocityVs = mesh_local.getModelVsOnElement(elementNumber);
            }

            constexpr int numNodesPerFace = (ORDER + 1) * (ORDER + 1);
            for (int q = 0; q < numNodesPerFace; ++q) {
              int const globalNodeIndex = mesh_local.getGlobalNodeFromFace(f, q);

              // Free-surface nodes are not damped.
              if (mesh_local.isFreeSurface(globalNodeIndex)) {
                continue;
              }

              if constexpr (IS_MODEL_ON_NODES) {
                density = mesh_local.getModelRhoOnNodes(globalNodeIndex);
                velocityVp = mesh_local.getModelVpOnNodes(globalNodeIndex);
                velocityVs = mesh_local.getModelVsOnNodes(globalNodeIndex);
              }

              real_t aux = density * INTEGRAL_TYPE::computeDampingTerm(q, coords);
              real_t localIncrementx = aux * (velocityVp * std::fabs(nx) + velocityVs * std::sqrt(ny * ny + nz * nz));
              real_t localIncrementy = aux * (velocityVp * std::fabs(ny) + velocityVs * std::sqrt(nx * nx + nz * nz));
              real_t localIncrementz = aux * (velocityVp * std::fabs(nz) + velocityVs * std::sqrt(nx * nx + ny * ny));

              ATOMICADD(local_dampingMatrixGlobal[0][globalNodeIndex], localIncrementx);
              ATOMICADD(local_dampingMatrixGlobal[1][globalNodeIndex], localIncrementy);
              ATOMICADD(local_dampingMatrixGlobal[2][globalNodeIndex], localIncrementz);
            }
          }
        }
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::allocateFEarrays() {
  massMatrixGlobal_ = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), "massMatrixGlobal");

  static constexpr const char* dampingNames[3] = {"dampingX", "dampingY", "dampingZ"};
  for (int f = 0; f < kNumFields; ++f) {
    dampingMatrixGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), dampingNames[f]);
  }

  static constexpr const char* workVectorNames[3] = {"workVec0", "workVec1", "workVec2"};
  for (int f = 0; f < kNumFields; ++f) {
    workVectorsGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), workVectorNames[f]);
  }

  if (attenuationEnabled_ && nSls_ > 0) {
    static constexpr const char* attWorkNames[3] = {"attWorkVec0", "attWorkVec1", "attWorkVec2"};
    static constexpr const char* attMemNames[3] = {"attMemory0", "attMemory1", "attMemory2"};
    for (int f = 0; f < kNumFields; ++f) {
      attenuationWorkVectorsGlobal_[f] = allocateVector<vectorReal>(m_mesh.getNumberOfNodes(), attWorkNames[f]);
      attenuationMemoryVariables_[f] = allocateArray2D<arrayReal>(m_mesh.getNumberOfNodes(), nSls_, attMemNames[f]);
    }
  }

  spongeTaperCoeff_ = allocateVector<dualVectorReal>(m_mesh.getNumberOfNodes(), "spongeTaperCoeff");
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initFEarrays() {
  initSpongeValues();

  if (attenuationEnabled_ && nSls_ > 0) {
    // Device-only work arrays: zeroed in place.
    for (int f = 0; f < kNumFields; ++f) {
      Kokkos::deep_copy(attenuationWorkVectorsGlobal_[f], 0.0f);
      Kokkos::deep_copy(attenuationMemoryVariables_[f], 0.0f);
    }
    FENCE
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::initSpongeValues() {
  // @todo VERIFY: unit and origin of the peak damping factor sigma_max = 0.15 (dimensionless?).
  const double sigma_max = 0.15;

  // Computed on the host (the mesh accessors read its host copy), then copied to the device.
  auto taper_host = spongeTaperCoeff_.view_host();
  for (int n = 0; n < m_mesh.getNumberOfNodes(); n++) {
    const double x = m_mesh.nodeCoord(n, 0);
    const double y = m_mesh.nodeCoord(n, 1);
    const double z = m_mesh.nodeCoord(n, 2);

    const double distToFrontierX = (surface_sponge_) ? m_mesh.domainSize(0) - x : std::min(m_mesh.domainSize(0) - x, x);
    const double distToFrontierY = std::min(m_mesh.domainSize(1) - y, y);
    const double distToFrontierZ = std::min(m_mesh.domainSize(2) - z, z);

    double minDistToFrontier = std::max(m_mesh.domainSize(0), std::max(m_mesh.domainSize(1), m_mesh.domainSize(2)));

    bool is_sponge = false;
    if (distToFrontierX < sponge_size_[0]) {
      is_sponge = true;
      minDistToFrontier = std::min(minDistToFrontier, distToFrontierX);
    }
    if (distToFrontierY < sponge_size_[1]) {
      is_sponge = true;
      minDistToFrontier = std::min(minDistToFrontier, distToFrontierY);
    }
    if (distToFrontierZ < sponge_size_[2]) {
      is_sponge = true;
      minDistToFrontier = std::min(minDistToFrontier, distToFrontierZ);
    }

    // Gaussian taper in the distance to the closest sponge frontier.
    if (is_sponge) {
      double d = minDistToFrontier;
      double delta = taper_delta_;
      double sigma = sigma_max * std::exp(-((d / delta) * (d / delta)));
      taper_host(n) = 1.0 / (1.0 + sigma);
    } else {
      taper_host(n) = 1.0;
    }
  }

  spongeTaperCoeff_.modify_host();
  spongeTaperCoeff_.sync_device();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES, physicType PHYSICS>
void SEMsolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES, PHYSICS>::outputSolutionValues(
    const int& t, int& e, const vectorReal& fieldGlobal, const char* fieldName) {
  // fieldGlobal lives on the device: only the printed value is copied to the host.
  float value = 0.0f;
  Kokkos::deep_copy(value, Kokkos::subview(fieldGlobal, m_mesh.globalNodeIndex(e, 0, 0, 0)));
  std::cout << "TimeStep=" << t << ";  " << fieldName << " @ elementSource location " << e
            << " after computeOneStep = " << value << std::endl;
}


}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_SEM_COMMON_INCLUDE_SEM_SOLVER_INIT_H_
