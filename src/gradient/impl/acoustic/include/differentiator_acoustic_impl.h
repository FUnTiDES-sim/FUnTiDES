#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_IMPL_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_IMPL_H_

#include <Kokkos_Core.hpp>
#include <typeinfo>

#include "differentiator_acoustic.h"

namespace gradient {

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::compute(
    model::ModelApi<float, int>& mesh, DataStruct& data, float dt) const {
  auto& myData = dynamic_cast<DifferentiatorDataAcoustic&>(data);
  auto& myMesh = dynamic_cast<MESH_TYPE&>(mesh);

  VECTOR_REAL_VIEW const pn = myData.getForwardField(0);
  VECTOR_REAL_VIEW const qn = myData.getBackwardField(0);
  VECTOR_REAL_VIEW const qnPrev = myData.getBackwardField(1);
  VECTOR_REAL_VIEW const qnPrevPrev = myData.getBackwardField(2);
  VECTOR_REAL_VIEW const gradKappa = myData.getGradient(0);
  VECTOR_REAL_VIEW const gradBuoyancy = myData.getGradient(1);

  if constexpr (!IS_MODEL_ON_NODES)
    computeOnElements(myMesh, dt, pn, qn, qnPrev, qnPrevPrev, gradKappa, gradBuoyancy);
  else
    computeOnNodes(myMesh, dt, pn, qn, qnPrev, qnPrevPrev, gradKappa, gradBuoyancy);
  Kokkos::fence();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
int DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::getOrder() const {
  return kOrder;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
bool DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::isModelOnNodes() const {
  return kIsModelOnNodes;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::print() const {
  std::cout << "DifferentiatorAcoustic<ORDER=" << kOrder << ", INTEGRAL_TYPE=" << typeid(INTEGRAL_TYPE).name()
            << ", MESH_TYPE=" << typeid(MESH_TYPE).name()
            << ", IS_MODEL_ON_NODES=" << (kIsModelOnNodes ? "true" : "false") << ">\n";
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeOnElements(
    MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn, VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qnPrev,
    VECTOR_REAL_VIEW const qnPrevPrev, VECTOR_REAL_VIEW const gradKappa, VECTOR_REAL_VIEW const gradBuoyancy) const {
  constexpr int nPerElem = kPointsPerElement;
  float const invDt2 = 1.0f / (dt * dt);

  Kokkos::parallel_for(
      "Compute Acoustic Gradient on Elements",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh.getNumberOfElements()) return;

        int const dim = mesh.getOrder() + 1;

        float X[8][3];
        {
          auto const elementIndex = mesh.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv) {
                auto const vertexIndex = mesh.globalVertexIndex(elementIndex, iv, jv, kv);
                mesh.vertexCoords(vertexIndex, X[I]);
                ++I;
              }
        }

        Kokkos::Array<float, nPerElem> localPn, localQn, localQnPrev, localQnPrevPrev;

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
              int const lIdx = i + j * dim + k * dim * dim;
              localPn[lIdx] = pn(gIdx);
              localQn[lIdx] = qn(gIdx);
              localQnPrev[lIdx] = qnPrev(gIdx);
              localQnPrevPrev[lIdx] = qnPrevPrev(gIdx);
            }

        float const invDt2 = 1.0f / (dt * dt);

        float localGradKappa = 0.0f;
        INTEGRAL_TYPE::computeMassTerm(X, [&](const int q, const real_t val) {
          float const qdt2 = (localQnPrevPrev[q] - 2.0f * localQnPrev[q] + localQn[q]) * invDt2;
          localGradKappa += qdt2 * localPn[q] * val;
        });
        gradKappa(elementNumber) += localGradKappa;

        float localGradBuoyancy = 0.0f;
        INTEGRAL_TYPE::computeStiffnessTerm(
            X, [&](const int /*qa*/, const int /*qb*/, const int /*qc*/) {},
            [&](const int i, const int j, const real_t val) { localGradBuoyancy += val * localQn[j] * localPn[i]; });
        gradBuoyancy(elementNumber) += localGradBuoyancy;
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeOnNodes(
    MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn, VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qnPrev,
    VECTOR_REAL_VIEW const qnPrevPrev, VECTOR_REAL_VIEW const gradKappa, VECTOR_REAL_VIEW const gradBuoyancy) const {
  // Get number of nodes
  int const nNodes = mesh.getNumberOfNodes();

  // Allocate mass matrix diagonal
  VECTOR_REAL_VIEW massDiag = VECTOR_REAL_VIEW("massDiag", nNodes);

  // Local copy of class members to avoid capturing 'this' on the device
  constexpr int nPerElem = kPointsPerElement;

  // =====================================================
  // Compute mass matrix diagonal
  // =====================================================
  // TODO: This is currently computed on the fly within computeOnNodes, but it
  // could be precomputed and reused across iterations for efficiency.
  Kokkos::parallel_for(
      "Compute Mass Matrix Diagonal",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh.getNumberOfElements()) return;

        int const dim = mesh.getOrder() + 1;

        float X[8][3];
        {
          auto const elementIndex = mesh.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv) {
                auto const vertexIndex = mesh.globalVertexIndex(elementIndex, iv, jv, kv);
                mesh.vertexCoords(vertexIndex, X[I]);
                ++I;
              }
        }

        int localGIdx[nPerElem] = {0};
        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
              int const lIdx = i + j * dim + k * dim * dim;
              localGIdx[lIdx] = gIdx;
            }

        // Accumulate mass matrix diagonal: M_ii = sum_K val_i
        INTEGRAL_TYPE::computeMassTerm(X,
                                       [&](const int q, const real_t val) { ATOMICADD(massDiag(localGIdx[q]), val); });
      });

  Kokkos::fence();  // Ensure mass diagonal computation is complete before
                    // proceeding

  // =====================================================
  // Compute gradients and normalize
  // =====================================================
  Kokkos::parallel_for(
      "Compute and Distribute Element Gradients to Nodes",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh.getNumberOfElements()) return;

        int const dim = mesh.getOrder() + 1;

        float X[8][3];
        {
          auto const elementIndex = mesh.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv) {
                auto const vertexIndex = mesh.globalVertexIndex(elementIndex, iv, jv, kv);
                mesh.vertexCoords(vertexIndex, X[I]);
                ++I;
              }
        }

        float localPn[kPointsPerElement] = {0};
        float localQn[kPointsPerElement] = {0};
        float localQnPrev[kPointsPerElement] = {0};
        float localQnPrevPrev[kPointsPerElement] = {0};
        int localGIdx[kPointsPerElement] = {0};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k) {
              int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
              int const lIdx = i + j * dim + k * dim * dim;
              localGIdx[lIdx] = gIdx;
              localPn[lIdx] = pn(gIdx);
              localQn[lIdx] = qn(gIdx);
              localQnPrev[lIdx] = qnPrev(gIdx);
              localQnPrevPrev[lIdx] = qnPrevPrev(gIdx);
            }

        float const invDt2 = 1.0f / (dt * dt);

        // =====================================================
        // Compute element gradient for kappa
        // =====================================================
        float localGradKappa = 0.0f;
        INTEGRAL_TYPE::computeMassTerm(X, [&](const int q, const real_t val) {
          float const qdt2 = (localQnPrevPrev[q] - 2.0f * localQnPrev[q] + localQn[q]) * invDt2;
          localGradKappa += qdt2 * localPn[q] * val;
        });

        // =====================================================
        // Distribute kappa gradient to nodes, weighted by mass matrix
        // =====================================================
        INTEGRAL_TYPE::computeMassTerm(X, [&](const int q, const real_t val) {
          int const gIdx = localGIdx[q];
          // Weight: local mass value / global mass diagonal (ensures proper
          // normalization)
          float const weight = val / massDiag(gIdx);
          float const contrib = localGradKappa * weight;
          ATOMICADD(gradKappa(gIdx), contrib);
        });

        // =====================================================
        // Compute element gradient for buoyancy
        // =====================================================
        float localGradBuoyancy = 0.0f;
        INTEGRAL_TYPE::computeStiffnessTerm(
            X, [&](const int /*qa*/, const int /*qb*/, const int /*qc*/) {},
            [&](const int i, const int j, const real_t val) { localGradBuoyancy += val * localQn[j] * localPn[i]; });

        // =====================================================
        // Distribute buoyancy gradient to nodes
        // For buoyancy, use the mass matrix weights of the test function nodes
        // =====================================================
        INTEGRAL_TYPE::computeMassTerm(X, [&](const int q, const real_t val) {
          int const gIdx = localGIdx[q];
          float const weight = val / massDiag(gIdx);
          float const contrib = localGradBuoyancy * weight;
          ATOMICADD(gradBuoyancy(gIdx), contrib);
        });
      });
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_IMPL_H_
