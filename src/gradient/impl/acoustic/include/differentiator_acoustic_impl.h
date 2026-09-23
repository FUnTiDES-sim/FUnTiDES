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

  vectorReal const pn = myData.getForwardField(0);
  vectorReal const qn = myData.getBackwardField(0);
  vectorReal const qnPrev = myData.getBackwardField(1);
  vectorReal const qnPrevPrev = myData.getBackwardField(2);
  vectorReal const gradKappa = myData.getGradient(0);
  vectorReal const gradBuoyancy = myData.getGradient(1);

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
vectorReal& DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::getGeometricMassMatrix() {
  return geometricMassMatrix_;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::initGeometricMassMatrix(
    model::ModelApi<float, int>& meshApi) {
  auto mesh = dynamic_cast<MESH_TYPE&>(meshApi);

  if (geometricMassMatrix_.extent(0) != mesh.getNumberOfNodes())
    geometricMassMatrix_ = allocateVector<vectorReal>(mesh.getNumberOfNodes(), "geometricMassMatrix");
  Kokkos::deep_copy(geometricMassMatrix_, 0.0f);

  auto local_geometricMassMatrix = geometricMassMatrix_;

  Kokkos::parallel_for(
      "Differentiator Compute Geometric Mass Matrix",
      Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        float massMatrixLocal[kPointsPerElement] = {0};
        int const dim = mesh.getOrder() + 1;

        float cornerCoords[8][3];
        {
          auto const eIdx = mesh.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
                mesh.vertexCoords(mesh.globalVertexIndex(eIdx, iv, jv, kv), cornerCoords[I++]);
        }

        // Geometric part only: no model factors.
        INTEGRAL_TYPE::computeMassTerm(cornerCoords, [&](const int j, const real_t val) { massMatrixLocal[j] += val; });

        for (int i = 0; i < mesh.getNumberOfPointsPerElement(); ++i) {
          int x = i % dim;
          int z = (i / dim) % dim;
          int y = i / (dim * dim);
          int const gIndex = mesh.globalNodeIndex(elementNumber, x, y, z);

          ATOMICADD(local_geometricMassMatrix[gIndex], massMatrixLocal[i]);
        }
      });
  Kokkos::fence();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE, bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::computeOnElements(
    MESH_TYPE mesh, float dt, vectorReal const pn, vectorReal const qn, vectorReal const qnPrev,
    vectorReal const qnPrevPrev, vectorReal const gradKappa, vectorReal const gradBuoyancy) const {
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
    MESH_TYPE mesh, float dt, vectorReal const pn, vectorReal const qn, vectorReal const qnPrev,
    vectorReal const qnPrevPrev, vectorReal const gradKappa, vectorReal const gradBuoyancy) const {
  // The geometric mass matrix (nodal volumes) normalizes the nodal gradients.
  // It is built here on first use if initGeometricMassMatrix() was not called,
  // and cached for later calls.
  if (geometricMassMatrix_.extent(0) != mesh.getNumberOfNodes())
    const_cast<DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>*>(this)
        ->initGeometricMassMatrix(mesh);
  auto massDiag = geometricMassMatrix_;

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

        float localGradKappa = 0.0f;
        INTEGRAL_TYPE::computeMassTerm(X, [&](const int q, const real_t val) {
          float const qdt2 = (localQnPrevPrev[q] - 2.0f * localQnPrev[q] + localQn[q]) * invDt2;
          localGradKappa += qdt2 * localPn[q] * val;
        });

        // Each node receives the element gradient weighted by its local mass
        // value divided by its global mass diagonal.
        INTEGRAL_TYPE::computeMassTerm(X, [&](const int q, const real_t val) {
          int const gIdx = localGIdx[q];
          float const weight = val / massDiag(gIdx);
          float const contrib = localGradKappa * weight;
          ATOMICADD(gradKappa(gIdx), contrib);
        });

        float localGradBuoyancy = 0.0f;
        INTEGRAL_TYPE::computeStiffnessTerm(
            X, [&](const int /*qa*/, const int /*qb*/, const int /*qc*/) {},
            [&](const int i, const int j, const real_t val) { localGradBuoyancy += val * localQn[j] * localPn[i]; });

        // Same mass-based weights as for kappa.
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
