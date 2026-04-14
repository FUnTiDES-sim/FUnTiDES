#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_IMPL_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_IMPL_H_

#include <Kokkos_Core.hpp>
#include <typeinfo>

#include "differentiator_acoustic.h"

namespace gradient
{

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<
    ORDER, INTEGRAL_TYPE, MESH_TYPE,
    IS_MODEL_ON_NODES>::compute(model::ModelApi<float, int>& mesh,
                                DataStruct& data, float dt) const
{
  auto& myData = dynamic_cast<DifferentiatorDataAcoustic&>(data);
  auto& myMesh = dynamic_cast<MESH_TYPE&>(mesh);

  VECTOR_REAL_VIEW const pn = myData.getForwardField(0);
  VECTOR_REAL_VIEW const qn = myData.getBackwardField(0);
  VECTOR_REAL_VIEW const qnPrev = myData.getBackwardField(1);
  VECTOR_REAL_VIEW const qnPrevPrev = myData.getBackwardField(2);
  VECTOR_REAL_VIEW const gradKappa = myData.getGradient(0);
  VECTOR_REAL_VIEW const gradBuoyancy = myData.getGradient(1);

  if constexpr (!IS_MODEL_ON_NODES)
    computeOnElements(myMesh, dt, pn, qn, qnPrev, qnPrevPrev, gradKappa,
                      gradBuoyancy);
  else
    computeOnNodes(myMesh, dt, pn, qn, qnPrev, qnPrevPrev, gradKappa,
                   gradBuoyancy);
  Kokkos::fence();
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
int DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                           IS_MODEL_ON_NODES>::getOrder() const
{
  return kOrder;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
bool DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                            IS_MODEL_ON_NODES>::isModelOnNodes() const
{
  return kIsModelOnNodes;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                            IS_MODEL_ON_NODES>::print() const
{
  std::cout << "DifferentiatorAcoustic<ORDER=" << kOrder
            << ", INTEGRAL_TYPE=" << typeid(INTEGRAL_TYPE).name()
            << ", MESH_TYPE=" << typeid(MESH_TYPE).name()
            << ", IS_MODEL_ON_NODES=" << (kIsModelOnNodes ? "true" : "false")
            << ">\n";
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                            IS_MODEL_ON_NODES>::
    computeOnElements(MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn,
                      VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qnPrev,
                      VECTOR_REAL_VIEW const qnPrevPrev,
                      VECTOR_REAL_VIEW const gradKappa,
                      VECTOR_REAL_VIEW const gradBuoyancy) const
{
  constexpr int nPerElem = kPointsPerElement;
  float const invDt2 = 1.0f / (dt * dt);

  Kokkos::parallel_for(
      "Compute Acoustic Gradient on Elements",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh.getNumberOfElements()) return;

        int const dim = mesh.getOrder() + 1;

        typename INTEGRAL_TYPE::TransformType transformData;
        {
          auto const elementIndex = mesh.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
              {
                auto const vertexIndex =
                    mesh.globalVertexIndex(elementIndex, iv, jv, kv);
                mesh.vertexCoords(vertexIndex, transformData.data[I]);
                ++I;
              }
        }

        Kokkos::Array<float, nPerElem> localPn, localQn, localQnPrev,
            localQnPrevPrev;

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k)
            {
              int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
              int const lIdx = i + j * dim + k * dim * dim;
              localPn[lIdx] = pn(gIdx);
              localQn[lIdx] = qn(gIdx);
              localQnPrev[lIdx] = qnPrev(gIdx);
              localQnPrevPrev[lIdx] = qnPrevPrev(gIdx);
            }

        // Accumulation directe dans les Views (pas de conflits entre éléments)
        INTEGRAL_TYPE::computeMassTerm(
            transformData, [=](const int q, const real_t val) {
              float const qdt2 =
                  (localQnPrevPrev[q] - 2.0f * localQnPrev[q] + localQn[q]) *
                  invDt2;
              gradKappa(elementNumber) += qdt2 * localPn[q] * val;
            });

        INTEGRAL_TYPE::computeStiffnessTerm(
            transformData, [=](const int, const int, const int) {},
            [=](const int i, const int j, const real_t val) {
              gradBuoyancy(elementNumber) += val * localQn[j] * localPn[i];
            });
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorAcoustic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                            IS_MODEL_ON_NODES>::
    computeOnNodes(MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn,
                   VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qnPrev,
                   VECTOR_REAL_VIEW const qnPrevPrev,
                   VECTOR_REAL_VIEW const gradKappa,
                   VECTOR_REAL_VIEW const gradBuoyancy) const
{
  constexpr int nPerElem = kPointsPerElement;
  float const invDt2 = 1.0f / (dt * dt);

  Kokkos::parallel_for(
      "Compute Acoustic Gradient on Nodes",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh.getNumberOfElements()) return;

        int const dim = mesh.getOrder() + 1;

        typename INTEGRAL_TYPE::TransformType transformData;
        {
          auto const elementIndex = mesh.elementIndex(elementNumber);
          int I = 0;
          for (int kv = 0; kv < 2; ++kv)
            for (int jv = 0; jv < 2; ++jv)
              for (int iv = 0; iv < 2; ++iv)
              {
                auto const vertexIndex =
                    mesh.globalVertexIndex(elementIndex, iv, jv, kv);
                mesh.vertexCoords(vertexIndex, transformData.data[I]);
                ++I;
              }
        }

        Kokkos::Array<float, nPerElem> localPn, localQn, localQnPrev,
            localQnPrevPrev;
        Kokkos::Array<int, nPerElem> localGIdx;

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k)
            {
              int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
              int const lIdx = i + j * dim + k * dim * dim;
              localGIdx[lIdx] = gIdx;
              localPn[lIdx] = pn(gIdx);
              localQn[lIdx] = qn(gIdx);
              localQnPrev[lIdx] = qnPrev(gIdx);
              localQnPrevPrev[lIdx] = qnPrevPrev(gIdx);
            }

        INTEGRAL_TYPE::computeMassTerm(
            transformData, [=](const int q, const real_t val) {
              float const qdt2 =
                  (localQnPrevPrev[q] - 2.0f * localQnPrev[q] + localQn[q]) *
                  invDt2;
              ATOMICADD(gradKappa(localGIdx[q]), qdt2 * localPn[q] * val);
            });

        INTEGRAL_TYPE::computeStiffnessTerm(
            transformData, [=](const int, const int, const int) {},
            [=](const int i, const int j, const real_t val) {
              ATOMICADD(gradBuoyancy(localGIdx[i]),
                        val * localQn[j] * localPn[i]);
            });
      });
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_IMPL_H_
