#ifndef FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_
#define FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_

#include <iostream>

#include "differentiator.h"
#include "differentiator_data_acoustic.h"
#include "model.h"

namespace gradient
{

/**
 * @brief Acoustic gradient computation for independent use.
 *
 * Computes model parameter gradients (grad_kappa, grad_buoyancy) from acoustic
 * forward and adjoint wavefields. Completely independent from the Solver.
 *
 * Features:
 * - Supports both node-based and element-based model discretization
 * - Uses standard SEM assembly with mass and stiffness matrices
 *
 * Template Parameters:
 *   ORDER                 - Polynomial order (1, 2, 3, ...)
 *   INTEGRAL_TYPE         - Integration kernel (e.g., makutu)
 *   MESH_TYPE             - Mesh topology (e.g., Cartesian)
 *   IS_MODEL_ON_NODES     - Model discretization (true=nodes, false=elements)
 */
template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
class DifferentiatorAcoustic : public Differentiator
{
 public:

  static constexpr int kOrder = ORDER;
  static constexpr bool kIsModelOnNodes = IS_MODEL_ON_NODES;
  static constexpr int kPointsPerElement = (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  ~DifferentiatorAcoustic() override = default;

  /**
   * @brief Compute acoustic gradients (Kappa, Buoyancy).
   *
   * Computes:
   *   grad_kappa   = ∑_elements ∑_quadrature q_dt² * p * mass_term
   *   grad_buoyancy = ∑_elements ∑_stiffness stiffness_term * q * p
   */
  void compute(model::ModelApi<float, int>& mesh,
               DataStruct& data) const override
  {
    auto& myData =
        dynamic_cast<DifferentiatorDataAcoustic&>(data);
    auto mesh_copy = dynamic_cast<MESH_TYPE&>(
        mesh);  // value copy for device capture TODO check that

    VECTOR_REAL_VIEW const pn =
        myData.getForwardField(0);  // forward pressure, node-indexed
    VECTOR_REAL_VIEW const qn =
        myData.getBackwardField(0);  // adjoint pressure, node-indexed
    VECTOR_REAL_VIEW const qdt2 = myData.getBackwardField(1);  // d²q/dt², node-indexed
    VECTOR_REAL_VIEW const gradKappa =
        myData.getGradient(0);  // grad_kappa, node- or element-indexed
    VECTOR_REAL_VIEW const gradBuoyancy =
        myData.getGradient(1);  // grad_buoyancy, node- or element-indexed

    if constexpr (!IS_MODEL_ON_NODES)
      computeOnElements(mesh_copy, pn, qn, qdt2, gradKappa, gradBuoyancy);
    else
      computeOnNodes(mesh_copy, pn, qn, qdt2, gradKappa, gradBuoyancy);
  }

  int getOrder() const override { return kOrder; }

  bool isModelOnNodes() const override { return kIsModelOnNodes; }

  void print() const override
  {
    std::cout << "DifferentiatorAcoustic<ORDER=" << kOrder
              << ", INTEGRAL_TYPE=" << typeid(INTEGRAL_TYPE).name()
              << ", MESH_TYPE=" << typeid(MESH_TYPE).name()
              << ", IS_MODEL_ON_NODES=" << (kIsModelOnNodes ? "true" : "false")
              << ">\n";
  }

 private:

  /**
   * @brief Each element writes to a unique index — no atomic add required.
   */
  void computeOnElements(MESH_TYPE mesh_copy, VECTOR_REAL_VIEW const pn,
                         VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qdt2,
                         VECTOR_REAL_VIEW const gradKappa,
                         VECTOR_REAL_VIEW const gradBuoyancy) const
  {
    MAINLOOPHEAD(mesh_copy.getNumberOfElements(), elementNumber)

    if (elementNumber >= mesh_copy.getNumberOfElements()) return;

    int const dim = mesh_copy.getOrder() + 1;

    typename INTEGRAL_TYPE::TransformType transformData;
    {
      auto const elementIndex = mesh_copy.elementIndex(elementNumber);
      int I = 0;
      for (int kv = 0; kv < 2; ++kv)
        for (int jv = 0; jv < 2; ++jv)
          for (int iv = 0; iv < 2; ++iv)
          {
            auto const vertexIndex =
                mesh_copy.globalVertexIndex(elementIndex, iv, jv, kv);
            mesh_copy.vertexCoords(vertexIndex, transformData.data[I]);
            ++I;
          }
    }

    float localPn[kPointsPerElement] = {0};
    float localQn[kPointsPerElement] = {0};
    float localQdt2[kPointsPerElement] = {0};
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
        {
          int const gIdx = mesh_copy.globalNodeIndex(elementNumber, i, j, k);
          int const lIdx = i + j * dim + k * dim * dim;
          localPn[lIdx] = pn(gIdx);
          localQn[lIdx] = qn(gIdx);
          localQdt2[lIdx] = qdt2(gIdx);
        }

    // grad_kappa: element-indexed — unique per thread, no atomic
    float localGradKappa = 0.0f;
    INTEGRAL_TYPE::computeMassTerm(
        transformData, [&](const int q, const real_t val) {
          localGradKappa += localQdt2[q] * localPn[q] * val;
        });
    gradKappa(elementNumber) += localGradKappa;

    // grad_buoyancy: element-indexed — unique per thread, no atomic
    float localGradBuoyancy = 0.0f;
    INTEGRAL_TYPE::computeStiffnessTerm(
        transformData,
        [&](const int /*qa*/, const int /*qb*/, const int /*qc*/) {},
        [&](const int i, const int j, const real_t val) {
          localGradBuoyancy += val * localQn[j] * localPn[i];
        });
    gradBuoyancy(elementNumber) += localGradBuoyancy;

    MAINLOOPEND
  }

  /**
   * @brief Multiple elements share boundary nodes — ATOMICADD required.
   */
  void computeOnNodes(MESH_TYPE mesh_copy, VECTOR_REAL_VIEW const pn,
                      VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qdt2,
                      VECTOR_REAL_VIEW const gradKappa,
                      VECTOR_REAL_VIEW const gradBuoyancy) const
  {
    MAINLOOPHEAD(mesh_copy.getNumberOfElements(), elementNumber)

    if (elementNumber >= mesh_copy.getNumberOfElements()) return;

    int const dim = mesh_copy.getOrder() + 1;

    typename INTEGRAL_TYPE::TransformType transformData;
    {
      auto const elementIndex = mesh_copy.elementIndex(elementNumber);
      int I = 0;
      for (int kv = 0; kv < 2; ++kv)
        for (int jv = 0; jv < 2; ++jv)
          for (int iv = 0; iv < 2; ++iv)
          {
            auto const vertexIndex =
                mesh_copy.globalVertexIndex(elementIndex, iv, jv, kv);
            mesh_copy.vertexCoords(vertexIndex, transformData.data[I]);
            ++I;
          }
    }

    float localPn[kPointsPerElement] = {0};
    float localQn[kPointsPerElement] = {0};
    float localQdt2[kPointsPerElement] = {0};
    int localGIdx[kPointsPerElement] = {0};
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
        {
          int const gIdx = mesh_copy.globalNodeIndex(elementNumber, i, j, k);
          int const lIdx = i + j * dim + k * dim * dim;
          localGIdx[lIdx] = gIdx;
          localPn[lIdx] = pn(gIdx);
          localQn[lIdx] = qn(gIdx);
          localQdt2[lIdx] = qdt2(gIdx);
        }

    // grad_kappa: scatter per quadrature point to its global node
    INTEGRAL_TYPE::computeMassTerm(
        transformData, [&](const int q, const real_t val) {
          ATOMICADD(gradKappa(localGIdx[q]), localQdt2[q] * localPn[q] * val);
        });

    // grad_buoyancy: scatter test-function node (i) contributions to global
    // node
    INTEGRAL_TYPE::computeStiffnessTerm(
        transformData,
        [&](const int /*qa*/, const int /*qb*/, const int /*qc*/) {},
        [&](const int i, const int j, const real_t val) {
          ATOMICADD(gradBuoyancy(localGIdx[i]), val * localQn[j] * localPn[i]);
        });

    MAINLOOPEND
  }

};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_
