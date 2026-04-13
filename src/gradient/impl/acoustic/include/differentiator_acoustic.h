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
  static constexpr int kPointsPerElement =
      (ORDER + 1) * (ORDER + 1) * (ORDER + 1);

  ~DifferentiatorAcoustic() override = default;

  /**
   * @brief Compute acoustic gradients (Kappa, Buoyancy).
   *
   * Computes:
   *   grad_kappa   = ∑_elements ∑_quadrature q_dt² * p * mass_term
   *   grad_buoyancy = ∑_elements ∑_stiffness stiffness_term * q * p
   */
  void compute(model::ModelApi<float, int>& mesh, DataStruct& data,
               float dt) const override
  {
    auto& myData = dynamic_cast<DifferentiatorDataAcoustic&>(data);
    auto& myMesh = dynamic_cast<MESH_TYPE&>(mesh);

    VECTOR_REAL_VIEW const pn =
        myData.getForwardField(0);  // forward pressure, node-indexed
    VECTOR_REAL_VIEW const qn =
        myData.getBackwardField(0);  // adjoint pressure at n, node-indexed
    VECTOR_REAL_VIEW const qnPrev =
        myData.getBackwardField(1);  // adjoint pressure at n-1, node-indexed
    VECTOR_REAL_VIEW const qnPrevPrev =
        myData.getBackwardField(2);  // adjoint pressure at n-2, node-indexed
    VECTOR_REAL_VIEW const gradKappa =
        myData.getGradient(0);  // grad_kappa, node- or element-indexed
    VECTOR_REAL_VIEW const gradBuoyancy =
        myData.getGradient(1);  // grad_buoyancy, node- or element-indexed

    if constexpr (!IS_MODEL_ON_NODES)
      computeOnElements(myMesh, dt, pn, qn, qnPrev, qnPrevPrev, gradKappa,
                        gradBuoyancy);
    else
      computeOnNodes(myMesh, dt, pn, qn, qnPrev, qnPrevPrev, gradKappa,
                     gradBuoyancy);
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
   *
   * Computes qdt2 = (qnPrevPrev - 2*qnPrev + qn) / dt² on the fly.
   */
  void computeOnElements(MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn,
                         VECTOR_REAL_VIEW const qn,
                         VECTOR_REAL_VIEW const qnPrev,
                         VECTOR_REAL_VIEW const qnPrevPrev,
                         VECTOR_REAL_VIEW const gradKappa,
                         VECTOR_REAL_VIEW const gradBuoyancy) const
  {
    Kokkos::parallel_for(
        "Compute Acoustic Gradient on Elements",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock,
                                                 LaunchMinBlocksPerSM>>(
            0, mesh.getNumberOfElements()),
        KOKKOS_CLASS_LAMBDA(const int elementNumber) {
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

          float localPn[kPointsPerElement] = {0};
          float localQn[kPointsPerElement] = {0};
          float localQnPrev[kPointsPerElement] = {0};
          float localQnPrevPrev[kPointsPerElement] = {0};
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

          float const invDt2 = 1.0f / (dt * dt);

          // grad_kappa: element-indexed — unique per thread, no atomic
          float localGradKappa = 0.0f;
          INTEGRAL_TYPE::computeMassTerm(
              transformData, [&](const int q, const real_t val) {
                float const qdt2 =
                    (localQnPrevPrev[q] - 2.0f * localQnPrev[q] + localQn[q]) *
                    invDt2;
                localGradKappa += qdt2 * localPn[q] * val;
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
        });
  }

  /**
   * @brief Multiple elements share boundary nodes — ATOMICADD required.
   *
   * Computes qdt2 = (qnPrevPrev - 2*qnPrev + qn) / dt² on the fly.
   * Gradients are normalized by the mass matrix diagonal to account for
   * multiple elements sharing nodes at the domain interior.
   */
  void computeOnNodes(MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const pn,
                      VECTOR_REAL_VIEW const qn, VECTOR_REAL_VIEW const qnPrev,
                      VECTOR_REAL_VIEW const qnPrevPrev,
                      VECTOR_REAL_VIEW const gradKappa,
                      VECTOR_REAL_VIEW const gradBuoyancy) const
  {
    // Get number of nodes
    int const nNodes =
        mesh.globalNodeIndex(mesh.getNumberOfElements() - 1, mesh.getOrder(),
                             mesh.getOrder(), mesh.getOrder()) +
        1;

    // Allocate mass matrix diagonal
    VECTOR_REAL_VIEW massDiag = VECTOR_REAL_VIEW("massDiag", nNodes);

    // =====================================================
    // Compute mass matrix diagonal
    // =====================================================
    // TODO: This is currently computed on the fly within computeOnNodes, but it
    // could be precomputed and reused across iterations for efficiency.
    Kokkos::parallel_for(
        "Compute Mass Matrix Diagonal",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock,
                                                 LaunchMinBlocksPerSM>>(
            0, mesh.getNumberOfElements()),
        KOKKOS_CLASS_LAMBDA(const int elementNumber) {
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

          int localGIdx[kPointsPerElement] = {0};
          for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
              for (int k = 0; k < dim; ++k)
              {
                int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
                int const lIdx = i + j * dim + k * dim * dim;
                localGIdx[lIdx] = gIdx;
              }

          // Accumulate mass matrix diagonal: M_ii = sum_K val_i
          INTEGRAL_TYPE::computeMassTerm(
              transformData, [&](const int q, const real_t val) {
                ATOMICADD(massDiag(localGIdx[q]), val);
              });
        });

    // =====================================================
    // Compute gradients and normalize
    // =====================================================
    Kokkos::parallel_for(
        "Compute Acoustic Gradient on Nodes (Normalized)",
        Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock,
                                                 LaunchMinBlocksPerSM>>(
            0, mesh.getNumberOfElements()),
        KOKKOS_CLASS_LAMBDA(const int elementNumber) {
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

          float localPn[kPointsPerElement] = {0};
          float localQn[kPointsPerElement] = {0};
          float localQnPrev[kPointsPerElement] = {0};
          float localQnPrevPrev[kPointsPerElement] = {0};
          int localGIdx[kPointsPerElement] = {0};
          bool isInteriorNode[kPointsPerElement] = {false};
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
                // Node is interior if not on any boundary of the element
                isInteriorNode[lIdx] =
                    (i > 0 && i < ORDER) && (j > 0 && j < ORDER) &&
                    (k > 0 && k < ORDER);
              }

          float const invDt2 = 1.0f / (dt * dt);
          bool shouldNormalize = (mesh.getNumberOfElements() > 1);

          // grad_kappa: scatter per quadrature point to its global node.
          // Normalize by M_ii only for interior nodes in multi-element meshes.
          INTEGRAL_TYPE::computeMassTerm(
              transformData, [&](const int q, const real_t val) {
                float const qdt2 =
                    (localQnPrevPrev[q] - 2.0f * localQnPrev[q] + localQn[q]) *
                    invDt2;
                int const gIdx = localGIdx[q];
                float contrib = qdt2 * localPn[q] * val;
                // Normalize only for interior nodes in multi-element meshes
                if (shouldNormalize && isInteriorNode[q])
                {
                  contrib /= massDiag(gIdx);
                }
                ATOMICADD(gradKappa(gIdx), contrib);
              });

          // grad_buoyancy: scatter test-function node (i) contributions to
          // global node. Normalize only for interior nodes in multi-element meshes.
          INTEGRAL_TYPE::computeStiffnessTerm(
              transformData,
              [&](const int /*qa*/, const int /*qb*/, const int /*qc*/) {},
              [&](const int i, const int j, const real_t val) {
                int const gIdx = localGIdx[i];
                float contrib = val * localQn[j] * localPn[i];
                // Normalize only for interior nodes in multi-element meshes
                if (shouldNormalize && isInteriorNode[i])
                {
                  contrib /= massDiag(gIdx);
                }
                ATOMICADD(gradBuoyancy(gIdx), contrib);
              });
        });
  }
};

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ACOUSTIC_INCLUDE_DIFFERENTIATOR_ACOUSTIC_H_
