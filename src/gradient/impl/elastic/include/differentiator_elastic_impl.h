#ifndef FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_IMPL_H_
#define FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_IMPL_H_

#include <Kokkos_Core.hpp>
#include <typeinfo>

#include "differentiator_elastic.h"

namespace gradient
{

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    compute(model::ModelApi<float, int>& mesh, DataStruct& data, float dt) const
{
  auto& myData = dynamic_cast<DifferentiatorDataElastic&>(data);
  auto& myMesh = dynamic_cast<MESH_TYPE&>(mesh);

  VECTOR_REAL_VIEW const ux_fwd = myData.getForwardField(0);
  VECTOR_REAL_VIEW const uy_fwd = myData.getForwardField(1);
  VECTOR_REAL_VIEW const uz_fwd = myData.getForwardField(2);

  VECTOR_REAL_VIEW const ux_adj = myData.getBackwardField(0);
  VECTOR_REAL_VIEW const uy_adj = myData.getBackwardField(1);
  VECTOR_REAL_VIEW const uz_adj = myData.getBackwardField(2);

  VECTOR_REAL_VIEW const ux_dt2 = myData.getBackwardField(3);
  VECTOR_REAL_VIEW const uy_dt2 = myData.getBackwardField(4);
  VECTOR_REAL_VIEW const uz_dt2 = myData.getBackwardField(5);

  VECTOR_REAL_VIEW const gradRho = myData.getGradient(0);
  VECTOR_REAL_VIEW const gradLambda = myData.getGradient(1);
  VECTOR_REAL_VIEW const gradMu = myData.getGradient(2);

  if constexpr (!IS_MODEL_ON_NODES)
    computeOnElements(myMesh, dt, ux_fwd, uy_fwd, uz_fwd, ux_adj, uy_adj,
                      uz_adj, ux_dt2, uy_dt2, uz_dt2, gradRho, gradLambda,
                      gradMu);
  else
    computeOnNodes(myMesh, dt, ux_fwd, uy_fwd, uz_fwd, ux_adj, uy_adj, uz_adj,
                   ux_dt2, uy_dt2, uz_dt2, gradRho, gradLambda, gradMu);
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
int DifferentiatorElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                          IS_MODEL_ON_NODES>::getOrder() const
{
  return kOrder;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
bool DifferentiatorElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                           IS_MODEL_ON_NODES>::isModelOnNodes() const
{
  return kIsModelOnNodes;
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE,
                           IS_MODEL_ON_NODES>::print() const
{
  std::cout << "DifferentiatorElastic<ORDER=" << kOrder
            << ", INTEGRAL_TYPE=" << typeid(INTEGRAL_TYPE).name()
            << ", MESH_TYPE=" << typeid(MESH_TYPE).name()
            << ", IS_MODEL_ON_NODES=" << (kIsModelOnNodes ? "true" : "false")
            << ">\n";
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
KOKKOS_INLINE_FUNCTION void
DifferentiatorElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeDisplacementGradient(int qa, int qb, int qc, float const (&J)[3][3],
                                float const* localUx, float const* localUy,
                                float const* localUz, float (&grad)[3][3])
{
  int const dim = kOrder + 1;
  for (int c = 0; c < 3; ++c)
    for (int d = 0; d < 3; ++d) grad[c][d] = 0.0f;

  for (int n = 0; n < dim; ++n)
  {
    int const nIdx_a = n + qb * dim + qc * dim * dim;
    float const dNdxi0 = INTEGRAL_TYPE::basisGradientAt(n, qa);
    for (int d = 0; d < 3; ++d)
    {
      float const JdN = dNdxi0 * J[0][d];
      grad[0][d] += JdN * localUx[nIdx_a];
      grad[1][d] += JdN * localUy[nIdx_a];
      grad[2][d] += JdN * localUz[nIdx_a];
    }

    int const nIdx_b = qa + n * dim + qc * dim * dim;
    float const dNdxi1 = INTEGRAL_TYPE::basisGradientAt(n, qb);
    for (int d = 0; d < 3; ++d)
    {
      float const JdN = dNdxi1 * J[1][d];
      grad[0][d] += JdN * localUx[nIdx_b];
      grad[1][d] += JdN * localUy[nIdx_b];
      grad[2][d] += JdN * localUz[nIdx_b];
    }

    int const nIdx_c = qa + qb * dim + n * dim * dim;
    float const dNdxi2 = INTEGRAL_TYPE::basisGradientAt(n, qc);
    for (int d = 0; d < 3; ++d)
    {
      float const JdN = dNdxi2 * J[2][d];
      grad[0][d] += JdN * localUx[nIdx_c];
      grad[1][d] += JdN * localUy[nIdx_c];
      grad[2][d] += JdN * localUz[nIdx_c];
    }
  }
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeOnElements(
        MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const ux_fwd,
        VECTOR_REAL_VIEW const uy_fwd, VECTOR_REAL_VIEW const uz_fwd,
        VECTOR_REAL_VIEW const ux_adj, VECTOR_REAL_VIEW const uy_adj,
        VECTOR_REAL_VIEW const uz_adj, VECTOR_REAL_VIEW const ux_dt2,
        VECTOR_REAL_VIEW const uy_dt2, VECTOR_REAL_VIEW const uz_dt2,
        VECTOR_REAL_VIEW const gradRho, VECTOR_REAL_VIEW const gradLambda,
        VECTOR_REAL_VIEW const gradMu) const
{
  Kokkos::parallel_for(
      "Compute Elastic Gradient on Elements",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_CLASS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh.getNumberOfElements()) return;

        int const dim = mesh.getOrder() + 1;

        float localUxFwd[kPointsPerElement] = {0};
        float localUyFwd[kPointsPerElement] = {0};
        float localUzFwd[kPointsPerElement] = {0};
        float localUxAdj[kPointsPerElement] = {0};
        float localUyAdj[kPointsPerElement] = {0};
        float localUzAdj[kPointsPerElement] = {0};
        float localUxDt2[kPointsPerElement] = {0};
        float localUyDt2[kPointsPerElement] = {0};
        float localUzDt2[kPointsPerElement] = {0};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k)
            {
              int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
              int const lIdx = i + j * dim + k * dim * dim;
              localUxFwd[lIdx] = ux_fwd(gIdx);
              localUyFwd[lIdx] = uy_fwd(gIdx);
              localUzFwd[lIdx] = uz_fwd(gIdx);
              localUxAdj[lIdx] = ux_adj(gIdx);
              localUyAdj[lIdx] = uy_adj(gIdx);
              localUzAdj[lIdx] = uz_adj(gIdx);
              localUxDt2[lIdx] = ux_dt2(gIdx);
              localUyDt2[lIdx] = uy_dt2(gIdx);
              localUzDt2[lIdx] = uz_dt2(gIdx);
            }

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

        float strainDiv[kPointsPerElement] = {0};
        float strainEps[kPointsPerElement] = {0};

        INTEGRAL_TYPE::computeStiffNessTermwithJac(
            transformData,
            [&](int qa, int qb, int qc, float const(&J)[3][3]) {
              float grad_fwd[3][3] = {{0}};
              float grad_adj[3][3] = {{0}};
              computeDisplacementGradient(qa, qb, qc, J, localUxFwd, localUyFwd,
                                          localUzFwd, grad_fwd);
              computeDisplacementGradient(qa, qb, qc, J, localUxAdj, localUyAdj,
                                          localUzAdj, grad_adj);

              int const qIdx = qa + qb * dim + qc * dim * dim;

              float const div_fwd =
                  grad_fwd[0][0] + grad_fwd[1][1] + grad_fwd[2][2];
              float const div_adj =
                  grad_adj[0][0] + grad_adj[1][1] + grad_adj[2][2];
              strainDiv[qIdx] = div_adj * div_fwd;

              float const exx_f = grad_fwd[0][0];
              float const eyy_f = grad_fwd[1][1];
              float const ezz_f = grad_fwd[2][2];
              float const exy_f = 0.5f * (grad_fwd[0][1] + grad_fwd[1][0]);
              float const exz_f = 0.5f * (grad_fwd[0][2] + grad_fwd[2][0]);
              float const eyz_f = 0.5f * (grad_fwd[1][2] + grad_fwd[2][1]);

              float const exx_a = grad_adj[0][0];
              float const eyy_a = grad_adj[1][1];
              float const ezz_a = grad_adj[2][2];
              float const exy_a = 0.5f * (grad_adj[0][1] + grad_adj[1][0]);
              float const exz_a = 0.5f * (grad_adj[0][2] + grad_adj[2][0]);
              float const eyz_a = 0.5f * (grad_adj[1][2] + grad_adj[2][1]);

              strainEps[qIdx] =
                  2.0f * (exx_a * exx_f + eyy_a * eyy_f + ezz_a * ezz_f) +
                  4.0f * (exy_a * exy_f + exz_a * exz_f + eyz_a * eyz_f);
            },
            [&](int, int, float, const int, const int) {});

        float localGradRho = 0.0f;
        float localGradLambda = 0.0f;
        float localGradMu = 0.0f;

        INTEGRAL_TYPE::computeMassTerm(
            transformData, [&](const int q, const real_t wdetJ) {
              localGradRho += (localUxDt2[q] * localUxFwd[q] +
                               localUyDt2[q] * localUyFwd[q] +
                               localUzDt2[q] * localUzFwd[q]) *
                              wdetJ;
              localGradLambda += strainDiv[q] * wdetJ;
              localGradMu += strainEps[q] * wdetJ;
            });

        gradRho(elementNumber) += localGradRho;
        gradLambda(elementNumber) += localGradLambda;
        gradMu(elementNumber) += localGradMu;
      });
}

template <int ORDER, typename INTEGRAL_TYPE, typename MESH_TYPE,
          bool IS_MODEL_ON_NODES>
void DifferentiatorElastic<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>::
    computeOnNodes(MESH_TYPE mesh, float dt, VECTOR_REAL_VIEW const ux_fwd,
                   VECTOR_REAL_VIEW const uy_fwd, VECTOR_REAL_VIEW const uz_fwd,
                   VECTOR_REAL_VIEW const ux_adj, VECTOR_REAL_VIEW const uy_adj,
                   VECTOR_REAL_VIEW const uz_adj, VECTOR_REAL_VIEW const ux_dt2,
                   VECTOR_REAL_VIEW const uy_dt2, VECTOR_REAL_VIEW const uz_dt2,
                   VECTOR_REAL_VIEW const gradRho,
                   VECTOR_REAL_VIEW const gradLambda,
                   VECTOR_REAL_VIEW const gradMu) const
{
  Kokkos::parallel_for(
      "Compute Elastic Gradient on Nodes",
      Kokkos::RangePolicy<
          Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(
          0, mesh.getNumberOfElements()),
      KOKKOS_CLASS_LAMBDA(const int elementNumber) {
        if (elementNumber >= mesh.getNumberOfElements()) return;

        int const dim = mesh.getOrder() + 1;

        float localUxFwd[kPointsPerElement] = {0};
        float localUyFwd[kPointsPerElement] = {0};
        float localUzFwd[kPointsPerElement] = {0};
        float localUxAdj[kPointsPerElement] = {0};
        float localUyAdj[kPointsPerElement] = {0};
        float localUzAdj[kPointsPerElement] = {0};
        float localUxDt2[kPointsPerElement] = {0};
        float localUyDt2[kPointsPerElement] = {0};
        float localUzDt2[kPointsPerElement] = {0};
        int localGIdx[kPointsPerElement] = {0};

        for (int i = 0; i < dim; ++i)
          for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k)
            {
              int const gIdx = mesh.globalNodeIndex(elementNumber, i, j, k);
              int const lIdx = i + j * dim + k * dim * dim;
              localGIdx[lIdx] = gIdx;
              localUxFwd[lIdx] = ux_fwd(gIdx);
              localUyFwd[lIdx] = uy_fwd(gIdx);
              localUzFwd[lIdx] = uz_fwd(gIdx);
              localUxAdj[lIdx] = ux_adj(gIdx);
              localUyAdj[lIdx] = uy_adj(gIdx);
              localUzAdj[lIdx] = uz_adj(gIdx);
              localUxDt2[lIdx] = ux_dt2(gIdx);
              localUyDt2[lIdx] = uy_dt2(gIdx);
              localUzDt2[lIdx] = uz_dt2(gIdx);
            }

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

        float strainDiv[kPointsPerElement] = {0};
        float strainEps[kPointsPerElement] = {0};

        INTEGRAL_TYPE::computeStiffNessTermwithJac(
            transformData,
            [&](int qa, int qb, int qc, float const(&J)[3][3]) {
              float grad_fwd[3][3] = {{0}};
              float grad_adj[3][3] = {{0}};
              computeDisplacementGradient(qa, qb, qc, J, localUxFwd, localUyFwd,
                                          localUzFwd, grad_fwd);
              computeDisplacementGradient(qa, qb, qc, J, localUxAdj, localUyAdj,
                                          localUzAdj, grad_adj);

              int const qIdx = qa + qb * dim + qc * dim * dim;

              float const div_fwd =
                  grad_fwd[0][0] + grad_fwd[1][1] + grad_fwd[2][2];
              float const div_adj =
                  grad_adj[0][0] + grad_adj[1][1] + grad_adj[2][2];
              strainDiv[qIdx] = div_adj * div_fwd;

              float const exx_f = grad_fwd[0][0];
              float const eyy_f = grad_fwd[1][1];
              float const ezz_f = grad_fwd[2][2];
              float const exy_f = 0.5f * (grad_fwd[0][1] + grad_fwd[1][0]);
              float const exz_f = 0.5f * (grad_fwd[0][2] + grad_fwd[2][0]);
              float const eyz_f = 0.5f * (grad_fwd[1][2] + grad_fwd[2][1]);

              float const exx_a = grad_adj[0][0];
              float const eyy_a = grad_adj[1][1];
              float const ezz_a = grad_adj[2][2];
              float const exy_a = 0.5f * (grad_adj[0][1] + grad_adj[1][0]);
              float const exz_a = 0.5f * (grad_adj[0][2] + grad_adj[2][0]);
              float const eyz_a = 0.5f * (grad_adj[1][2] + grad_adj[2][1]);

              strainEps[qIdx] =
                  2.0f * (exx_a * exx_f + eyy_a * eyy_f + ezz_a * ezz_f) +
                  4.0f * (exy_a * exy_f + exz_a * exz_f + eyz_a * eyz_f);
            },
            [&](int, int, float, const int, const int) {});

        INTEGRAL_TYPE::computeMassTerm(
            transformData, [&](const int q, const real_t wdetJ) {
              float const dot = localUxDt2[q] * localUxFwd[q] +
                                localUyDt2[q] * localUyFwd[q] +
                                localUzDt2[q] * localUzFwd[q];
              ATOMICADD(gradRho(localGIdx[q]), dot * wdetJ);
              ATOMICADD(gradLambda(localGIdx[q]), strainDiv[q] * wdetJ);
              ATOMICADD(gradMu(localGIdx[q]), strainEps[q] * wdetJ);
            });
      });
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_IMPL_ELASTIC_INCLUDE_DIFFERENTIATOR_ELASTIC_IMPL_H_
