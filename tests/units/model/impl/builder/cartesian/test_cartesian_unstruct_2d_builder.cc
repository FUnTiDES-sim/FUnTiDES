// Tests of the 2-D Cartesian builder and of ModelUnstruct with dimension 2 (COMPILE_SEM_2D).

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include "Qk_Quadrilateral_Lagrange_GaussLobatto.h"
#include "cartesian_unstruct_2d_builder.h"
#include "gllpoints.h"

namespace model {
namespace test {

using Model2D = ModelUnstruct<float, int>;

struct Mesh2D {
  int order, ex, ez;
  float lx, lz, ox, oz;
};

static std::shared_ptr<ModelApi<float, int>> build(Mesh2D const& m, bool onNodes = false, bool freeSurface = true) {
  CartesianParams<float, int> p(m.order, m.ex, 1, m.ez, m.lx, 1.0f, m.lz, onNodes, false);
  p.origin_x = m.ox;
  p.origin_z = m.oz;
  return CartesianUnstruct2DBuilder<float, int>(p).getModel(freeSurface);
}

static Mesh2D const kMeshes[] = {{1, 3, 2, 300.f, 100.f, 0.f, 0.f},
                                 {2, 4, 3, 1000.f, 450.f, -50.f, 20.f},
                                 {3, 2, 5, 80.f, 200.f, 0.f, -200.f},
                                 {4, 5, 4, 1500.f, 1200.f, 10.f, 0.f}};

TEST(CartesianUnstruct2DBuilderTest, CountsAndDimension) {
  for (auto const& m : kMeshes) {
    auto model = build(m);
    int const n1d = m.order + 1;
    EXPECT_EQ(model->getDimension(), 2);
    EXPECT_EQ(model->getOrder(), m.order);
    EXPECT_EQ(model->getNumberOfElements(), m.ex * m.ez);
    EXPECT_EQ(model->getNumberOfNodes(), (m.ex * m.order + 1) * (m.ez * m.order + 1));
    EXPECT_EQ(model->getNumberOfPointsPerElement(), n1d * n1d);
    EXPECT_FLOAT_EQ(model->domainSize(0), m.lx);
    EXPECT_FLOAT_EQ(model->domainSize(1), 0.f);
    EXPECT_FLOAT_EQ(model->domainSize(2), m.lz);
  }
}

TEST(CartesianUnstruct2DBuilderTest, NodeCoordinatesFollowElementGeometry) {
  for (auto const& m : kMeshes) {
    auto model = build(m);
    float const hx = m.lx / m.ex, hz = m.lz / m.ez;
    float const tol = 1e-4f * std::max(m.lx, m.lz);
    for (int k = 0; k < m.ez; ++k)
      for (int i = 0; i < m.ex; ++i) {
        int const e = i + k * m.ex;
        for (int c = 0; c <= m.order; ++c)
          for (int a = 0; a <= m.order; ++a) {
            int const n = model->globalNodeIndex(e, a, 0, c);
            EXPECT_EQ(model->globalNodeIndex(e, a, 3, c), n) << "j must be ignored in 2-D";
            float const xa = m.ox + i * hx + (GLLPoints::get(m.order, a) + 1) * hx / 2;
            float const zc = m.oz + k * hz + (GLLPoints::get(m.order, c) + 1) * hz / 2;
            EXPECT_NEAR(model->nodeCoord(n, 0), xa, tol);
            EXPECT_NEAR(model->nodeCoord(n, 1), 0.f, tol);
            EXPECT_NEAR(model->nodeCoord(n, 2), zc, tol);
          }
      }
  }
}

TEST(CartesianUnstruct2DBuilderTest, NeighbourElementsShareEdgeNodes) {
  for (auto const& m : kMeshes) {
    auto model = build(m);
    int const o = m.order;
    for (int k = 0; k < m.ez; ++k)
      for (int i = 0; i < m.ex; ++i) {
        int const e = i + k * m.ex;
        for (int t = 0; t <= o; ++t) {
          if (i + 1 < m.ex) EXPECT_EQ(model->globalNodeIndex(e, o, 0, t), model->globalNodeIndex(e + 1, 0, 0, t));
          if (k + 1 < m.ez) EXPECT_EQ(model->globalNodeIndex(e, t, 0, o), model->globalNodeIndex(e + m.ex, t, 0, 0));
        }
      }
  }
}

template <int ORDER>
static void checkElementAreas(Mesh2D const& m) {
  using QK = typename Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<ORDER>::type;
  auto model = build(m);
  float const area = (m.lx / m.ex) * (m.lz / m.ez);
  for (int e = 0; e < model->getNumberOfElements(); ++e) {
    real_t X[4][2];
    for (int v = 0; v < 4; ++v) {
      int const n = model->globalNodeIndex(e, (v % 2) * ORDER, 0, (v / 2) * ORDER);
      X[v][0] = model->nodeCoord(n, 0);
      X[v][1] = model->nodeCoord(n, 2);
    }
    real_t total = 0;
    bool positive = true;
    QK::computeMassTerm(X, [&](int, real_t val) { total += val; });
    for (int qb = 0; qb <= ORDER; ++qb)
      for (int qa = 0; qa <= ORDER; ++qa) {
        real_t J[2][2];
        QK::jacobianTransformation(qa, qb, X, J);
        positive = positive && QK::determinant2x2(J) > 0;
      }
    EXPECT_TRUE(positive) << "element " << e << " must be counter-clockwise (det J > 0)";
    EXPECT_NEAR(total, area, 1e-4 * area);
  }
}

TEST(CartesianUnstruct2DBuilderTest, ElementsArePositivelyOrientedWithExactArea) {
  checkElementAreas<1>(kMeshes[0]);
  checkElementAreas<2>(kMeshes[1]);
  checkElementAreas<3>(kMeshes[2]);
  checkElementAreas<4>(kMeshes[3]);
}

TEST(CartesianUnstruct2DBuilderTest, BoundaryFlags) {
  for (bool freeSurface : {true, false})
    for (auto const& m : kMeshes) {
      auto model = build(m, false, freeSurface);
      float const tol = 1e-3f * std::min(m.lx / m.ex, m.lz / m.ez);
      int nBoundary = 0;
      for (int n = 0; n < model->getNumberOfNodes(); ++n) {
        float const x = model->nodeCoord(n, 0), z = model->nodeCoord(n, 2);
        bool const left = std::fabs(x - m.ox) < tol, right = std::fabs(x - (m.ox + m.lx)) < tol;
        bool const bottom = std::fabs(z - m.oz) < tol, top = std::fabs(z - (m.oz + m.lz)) < tol;
        BoundaryFlag expected = InteriorNode;
        if (top && freeSurface)
          expected = Surface;
        else if (left || right || bottom || top)
          expected = Damping;
        EXPECT_EQ(model->boundaryType(n), expected) << "node " << n << " at (" << x << ", " << z << ")";
        nBoundary += (expected != InteriorNode);
      }
      int const nx = m.ex * m.order + 1, nz = m.ez * m.order + 1;
      EXPECT_EQ(nBoundary, 2 * nx + 2 * nz - 4);
    }
}

TEST(CartesianUnstruct2DBuilderTest, MinSpacingIsSmallestGllGap) {
  for (auto const& m : kMeshes) {
    auto model = build(m);
    float const h = std::min(m.lx / m.ex, m.lz / m.ez);
    float const gap = (GLLPoints::get(m.order, 1) - GLLPoints::get(m.order, 0)) * h / 2;
    EXPECT_NEAR(model->getMinSpacing(), gap, 1e-4f * h);
  }
}

TEST(CartesianUnstruct2DBuilderTest, HomogeneousMaterial) {
  for (bool onNodes : {false, true}) {
    auto model = build(kMeshes[1], onNodes);
    if (onNodes) {
      for (int n = 0; n < model->getNumberOfNodes(); ++n) {
        EXPECT_FLOAT_EQ(model->getModelVpOnNodes(n), 1500.f);
        EXPECT_FLOAT_EQ(model->getModelRhoOnNodes(n), 1.f);
      }
    } else {
      for (int e = 0; e < model->getNumberOfElements(); ++e) {
        EXPECT_FLOAT_EQ(model->getModelVpOnElement(e), 1500.f);
        EXPECT_FLOAT_EQ(model->getModelRhoOnElement(e), 1.f);
      }
    }
    EXPECT_FLOAT_EQ(model->getMaxSpeed(), 1500.f);
  }
}

TEST(CartesianUnstruct2DBuilderTest, RejectsUnsupportedInput) {
  CartesianParams<float, int> p(2, 4, 1, 4, 100.f, 1.f, 100.f, false, false);
  p.isElastic = true;
  EXPECT_THROW((CartesianUnstruct2DBuilder<float, int>(p)), std::invalid_argument);
  p.isElastic = false;
  p.order = 0;
  EXPECT_THROW((CartesianUnstruct2DBuilder<float, int>(p)), std::invalid_argument);
  p.order = 2;
  p.ez = 0;
  EXPECT_THROW((CartesianUnstruct2DBuilder<float, int>(p)), std::invalid_argument);
}

TEST(CartesianUnstruct2DBuilderTest, NoFaceConnectivityIn2D) {
  auto model = build(kMeshes[0]);
  EXPECT_THROW(model->buildFaceConnectivity(), std::logic_error);
}

}  // namespace test
}  // namespace model
