#pragma once
#include "common.h"

template <typename QK_BASIS>
class IndexingTest : public ::testing::Test {};

template <typename QK_BASIS>
class BMatrixTest : public ::testing::Test {};

TYPED_TEST_SUITE(IndexingTest, TestedBases);
TYPED_TEST_SUITE(BMatrixTest, TestedBases);

// ============================================================================
// INDEXING TESTS — 3D
// ============================================================================

TYPED_TEST(IndexingTest, LinearIndex3DConsistency) {
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int c = 0; c < num1d; ++c) {
    for (int b = 0; b < num1d; ++b) {
      for (int a = 0; a < num1d; ++a) {
        int idx1 = QK::linearIndex3DVal(a, b, c);
        int idx2 = QK::BasisType::TensorProduct3D::linearIndex(a, b, c);

        EXPECT_EQ(idx1, idx2) << "Linear indexing should be consistent";
      }
    }
  }
}

TYPED_TEST(IndexingTest, MeshIndexToLinearIndexBijection) {
  using QK = TypeParam;

  std::set<int> nodeIndices;
  for (int k = 0; k < 8; ++k) {
    int nodeIdx = QK::meshIndexToLinearIndex3D(k);
    nodeIndices.insert(nodeIdx);

    EXPECT_GE(nodeIdx, 0);
    EXPECT_LT(nodeIdx, QK::numNodes);
  }

  EXPECT_EQ(nodeIndices.size(), 8) << "Mesh corners should map to 8 distinct node indices";
}

TYPED_TEST(IndexingTest, LinearIndex2DConsistency) {
  using QK = TypeParam;
  constexpr int num1d = QK::num1dNodes;

  for (int b = 0; b < num1d; ++b) {
    for (int a = 0; a < num1d; ++a) {
      int idx1 = QK::linearIndex2DVal(a, b);
      int idx2 = QK::BasisType::TensorProduct2D::linearIndex(a, b);

      EXPECT_EQ(idx1, idx2) << "2D linear indexing should be consistent";
    }
  }
}

TYPED_TEST(IndexingTest, MeshIndexToLinearIndex2DBijection) {
  using QK = TypeParam;

  std::set<int> indices;
  for (int k = 0; k < 4; ++k) {
    int idx = QK::meshIndexToLinearIndex2D(k);
    indices.insert(idx);
    EXPECT_GE(idx, 0);
    EXPECT_LT(idx, QK::numNodesPerFace);
  }
  EXPECT_EQ(indices.size(), 4u) << "4 face corners must map to 4 distinct indices";
}

TYPED_TEST(IndexingTest, MeshIndexToLinearIndex2DCornerValues) {
  using QK = TypeParam;
  constexpr int n = QK::num1dNodes;
  // k=0 -> (0,0), k=1 -> (n-1,0), k=2 -> (0,n-1), k=3 -> (n-1,n-1)
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(0), 0);
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(1), n - 1);
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(2), (n - 1) * n);
  EXPECT_EQ(QK::meshIndexToLinearIndex2D(3), QK::numNodesPerFace - 1);
}

// ============================================================================
// B MATRIX TESTS
// ============================================================================

TYPED_TEST(BMatrixTest, BMatrixSymmetry_ArbitraryCube) {
  using QK = TypeParam;

  real_t X[8][3];
  createArbitraryCube<QK>(X, -3.0, 1.5, -0.5, 1.8);

  int numTestPoints = (QK::numQuadraturePoints < 5) ? QK::numQuadraturePoints : 5;
  for (int q = 0; q < numTestPoints; ++q) {
    int qa, qb, qc;
    QK::BasisType::TensorProduct3D::multiIndex(q, qa, qb, qc);

    real_t J[3][3] = {{0}};
    real_t B[6] = {0};

    QK::computeBMatrix(qa, qb, qc, X, J, B);

    for (int i = 0; i < 6; ++i) {
      EXPECT_TRUE(std::isfinite(B[i])) << "B matrix component " << i << " should be finite";
    }

    EXPECT_GT(B[0], 0.0) << "B[0] (xx component) should be positive";
    EXPECT_GT(B[1], 0.0) << "B[1] (yy component) should be positive";
    EXPECT_GT(B[2], 0.0) << "B[2] (zz component) should be positive";
  }
}
