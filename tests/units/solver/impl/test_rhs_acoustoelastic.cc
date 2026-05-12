#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "rhs_acoustoelastic.h"

namespace solver {
namespace fe {
namespace test {

class RhsAcoustoElasticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    numElements = 10;
    numNodesPerElement = 8;

    acoustic_term = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "acoustic_term");
    elastic_termx = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "elastic_termx");
    elastic_termy = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "elastic_termy");
    elastic_termz = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "elastic_termz");
    element = allocateVector<vectorInt>(numElements, "element");
    weights = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "weights");

    acoustic_term2 = allocateArray2D<arrayReal>(20, 12, "acoustic_term2");
    elastic_termx2 = allocateArray2D<arrayReal>(20, 12, "elastic_termx2");
    elastic_termy2 = allocateArray2D<arrayReal>(20, 12, "elastic_termy2");
    elastic_termz2 = allocateArray2D<arrayReal>(20, 12, "elastic_termz2");
    element2 = allocateVector<vectorInt>(20, "element2");
    weights2 = allocateArray2D<arrayReal>(20, 12, "weights2");

    for (size_t i = 0; i < numElements; ++i) {
      element(i) = i * 10;
      for (size_t j = 0; j < numNodesPerElement; ++j) {
        acoustic_term(i, j) = i * 100 + j;
        elastic_termx(i, j) = i * 200 + j;
        elastic_termy(i, j) = i * 300 + j;
        elastic_termz(i, j) = i * 400 + j;
        weights(i, j) = i + j * 0.1;
      }
    }

    for (size_t i = 0; i < 20; ++i) {
      element2(i) = i * 5;
      for (size_t j = 0; j < 12; ++j) {
        acoustic_term2(i, j) = i * 50 + j;
        elastic_termx2(i, j) = i * 60 + j;
        elastic_termy2(i, j) = i * 70 + j;
        elastic_termz2(i, j) = i * 80 + j;
        weights2(i, j) = i * 2 + j * 0.2;
      }
    }
  }

  size_t numElements;
  size_t numNodesPerElement;
  arrayReal acoustic_term, elastic_termx, elastic_termy, elastic_termz;
  arrayReal acoustic_term2, elastic_termx2, elastic_termy2, elastic_termz2;
  vectorInt element, element2;
  arrayReal weights, weights2;
};

TEST_F(RhsAcoustoElasticTest, Constructor) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  EXPECT_EQ(rhs.m_rhs_acoustic.m_term.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_acoustic.m_term.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_rhs_acoustic.m_element.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_acoustic.m_weights.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termx.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termy.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termz.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_elastic.m_element.extent(0), numElements);

  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(rhs.m_rhs_acoustic.m_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(rhs.m_rhs_acoustic.m_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(rhs.m_rhs_elastic.m_termx(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(rhs.m_rhs_elastic.m_termy(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(rhs.m_rhs_elastic.m_termz(i, j), i * 400 + j);
      EXPECT_FLOAT_EQ(rhs.m_rhs_acoustic.m_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsAcoustoElasticTest, GetNumRhsComponents) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);
  EXPECT_EQ(rhs.getNumRhsComponents(), 4);
}

TEST_F(RhsAcoustoElasticTest, GetTermAcoustic) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  auto t = rhs.getTerm(0);
  EXPECT_EQ(t.extent(0), numElements);
  EXPECT_EQ(t.extent(1), numNodesPerElement);
  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) EXPECT_FLOAT_EQ(t(i, j), i * 100 + j);
}

TEST_F(RhsAcoustoElasticTest, GetTermElastic) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  auto tx = rhs.getTerm(1);
  auto ty = rhs.getTerm(2);
  auto tz = rhs.getTerm(3);

  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(tx(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(ty(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(tz(i, j), i * 400 + j);
    }
}

TEST_F(RhsAcoustoElasticTest, GetElement) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  auto elem = rhs.getElement();
  EXPECT_EQ(elem.extent(0), numElements);
  for (size_t i = 0; i < numElements; ++i) EXPECT_EQ(elem(i), i * 10);
}

TEST_F(RhsAcoustoElasticTest, GetWeights) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  auto w = rhs.getWeights();
  EXPECT_EQ(w.extent(0), numElements);
  EXPECT_EQ(w.extent(1), numNodesPerElement);
  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) EXPECT_FLOAT_EQ(w(i, j), i + j * 0.1);
}

TEST_F(RhsAcoustoElasticTest, CopyConstructor) {
  RhsAcoustoElastic original(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);
  RhsAcoustoElastic copy(original);

  EXPECT_EQ(copy.m_rhs_acoustic.m_term.extent(0), numElements);
  EXPECT_EQ(copy.m_rhs_elastic.m_termx.extent(0), numElements);

  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(copy.m_rhs_acoustic.m_term(i, j), original.m_rhs_acoustic.m_term(i, j));
      EXPECT_FLOAT_EQ(copy.m_rhs_elastic.m_termx(i, j), original.m_rhs_elastic.m_termx(i, j));
    }

  // Kokkos views are shallow: modifying original affects copy
  original.m_rhs_acoustic.m_term(0, 0) = 999.0f;
  EXPECT_FLOAT_EQ(copy.m_rhs_acoustic.m_term(0, 0), 999.0f);
  original.m_rhs_elastic.m_termx(1, 2) = 888.0f;
  EXPECT_FLOAT_EQ(copy.m_rhs_elastic.m_termx(1, 2), 888.0f);
  original.m_rhs_elastic.m_termy(3, 4) = 777.0f;
  EXPECT_FLOAT_EQ(copy.m_rhs_elastic.m_termy(3, 4), 777.0f);
  original.m_rhs_elastic.m_termz(5, 6) = 666.0f;
  EXPECT_FLOAT_EQ(copy.m_rhs_elastic.m_termz(5, 6), 666.0f);
}

TEST_F(RhsAcoustoElasticTest, CopyAssignmentOperator) {
  RhsAcoustoElastic rhs1(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);
  RhsAcoustoElastic rhs2(acoustic_term2, element2, weights2, elastic_termx2, elastic_termy2, elastic_termz2);

  EXPECT_EQ(rhs2.m_rhs_acoustic.m_term.extent(0), 20);
  EXPECT_EQ(rhs2.m_rhs_elastic.m_termx.extent(0), 20);

  rhs2 = rhs1;

  EXPECT_EQ(rhs2.m_rhs_acoustic.m_term.extent(0), numElements);
  EXPECT_EQ(rhs2.m_rhs_elastic.m_termx.extent(0), numElements);
  EXPECT_EQ(rhs2.m_rhs_elastic.m_termy.extent(0), numElements);
  EXPECT_EQ(rhs2.m_rhs_elastic.m_termz.extent(0), numElements);

  // Shallow copy after assignment
  rhs1.m_rhs_acoustic.m_term(0, 0) = 777.0f;
  EXPECT_FLOAT_EQ(rhs2.m_rhs_acoustic.m_term(0, 0), 777.0f);
  rhs1.m_rhs_elastic.m_termz(1, 1) = 555.0f;
  EXPECT_FLOAT_EQ(rhs2.m_rhs_elastic.m_termz(1, 1), 555.0f);
}

TEST_F(RhsAcoustoElasticTest, CopyAssignmentSelfAssignment) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  rhs = rhs;

  EXPECT_EQ(rhs.m_rhs_acoustic.m_term.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termx.extent(0), numElements);

  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(rhs.m_rhs_acoustic.m_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(rhs.m_rhs_elastic.m_termx(i, j), i * 200 + j);
    }
}

TEST_F(RhsAcoustoElasticTest, EmptyFields) {
  auto emptyTerm = allocateArray2D<arrayReal>(0, 0, "emptyTerm");
  auto emptyElement = allocateVector<vectorInt>(0, "emptyElement");
  auto emptyWeights = allocateArray2D<arrayReal>(0, 0, "emptyWeights");

  RhsAcoustoElastic rhs(emptyTerm, emptyElement, emptyWeights, emptyTerm, emptyTerm, emptyTerm);

  EXPECT_EQ(rhs.m_rhs_acoustic.m_term.extent(0), 0);
  EXPECT_EQ(rhs.m_rhs_acoustic.m_element.extent(0), 0);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termx.extent(0), 0);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termy.extent(0), 0);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termz.extent(0), 0);
}

TEST_F(RhsAcoustoElasticTest, ModifyAfterConstruction) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  rhs.m_rhs_acoustic.m_term(5, 3) = 123.456f;
  rhs.m_rhs_elastic.m_termx(2, 1) = 234.567f;
  rhs.m_rhs_elastic.m_termy(3, 2) = 345.678f;
  rhs.m_rhs_elastic.m_termz(4, 3) = 456.789f;
  rhs.m_rhs_acoustic.m_element(7) = 789;

  EXPECT_FLOAT_EQ(rhs.m_rhs_acoustic.m_term(5, 3), 123.456f);
  EXPECT_FLOAT_EQ(rhs.m_rhs_elastic.m_termx(2, 1), 234.567f);
  EXPECT_FLOAT_EQ(rhs.m_rhs_elastic.m_termy(3, 2), 345.678f);
  EXPECT_FLOAT_EQ(rhs.m_rhs_elastic.m_termz(4, 3), 456.789f);
  EXPECT_EQ(rhs.m_rhs_acoustic.m_element(7), 789);

  // Changes visible through original views (shallow)
  EXPECT_FLOAT_EQ(acoustic_term(5, 3), 123.456f);
  EXPECT_FLOAT_EQ(elastic_termx(2, 1), 234.567f);
  EXPECT_EQ(element(7), 789);
}

}  // namespace test
}  // namespace fe
}  // namespace solver
