#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "device_view_access.h"
#include "rhs_acoustoelastic.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromHost;

using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

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

    fillFromHost(element, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i) h(i) = i * 10;
    });
    fillFromHost(acoustic_term, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 100 + j;
    });
    fillFromHost(elastic_termx, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 200 + j;
    });
    fillFromHost(elastic_termy, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 300 + j;
    });
    fillFromHost(elastic_termz, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 400 + j;
    });
    fillFromHost(weights, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i + j * 0.1;
    });

    fillFromHost(element2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i) h(i) = i * 5;
    });
    fillFromHost(acoustic_term2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 50 + j;
    });
    fillFromHost(elastic_termx2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 60 + j;
    });
    fillFromHost(elastic_termy2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 70 + j;
    });
    fillFromHost(elastic_termz2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 80 + j;
    });
    fillFromHost(weights2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 2 + j * 0.2;
    });
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

  auto h_element = toHost(rhs.m_rhs_acoustic.m_element);
  auto h_term = toHost(rhs.m_rhs_acoustic.m_term);
  auto h_termx = toHost(rhs.m_rhs_elastic.m_termx);
  auto h_termy = toHost(rhs.m_rhs_elastic.m_termy);
  auto h_termz = toHost(rhs.m_rhs_elastic.m_termz);
  auto h_weights = toHost(rhs.m_rhs_acoustic.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_termx(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(h_termy(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(h_termz(i, j), i * 400 + j);
      EXPECT_FLOAT_EQ(h_weights(i, j), i + j * 0.1);
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
  auto h_t = toHost(t);
  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) EXPECT_FLOAT_EQ(h_t(i, j), i * 100 + j);
}

TEST_F(RhsAcoustoElasticTest, GetTermElastic) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  auto tx = rhs.getTerm(1);
  auto ty = rhs.getTerm(2);
  auto tz = rhs.getTerm(3);

  auto h_tx = toHost(tx);
  auto h_ty = toHost(ty);
  auto h_tz = toHost(tz);
  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_tx(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(h_ty(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(h_tz(i, j), i * 400 + j);
    }
}

TEST_F(RhsAcoustoElasticTest, GetElement) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  auto elem = rhs.getElement();
  EXPECT_EQ(elem.extent(0), numElements);
  auto h_elem = toHost(elem);
  for (size_t i = 0; i < numElements; ++i) EXPECT_EQ(h_elem(i), i * 10);
}

TEST_F(RhsAcoustoElasticTest, GetWeights) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  auto w = rhs.getWeights();
  EXPECT_EQ(w.extent(0), numElements);
  EXPECT_EQ(w.extent(1), numNodesPerElement);
  auto h_w = toHost(w);
  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) EXPECT_FLOAT_EQ(h_w(i, j), i + j * 0.1);
}

TEST_F(RhsAcoustoElasticTest, CopyConstructor) {
  RhsAcoustoElastic original(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);
  RhsAcoustoElastic copy(original);

  EXPECT_EQ(copy.m_rhs_acoustic.m_term.extent(0), numElements);
  EXPECT_EQ(copy.m_rhs_elastic.m_termx.extent(0), numElements);

  auto h_copy_term = toHost(copy.m_rhs_acoustic.m_term);
  auto h_original_term = toHost(original.m_rhs_acoustic.m_term);
  auto h_copy_termx = toHost(copy.m_rhs_elastic.m_termx);
  auto h_original_termx = toHost(original.m_rhs_elastic.m_termx);
  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_copy_term(i, j), h_original_term(i, j));
      EXPECT_FLOAT_EQ(h_copy_termx(i, j), h_original_termx(i, j));
    }

  // Kokkos views are shallow: modifying original affects copy
  setAt(original.m_rhs_acoustic.m_term, 999.0f, 0, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_rhs_acoustic.m_term, 0, 0), 999.0f);
  setAt(original.m_rhs_elastic.m_termx, 888.0f, 1, 2);
  EXPECT_FLOAT_EQ(getAt(copy.m_rhs_elastic.m_termx, 1, 2), 888.0f);
  setAt(original.m_rhs_elastic.m_termy, 777.0f, 3, 4);
  EXPECT_FLOAT_EQ(getAt(copy.m_rhs_elastic.m_termy, 3, 4), 777.0f);
  setAt(original.m_rhs_elastic.m_termz, 666.0f, 5, 6);
  EXPECT_FLOAT_EQ(getAt(copy.m_rhs_elastic.m_termz, 5, 6), 666.0f);
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
  setAt(rhs1.m_rhs_acoustic.m_term, 777.0f, 0, 0);
  EXPECT_FLOAT_EQ(getAt(rhs2.m_rhs_acoustic.m_term, 0, 0), 777.0f);
  setAt(rhs1.m_rhs_elastic.m_termz, 555.0f, 1, 1);
  EXPECT_FLOAT_EQ(getAt(rhs2.m_rhs_elastic.m_termz, 1, 1), 555.0f);
}

TEST_F(RhsAcoustoElasticTest, CopyAssignmentSelfAssignment) {
  RhsAcoustoElastic rhs(acoustic_term, element, weights, elastic_termx, elastic_termy, elastic_termz);

  rhs = rhs;

  EXPECT_EQ(rhs.m_rhs_acoustic.m_term.extent(0), numElements);
  EXPECT_EQ(rhs.m_rhs_elastic.m_termx.extent(0), numElements);

  auto h_term = toHost(rhs.m_rhs_acoustic.m_term);
  auto h_termx = toHost(rhs.m_rhs_elastic.m_termx);
  for (size_t i = 0; i < numElements; ++i)
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_termx(i, j), i * 200 + j);
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

  setAt(rhs.m_rhs_acoustic.m_term, 123.456f, 5, 3);
  setAt(rhs.m_rhs_elastic.m_termx, 234.567f, 2, 1);
  setAt(rhs.m_rhs_elastic.m_termy, 345.678f, 3, 2);
  setAt(rhs.m_rhs_elastic.m_termz, 456.789f, 4, 3);
  setAt(rhs.m_rhs_acoustic.m_element, 789, 7);

  EXPECT_FLOAT_EQ(getAt(rhs.m_rhs_acoustic.m_term, 5, 3), 123.456f);
  EXPECT_FLOAT_EQ(getAt(rhs.m_rhs_elastic.m_termx, 2, 1), 234.567f);
  EXPECT_FLOAT_EQ(getAt(rhs.m_rhs_elastic.m_termy, 3, 2), 345.678f);
  EXPECT_FLOAT_EQ(getAt(rhs.m_rhs_elastic.m_termz, 4, 3), 456.789f);
  EXPECT_EQ(getAt(rhs.m_rhs_acoustic.m_element, 7), 789);

  // Changes visible through original views (shallow)
  EXPECT_FLOAT_EQ(getAt(acoustic_term, 5, 3), 123.456f);
  EXPECT_FLOAT_EQ(getAt(elastic_termx, 2, 1), 234.567f);
  EXPECT_EQ(getAt(element, 7), 789);
}

}  // namespace test
}  // namespace fe
}  // namespace solver
