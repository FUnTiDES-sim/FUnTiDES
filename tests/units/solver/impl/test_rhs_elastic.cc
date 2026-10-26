#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "device_view_access.h"
#include "rhs_elastic.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromHost;

using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

class RhsElasticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create test data with specific sizes
    numElements = 10;
    numNodesPerElement = 8;

    termx = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "termx");
    termy = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "termy");
    termz = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "termz");
    element = allocateVector<vectorInt>(numElements, "element");
    weights = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "weights");

    termx2 = allocateArray2D<arrayReal>(20, 12, "termx2");
    termy2 = allocateArray2D<arrayReal>(20, 12, "termy2");
    termz2 = allocateArray2D<arrayReal>(20, 12, "termz2");
    element2 = allocateVector<vectorInt>(20, "element2");
    weights2 = allocateArray2D<arrayReal>(20, 12, "weights2");

    // Initialize with test values
    fillFromHost(element, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i) h(i) = i * 10;
    });
    fillFromHost(termx, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 100 + j;
    });
    fillFromHost(termy, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 200 + j;
    });
    fillFromHost(termz, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 300 + j;
    });
    fillFromHost(weights, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i + j * 0.1;
    });

    fillFromHost(element2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i) h(i) = i * 5;
    });
    fillFromHost(termx2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 50 + j;
    });
    fillFromHost(termy2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 60 + j;
    });
    fillFromHost(termz2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 70 + j;
    });
    fillFromHost(weights2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 2 + j * 0.2;
    });
  }

  size_t numElements;
  size_t numNodesPerElement;
  arrayReal termx, termy, termz;
  arrayReal termx2, termy2, termz2;
  vectorInt element, element2;
  arrayReal weights, weights2;
};

TEST_F(RhsElasticTest, Constructor) {
  RhsElastic rhs(termx, termy, termz, element, weights);

  EXPECT_EQ(rhs.m_termx.extent(0), numElements);
  EXPECT_EQ(rhs.m_termx.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_termy.extent(0), numElements);
  EXPECT_EQ(rhs.m_termy.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_termz.extent(0), numElements);
  EXPECT_EQ(rhs.m_termz.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_element.extent(0), numElements);
  EXPECT_EQ(rhs.m_weights.extent(0), numElements);
  EXPECT_EQ(rhs.m_weights.extent(1), numNodesPerElement);

  // Verify data is correctly stored
  auto h_element = toHost(rhs.m_element);
  auto h_termx = toHost(rhs.m_termx);
  auto h_termy = toHost(rhs.m_termy);
  auto h_termz = toHost(rhs.m_termz);
  auto h_weights = toHost(rhs.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_termx(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_termy(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(h_termz(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(h_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsElasticTest, CopyConstructor) {
  RhsElastic original(termx, termy, termz, element, weights);
  RhsElastic copy(original);

  // Check that copy has the same extent
  EXPECT_EQ(copy.m_termx.extent(0), original.m_termx.extent(0));
  EXPECT_EQ(copy.m_termx.extent(1), original.m_termx.extent(1));
  EXPECT_EQ(copy.m_termy.extent(0), original.m_termy.extent(0));
  EXPECT_EQ(copy.m_termy.extent(1), original.m_termy.extent(1));
  EXPECT_EQ(copy.m_termz.extent(0), original.m_termz.extent(0));
  EXPECT_EQ(copy.m_termz.extent(1), original.m_termz.extent(1));
  EXPECT_EQ(copy.m_element.extent(0), original.m_element.extent(0));
  EXPECT_EQ(copy.m_weights.extent(0), original.m_weights.extent(0));
  EXPECT_EQ(copy.m_weights.extent(1), original.m_weights.extent(1));

  // Check that copy has the same data
  auto h_copy_element = toHost(copy.m_element);
  auto h_original_element = toHost(original.m_element);
  auto h_copy_termx = toHost(copy.m_termx);
  auto h_original_termx = toHost(original.m_termx);
  auto h_copy_termy = toHost(copy.m_termy);
  auto h_original_termy = toHost(original.m_termy);
  auto h_copy_termz = toHost(copy.m_termz);
  auto h_original_termz = toHost(original.m_termz);
  auto h_copy_weights = toHost(copy.m_weights);
  auto h_original_weights = toHost(original.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_copy_element(i), h_original_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_copy_termx(i, j), h_original_termx(i, j));
      EXPECT_FLOAT_EQ(h_copy_termy(i, j), h_original_termy(i, j));
      EXPECT_FLOAT_EQ(h_copy_termz(i, j), h_original_termz(i, j));
      EXPECT_FLOAT_EQ(h_copy_weights(i, j), h_original_weights(i, j));
    }
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(original.m_termx, 999.0f, 0, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_termx, 0, 0), 999.0f);
  setAt(original.m_termy, 888.0f, 1, 1);
  EXPECT_FLOAT_EQ(getAt(copy.m_termy, 1, 1), 888.0f);
  setAt(original.m_termz, 777.0f, 2, 2);
  EXPECT_FLOAT_EQ(getAt(copy.m_termz, 2, 2), 777.0f);
  setAt(original.m_element, 666, 3);
  EXPECT_EQ(getAt(copy.m_element, 3), 666);
  setAt(original.m_weights, 555.0f, 4, 5);
  EXPECT_FLOAT_EQ(getAt(copy.m_weights, 4, 5), 555.0f);
}

TEST_F(RhsElasticTest, CopyAssignmentOperator) {
  RhsElastic rhs1(termx, termy, termz, element, weights);
  RhsElastic rhs2(termx2, termy2, termz2, element2, weights2);

  // Verify initial state
  EXPECT_EQ(rhs2.m_termx.extent(0), 20);
  EXPECT_EQ(rhs2.m_termx.extent(1), 12);
  EXPECT_EQ(rhs2.m_element.extent(0), 20);

  // Perform assignment
  rhs2 = rhs1;

  // Check that rhs2 now has the same extent as rhs1
  EXPECT_EQ(rhs2.m_termx.extent(0), numElements);
  EXPECT_EQ(rhs2.m_termx.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs2.m_termy.extent(0), numElements);
  EXPECT_EQ(rhs2.m_termy.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs2.m_termz.extent(0), numElements);
  EXPECT_EQ(rhs2.m_termz.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs2.m_element.extent(0), numElements);
  EXPECT_EQ(rhs2.m_weights.extent(0), numElements);
  EXPECT_EQ(rhs2.m_weights.extent(1), numNodesPerElement);

  // Check that data matches
  auto h_rhs2_element = toHost(rhs2.m_element);
  auto h_rhs1_element = toHost(rhs1.m_element);
  auto h_rhs2_termx = toHost(rhs2.m_termx);
  auto h_rhs1_termx = toHost(rhs1.m_termx);
  auto h_rhs2_termy = toHost(rhs2.m_termy);
  auto h_rhs1_termy = toHost(rhs1.m_termy);
  auto h_rhs2_termz = toHost(rhs2.m_termz);
  auto h_rhs1_termz = toHost(rhs1.m_termz);
  auto h_rhs2_weights = toHost(rhs2.m_weights);
  auto h_rhs1_weights = toHost(rhs1.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_rhs2_element(i), h_rhs1_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_rhs2_termx(i, j), h_rhs1_termx(i, j));
      EXPECT_FLOAT_EQ(h_rhs2_termy(i, j), h_rhs1_termy(i, j));
      EXPECT_FLOAT_EQ(h_rhs2_termz(i, j), h_rhs1_termz(i, j));
      EXPECT_FLOAT_EQ(h_rhs2_weights(i, j), h_rhs1_weights(i, j));
    }
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(rhs1.m_termx, 777.0f, 0, 0);
  EXPECT_FLOAT_EQ(getAt(rhs2.m_termx, 0, 0), 777.0f);
  setAt(rhs1.m_termy, 666.0f, 1, 1);
  EXPECT_FLOAT_EQ(getAt(rhs2.m_termy, 1, 1), 666.0f);
  setAt(rhs1.m_element, 555, 2);
  EXPECT_EQ(getAt(rhs2.m_element, 2), 555);
}

TEST_F(RhsElasticTest, CopyAssignmentSelfAssignment) {
  RhsElastic rhs(termx, termy, termz, element, weights);

  // Self-assignment should not cause issues
  rhs = rhs;

  // Verify data is unchanged
  EXPECT_EQ(rhs.m_termx.extent(0), numElements);
  EXPECT_EQ(rhs.m_termx.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_termy.extent(0), numElements);
  EXPECT_EQ(rhs.m_termz.extent(0), numElements);
  EXPECT_EQ(rhs.m_element.extent(0), numElements);

  auto h_element = toHost(rhs.m_element);
  auto h_termx = toHost(rhs.m_termx);
  auto h_termy = toHost(rhs.m_termy);
  auto h_termz = toHost(rhs.m_termz);
  auto h_weights = toHost(rhs.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_termx(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_termy(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(h_termz(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(h_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsElasticTest, GetTerm) {
  RhsElastic rhs(termx, termy, termz, element, weights);

  auto tx = rhs.getTerm(0);
  auto ty = rhs.getTerm(1);
  auto tz = rhs.getTerm(2);

  EXPECT_EQ(tx.extent(0), numElements);
  EXPECT_EQ(tx.extent(1), numNodesPerElement);
  EXPECT_EQ(ty.extent(0), numElements);
  EXPECT_EQ(ty.extent(1), numNodesPerElement);
  EXPECT_EQ(tz.extent(0), numElements);
  EXPECT_EQ(tz.extent(1), numNodesPerElement);

  auto h_tx = toHost(tx);
  auto h_ty = toHost(ty);
  auto h_tz = toHost(tz);
  for (size_t i = 0; i < numElements; ++i) {
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_tx(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_ty(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(h_tz(i, j), i * 300 + j);
    }
  }
}

TEST_F(RhsElasticTest, GetElement) {
  RhsElastic rhs(termx, termy, termz, element, weights);

  auto elem = rhs.getElement();

  EXPECT_EQ(elem.extent(0), numElements);
  auto h_elem = toHost(elem);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_elem(i), i * 10);
  }
}

TEST_F(RhsElasticTest, GetWeights) {
  RhsElastic rhs(termx, termy, termz, element, weights);

  auto w = rhs.getWeights();

  EXPECT_EQ(w.extent(0), numElements);
  EXPECT_EQ(w.extent(1), numNodesPerElement);
  auto h_w = toHost(w);
  for (size_t i = 0; i < numElements; ++i) {
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_w(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsElasticTest, EmptyFields) {
  auto emptyTermx = allocateArray2D<arrayReal>(0, 0, "emptyTermx");
  auto emptyTermy = allocateArray2D<arrayReal>(0, 0, "emptyTermy");
  auto emptyTermz = allocateArray2D<arrayReal>(0, 0, "emptyTermz");
  auto emptyElement = allocateVector<vectorInt>(0, "emptyElement");
  auto emptyWeights = allocateArray2D<arrayReal>(0, 0, "emptyWeights");

  RhsElastic rhs(emptyTermx, emptyTermy, emptyTermz, emptyElement, emptyWeights);

  EXPECT_EQ(rhs.m_termx.extent(0), 0);
  EXPECT_EQ(rhs.m_termy.extent(0), 0);
  EXPECT_EQ(rhs.m_termz.extent(0), 0);
  EXPECT_EQ(rhs.m_element.extent(0), 0);
  EXPECT_EQ(rhs.m_weights.extent(0), 0);
}

TEST_F(RhsElasticTest, ModifyAfterConstruction) {
  RhsElastic rhs(termx, termy, termz, element, weights);

  // Modify data through the RHS object
  setAt(rhs.m_termx, 123.456f, 5, 3);
  setAt(rhs.m_termy, 234.567f, 6, 2);
  setAt(rhs.m_termz, 345.678f, 7, 1);
  setAt(rhs.m_element, 789, 8);
  setAt(rhs.m_weights, 0.999f, 2, 1);

  // Verify modifications
  EXPECT_FLOAT_EQ(getAt(rhs.m_termx, 5, 3), 123.456f);
  EXPECT_FLOAT_EQ(getAt(rhs.m_termy, 6, 2), 234.567f);
  EXPECT_FLOAT_EQ(getAt(rhs.m_termz, 7, 1), 345.678f);
  EXPECT_EQ(getAt(rhs.m_element, 8), 789);
  EXPECT_FLOAT_EQ(getAt(rhs.m_weights, 2, 1), 0.999f);

  // Original views should reflect changes (shallow copy)
  EXPECT_FLOAT_EQ(getAt(termx, 5, 3), 123.456f);
  EXPECT_FLOAT_EQ(getAt(termy, 6, 2), 234.567f);
  EXPECT_FLOAT_EQ(getAt(termz, 7, 1), 345.678f);
  EXPECT_EQ(getAt(element, 8), 789);
  EXPECT_FLOAT_EQ(getAt(weights, 2, 1), 0.999f);
}

TEST_F(RhsElasticTest, CopyInContainerClass) {
  // Create a simple container class that stores rhs by copy
  struct RhsContainer {
    RhsElastic rhs;

    RhsContainer(const RhsElastic& r) : rhs(r) {}
  };

  RhsElastic original(termx, termy, termz, element, weights);

  // Create container with rhs copy
  RhsContainer container(original);

  // Verify container has correct initial state
  EXPECT_EQ(container.rhs.m_termx.extent(0), numElements);
  EXPECT_EQ(container.rhs.m_termy.extent(0), numElements);
  EXPECT_EQ(container.rhs.m_termz.extent(0), numElements);
  EXPECT_EQ(container.rhs.m_element.extent(0), numElements);

  // Modify original
  setAt(original.m_termx, 999.0f, 0, 0);
  setAt(original.m_termy, 888.0f, 1, 1);
  setAt(original.m_termz, 777.0f, 2, 2);
  setAt(original.m_element, 666, 3);

  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(getAt(container.rhs.m_termx, 0, 0), 999.0f);
  EXPECT_FLOAT_EQ(getAt(container.rhs.m_termy, 1, 1), 888.0f);
  EXPECT_FLOAT_EQ(getAt(container.rhs.m_termz, 2, 2), 777.0f);
  EXPECT_EQ(getAt(container.rhs.m_element, 3), 666);

  // Modify through container
  setAt(container.rhs.m_weights, 0.555f, 4, 5);

  // Original should reflect change
  EXPECT_FLOAT_EQ(getAt(original.m_weights, 4, 5), 0.555f);
}

TEST_F(RhsElasticTest, GetNumRhsComponentsReturnsThree) {
  RhsElastic rhs(termx, termy, termz, element, weights);
  EXPECT_EQ(rhs.getNumRhsComponents(), 3);
}

TEST_F(RhsElasticTest, PrintDoesNotCrash) {
  RhsElastic rhs(termx, termy, termz, element, weights);
  testing::internal::CaptureStdout();
  EXPECT_NO_THROW(rhs.print());
  testing::internal::GetCapturedStdout();
}

}  // namespace test
}  // namespace fe
}  // namespace solver
