#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "device_view_access.h"
#include "rhs_acoustic.h"

namespace solver {
namespace fe {
namespace test {

using test_utils::fillFromHost;
using test_utils::getAt;
using test_utils::setAt;
using test_utils::toHost;

class RhsAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create test data with specific sizes
    numElements = 10;
    numNodesPerElement = 8;

    term = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "term");
    element = allocateVector<vectorInt>(numElements, "element");
    weights = allocateArray2D<arrayReal>(numElements, numNodesPerElement, "weights");

    term2 = allocateArray2D<arrayReal>(20, 12, "term2");
    element2 = allocateVector<vectorInt>(20, "element2");
    weights2 = allocateArray2D<arrayReal>(20, 12, "weights2");

    // Initialize with test values
    fillFromHost(element, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i) h(i) = i * 10;
    });
    fillFromHost(term, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i * 100 + j;
    });
    fillFromHost(weights, [&](const auto& h) {
      for (size_t i = 0; i < numElements; ++i)
        for (size_t j = 0; j < numNodesPerElement; ++j) h(i, j) = i + j * 0.1;
    });

    fillFromHost(element2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i) h(i) = i * 5;
    });
    fillFromHost(term2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 50 + j;
    });
    fillFromHost(weights2, [&](const auto& h) {
      for (size_t i = 0; i < 20; ++i)
        for (size_t j = 0; j < 12; ++j) h(i, j) = i * 2 + j * 0.2;
    });
  }

  size_t numElements;
  size_t numNodesPerElement;
  arrayReal term, term2;
  vectorInt element, element2;
  arrayReal weights, weights2;
};

TEST_F(RhsAcousticTest, Constructor) {
  RhsAcoustic rhs(term, element, weights);

  EXPECT_EQ(rhs.m_term.extent(0), numElements);
  EXPECT_EQ(rhs.m_term.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_element.extent(0), numElements);
  EXPECT_EQ(rhs.m_weights.extent(0), numElements);
  EXPECT_EQ(rhs.m_weights.extent(1), numNodesPerElement);

  // Verify data is correctly stored
  auto h_element = toHost(rhs.m_element);
  auto h_term = toHost(rhs.m_term);
  auto h_weights = toHost(rhs.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsAcousticTest, CopyConstructor) {
  RhsAcoustic original(term, element, weights);
  RhsAcoustic copy(original);

  // Check that copy has the same extent
  EXPECT_EQ(copy.m_term.extent(0), original.m_term.extent(0));
  EXPECT_EQ(copy.m_term.extent(1), original.m_term.extent(1));
  EXPECT_EQ(copy.m_element.extent(0), original.m_element.extent(0));
  EXPECT_EQ(copy.m_weights.extent(0), original.m_weights.extent(0));
  EXPECT_EQ(copy.m_weights.extent(1), original.m_weights.extent(1));

  // Check that copy has the same data
  auto h_copy_element = toHost(copy.m_element);
  auto h_original_element = toHost(original.m_element);
  auto h_copy_term = toHost(copy.m_term);
  auto h_original_term = toHost(original.m_term);
  auto h_copy_weights = toHost(copy.m_weights);
  auto h_original_weights = toHost(original.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_copy_element(i), h_original_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_copy_term(i, j), h_original_term(i, j));
      EXPECT_FLOAT_EQ(h_copy_weights(i, j), h_original_weights(i, j));
    }
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(original.m_term, 999.0f, 0, 0);
  EXPECT_FLOAT_EQ(getAt(copy.m_term, 0, 0), 999.0f);
  setAt(original.m_element, 888, 1);
  EXPECT_EQ(getAt(copy.m_element, 1), 888);
  setAt(original.m_weights, 777.0f, 2, 3);
  EXPECT_FLOAT_EQ(getAt(copy.m_weights, 2, 3), 777.0f);
}

TEST_F(RhsAcousticTest, CopyAssignmentOperator) {
  RhsAcoustic rhs1(term, element, weights);
  RhsAcoustic rhs2(term2, element2, weights2);

  // Verify initial state
  EXPECT_EQ(rhs2.m_term.extent(0), 20);
  EXPECT_EQ(rhs2.m_term.extent(1), 12);
  EXPECT_EQ(rhs2.m_element.extent(0), 20);

  // Perform assignment
  rhs2 = rhs1;

  // Check that rhs2 now has the same extent as rhs1
  EXPECT_EQ(rhs2.m_term.extent(0), numElements);
  EXPECT_EQ(rhs2.m_term.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs2.m_element.extent(0), numElements);
  EXPECT_EQ(rhs2.m_weights.extent(0), numElements);
  EXPECT_EQ(rhs2.m_weights.extent(1), numNodesPerElement);

  // Check that data matches
  auto h_rhs2_element = toHost(rhs2.m_element);
  auto h_rhs1_element = toHost(rhs1.m_element);
  auto h_rhs2_term = toHost(rhs2.m_term);
  auto h_rhs1_term = toHost(rhs1.m_term);
  auto h_rhs2_weights = toHost(rhs2.m_weights);
  auto h_rhs1_weights = toHost(rhs1.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_rhs2_element(i), h_rhs1_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_rhs2_term(i, j), h_rhs1_term(i, j));
      EXPECT_FLOAT_EQ(h_rhs2_weights(i, j), h_rhs1_weights(i, j));
    }
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  setAt(rhs1.m_term, 777.0f, 0, 0);
  EXPECT_FLOAT_EQ(getAt(rhs2.m_term, 0, 0), 777.0f);
  setAt(rhs1.m_element, 666, 1);
  EXPECT_EQ(getAt(rhs2.m_element, 1), 666);
}

TEST_F(RhsAcousticTest, CopyAssignmentSelfAssignment) {
  RhsAcoustic rhs(term, element, weights);

  // Self-assignment should not cause issues
  rhs = rhs;

  // Verify data is unchanged
  EXPECT_EQ(rhs.m_term.extent(0), numElements);
  EXPECT_EQ(rhs.m_term.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_element.extent(0), numElements);

  auto h_element = toHost(rhs.m_element);
  auto h_term = toHost(rhs.m_term);
  auto h_weights = toHost(rhs.m_weights);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsAcousticTest, GetTerm) {
  RhsAcoustic rhs(term, element, weights);

  // For acoustic, getTerm should return the same term regardless of index
  auto term0 = rhs.getTerm(0);
  auto term1 = rhs.getTerm(1);
  auto term2 = rhs.getTerm(2);

  EXPECT_EQ(term0.extent(0), numElements);
  EXPECT_EQ(term0.extent(1), numNodesPerElement);
  EXPECT_EQ(term1.extent(0), numElements);
  EXPECT_EQ(term2.extent(0), numElements);

  auto h_term0 = toHost(term0);
  auto h_term1 = toHost(term1);
  auto h_term2 = toHost(term2);
  for (size_t i = 0; i < numElements; ++i) {
    for (size_t j = 0; j < numNodesPerElement; ++j) {
      EXPECT_FLOAT_EQ(h_term0(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_term1(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(h_term2(i, j), i * 100 + j);
    }
  }
}

TEST_F(RhsAcousticTest, GetElement) {
  RhsAcoustic rhs(term, element, weights);

  auto elem = rhs.getElement();

  EXPECT_EQ(elem.extent(0), numElements);
  auto h_elem = toHost(elem);
  for (size_t i = 0; i < numElements; ++i) {
    EXPECT_EQ(h_elem(i), i * 10);
  }
}

TEST_F(RhsAcousticTest, GetWeights) {
  RhsAcoustic rhs(term, element, weights);

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

TEST_F(RhsAcousticTest, EmptyFields) {
  auto emptyTerm = allocateArray2D<arrayReal>(0, 0, "emptyTerm");
  auto emptyElement = allocateVector<vectorInt>(0, "emptyElement");
  auto emptyWeights = allocateArray2D<arrayReal>(0, 0, "emptyWeights");

  RhsAcoustic rhs(emptyTerm, emptyElement, emptyWeights);

  EXPECT_EQ(rhs.m_term.extent(0), 0);
  EXPECT_EQ(rhs.m_element.extent(0), 0);
  EXPECT_EQ(rhs.m_weights.extent(0), 0);
}

TEST_F(RhsAcousticTest, ModifyAfterConstruction) {
  RhsAcoustic rhs(term, element, weights);

  // Modify data through the RHS object
  setAt(rhs.m_term, 123.456f, 5, 3);
  setAt(rhs.m_element, 789, 7);
  setAt(rhs.m_weights, 0.999f, 2, 1);

  // Verify modifications
  EXPECT_FLOAT_EQ(getAt(rhs.m_term, 5, 3), 123.456f);
  EXPECT_EQ(getAt(rhs.m_element, 7), 789);
  EXPECT_FLOAT_EQ(getAt(rhs.m_weights, 2, 1), 0.999f);

  // Original views should reflect changes (shallow copy)
  EXPECT_FLOAT_EQ(getAt(term, 5, 3), 123.456f);
  EXPECT_EQ(getAt(element, 7), 789);
  EXPECT_FLOAT_EQ(getAt(weights, 2, 1), 0.999f);
}

TEST_F(RhsAcousticTest, CopyInContainerClass) {
  // Create a simple container class that stores rhs by copy
  struct RhsContainer {
    RhsAcoustic rhs;

    RhsContainer(const RhsAcoustic& r) : rhs(r) {}
  };

  RhsAcoustic original(term, element, weights);

  // Create container with rhs copy
  RhsContainer container(original);

  // Verify container has correct initial state
  EXPECT_EQ(container.rhs.m_term.extent(0), numElements);
  EXPECT_EQ(container.rhs.m_element.extent(0), numElements);

  // Modify original
  setAt(original.m_term, 999.0f, 0, 0);
  setAt(original.m_element, 888, 1);

  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(getAt(container.rhs.m_term, 0, 0), 999.0f);
  EXPECT_EQ(getAt(container.rhs.m_element, 1), 888);

  // Modify through container
  setAt(container.rhs.m_weights, 0.777f, 3, 4);

  // Original should reflect change
  EXPECT_FLOAT_EQ(getAt(original.m_weights, 3, 4), 0.777f);
}

TEST_F(RhsAcousticTest, GetNumRhsComponentsReturnsOne) {
  RhsAcoustic rhs(term, element, weights);
  EXPECT_EQ(rhs.getNumRhsComponents(), 1);
}

TEST_F(RhsAcousticTest, PrintDoesNotCrash) {
  RhsAcoustic rhs(term, element, weights);
  testing::internal::CaptureStdout();
  EXPECT_NO_THROW(rhs.print());
  testing::internal::GetCapturedStdout();
}

}  // namespace test
}  // namespace fe
}  // namespace solver
