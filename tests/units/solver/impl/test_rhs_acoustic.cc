#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "rhs_acoustic.h"

namespace solver
{
namespace fe
{
namespace test
{

class RhsAcousticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    // Create test data with specific sizes
    numElements = 10;
    numNodesPerElement = 8;
    
    term = allocateArray2D<ARRAY_REAL_VIEW>(numElements, numNodesPerElement, "term");
    element = allocateVector<VECTOR_INT_VIEW>(numElements, "element");
    weights = allocateArray2D<ARRAY_REAL_VIEW>(numElements, numNodesPerElement, "weights");
    
    term2 = allocateArray2D<ARRAY_REAL_VIEW>(20, 12, "term2");
    element2 = allocateVector<VECTOR_INT_VIEW>(20, "element2");
    weights2 = allocateArray2D<ARRAY_REAL_VIEW>(20, 12, "weights2");
    
    // Initialize with test values
    for (size_t i = 0; i < numElements; ++i)
    {
      element(i) = i * 10;
      for (size_t j = 0; j < numNodesPerElement; ++j)
      {
        term(i, j) = i * 100 + j;
        weights(i, j) = i + j * 0.1;
      }
    }
    
    for (size_t i = 0; i < 20; ++i)
    {
      element2(i) = i * 5;
      for (size_t j = 0; j < 12; ++j)
      {
        term2(i, j) = i * 50 + j;
        weights2(i, j) = i * 2 + j * 0.2;
      }
    }
  }

  size_t numElements;
  size_t numNodesPerElement;
  ARRAY_REAL_VIEW term, term2;
  VECTOR_INT_VIEW element, element2;
  ARRAY_REAL_VIEW weights, weights2;
};

TEST_F(RhsAcousticTest, Constructor)
{
  RhsAcoustic rhs(term, element, weights);
  
  EXPECT_EQ(rhs.m_term.extent(0), numElements);
  EXPECT_EQ(rhs.m_term.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_element.extent(0), numElements);
  EXPECT_EQ(rhs.m_weights.extent(0), numElements);
  EXPECT_EQ(rhs.m_weights.extent(1), numNodesPerElement);
  
  // Verify data is correctly stored
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(rhs.m_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(rhs.m_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(rhs.m_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsAcousticTest, CopyConstructor)
{
  RhsAcoustic original(term, element, weights);
  RhsAcoustic copy(original);
  
  // Check that copy has the same extent
  EXPECT_EQ(copy.m_term.extent(0), original.m_term.extent(0));
  EXPECT_EQ(copy.m_term.extent(1), original.m_term.extent(1));
  EXPECT_EQ(copy.m_element.extent(0), original.m_element.extent(0));
  EXPECT_EQ(copy.m_weights.extent(0), original.m_weights.extent(0));
  EXPECT_EQ(copy.m_weights.extent(1), original.m_weights.extent(1));
  
  // Check that copy has the same data
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(copy.m_element(i), original.m_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(copy.m_term(i, j), original.m_term(i, j));
      EXPECT_FLOAT_EQ(copy.m_weights(i, j), original.m_weights(i, j));
    }
  }
  
  // Modify original data to verify shallow copy behavior (Kokkos views are shallow by default)
  original.m_term(0, 0) = 999.0f;
  EXPECT_FLOAT_EQ(copy.m_term(0, 0), 999.0f);
  original.m_element(1) = 888;
  EXPECT_EQ(copy.m_element(1), 888);
  original.m_weights(2, 3) = 777.0f;
  EXPECT_FLOAT_EQ(copy.m_weights(2, 3), 777.0f);
}

TEST_F(RhsAcousticTest, CopyAssignmentOperator)
{
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
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(rhs2.m_element(i), rhs1.m_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(rhs2.m_term(i, j), rhs1.m_term(i, j));
      EXPECT_FLOAT_EQ(rhs2.m_weights(i, j), rhs1.m_weights(i, j));
    }
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are shallow by default)
  rhs1.m_term(0, 0) = 777.0f;
  EXPECT_FLOAT_EQ(rhs2.m_term(0, 0), 777.0f);
  rhs1.m_element(1) = 666;
  EXPECT_EQ(rhs2.m_element(1), 666);
}

TEST_F(RhsAcousticTest, CopyAssignmentSelfAssignment)
{
  RhsAcoustic rhs(term, element, weights);
  
  // Self-assignment should not cause issues
  rhs = rhs;
  
  // Verify data is unchanged
  EXPECT_EQ(rhs.m_term.extent(0), numElements);
  EXPECT_EQ(rhs.m_term.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_element.extent(0), numElements);
  
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(rhs.m_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(rhs.m_term(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(rhs.m_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsAcousticTest, GetTerm)
{
  RhsAcoustic rhs(term, element, weights);
  
  // For acoustic, getTerm should return the same term regardless of index
  auto term0 = rhs.getTerm(0);
  auto term1 = rhs.getTerm(1);
  auto term2 = rhs.getTerm(2);
  
  EXPECT_EQ(term0.extent(0), numElements);
  EXPECT_EQ(term0.extent(1), numNodesPerElement);
  EXPECT_EQ(term1.extent(0), numElements);
  EXPECT_EQ(term2.extent(0), numElements);
  
  for (size_t i = 0; i < numElements; ++i)
  {
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(term0(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(term1(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(term2(i, j), i * 100 + j);
    }
  }
}

TEST_F(RhsAcousticTest, GetElement)
{
  RhsAcoustic rhs(term, element, weights);
  
  auto elem = rhs.getElement();
  
  EXPECT_EQ(elem.extent(0), numElements);
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(elem(i), i * 10);
  }
}

TEST_F(RhsAcousticTest, GetWeights)
{
  RhsAcoustic rhs(term, element, weights);
  
  auto w = rhs.getWeights();
  
  EXPECT_EQ(w.extent(0), numElements);
  EXPECT_EQ(w.extent(1), numNodesPerElement);
  for (size_t i = 0; i < numElements; ++i)
  {
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(w(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsAcousticTest, EmptyFields)
{
  auto emptyTerm = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "emptyTerm");
  auto emptyElement = allocateVector<VECTOR_INT_VIEW>(0, "emptyElement");
  auto emptyWeights = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "emptyWeights");
  
  RhsAcoustic rhs(emptyTerm, emptyElement, emptyWeights);
  
  EXPECT_EQ(rhs.m_term.extent(0), 0);
  EXPECT_EQ(rhs.m_element.extent(0), 0);
  EXPECT_EQ(rhs.m_weights.extent(0), 0);
}

TEST_F(RhsAcousticTest, ModifyAfterConstruction)
{
  RhsAcoustic rhs(term, element, weights);
  
  // Modify data through the RHS object
  rhs.m_term(5, 3) = 123.456f;
  rhs.m_element(7) = 789;
  rhs.m_weights(2, 1) = 0.999f;
  
  // Verify modifications
  EXPECT_FLOAT_EQ(rhs.m_term(5, 3), 123.456f);
  EXPECT_EQ(rhs.m_element(7), 789);
  EXPECT_FLOAT_EQ(rhs.m_weights(2, 1), 0.999f);
  
  // Original views should reflect changes (shallow copy)
  EXPECT_FLOAT_EQ(term(5, 3), 123.456f);
  EXPECT_EQ(element(7), 789);
  EXPECT_FLOAT_EQ(weights(2, 1), 0.999f);
}

TEST_F(RhsAcousticTest, CopyInContainerClass)
{
  // Create a simple container class that stores rhs by copy
  struct RhsContainer
  {
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
  original.m_term(0, 0) = 999.0f;
  original.m_element(1) = 888;
  
  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(container.rhs.m_term(0, 0), 999.0f);
  EXPECT_EQ(container.rhs.m_element(1), 888);
  
  // Modify through container
  container.rhs.m_weights(3, 4) = 0.777f;
  
  // Original should reflect change
  EXPECT_FLOAT_EQ(original.m_weights(3, 4), 0.777f);
}

} // namespace test
}  // namespace fe
}  // namespace solver
