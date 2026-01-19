#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "rhs_elastic.h"

namespace solver
{
namespace fe
{
namespace test
{

class RhsElasticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    // Create test data with specific sizes
    numElements = 10;
    numNodesPerElement = 8;

    termx = allocateArray2D<ARRAY_REAL_VIEW>(numElements, numNodesPerElement,
                                             "termx");
    termy = allocateArray2D<ARRAY_REAL_VIEW>(numElements, numNodesPerElement,
                                             "termy");
    termz = allocateArray2D<ARRAY_REAL_VIEW>(numElements, numNodesPerElement,
                                             "termz");
    element = allocateVector<VECTOR_INT_VIEW>(numElements, "element");
    weights = allocateArray2D<ARRAY_REAL_VIEW>(numElements, numNodesPerElement,
                                               "weights");

    termx2 = allocateArray2D<ARRAY_REAL_VIEW>(20, 12, "termx2");
    termy2 = allocateArray2D<ARRAY_REAL_VIEW>(20, 12, "termy2");
    termz2 = allocateArray2D<ARRAY_REAL_VIEW>(20, 12, "termz2");
    element2 = allocateVector<VECTOR_INT_VIEW>(20, "element2");
    weights2 = allocateArray2D<ARRAY_REAL_VIEW>(20, 12, "weights2");

    // Initialize with test values
    for (size_t i = 0; i < numElements; ++i)
    {
      element(i) = i * 10;
      for (size_t j = 0; j < numNodesPerElement; ++j)
      {
        termx(i, j) = i * 100 + j;
        termy(i, j) = i * 200 + j;
        termz(i, j) = i * 300 + j;
        weights(i, j) = i + j * 0.1;
      }
    }

    for (size_t i = 0; i < 20; ++i)
    {
      element2(i) = i * 5;
      for (size_t j = 0; j < 12; ++j)
      {
        termx2(i, j) = i * 50 + j;
        termy2(i, j) = i * 60 + j;
        termz2(i, j) = i * 70 + j;
        weights2(i, j) = i * 2 + j * 0.2;
      }
    }
  }

  size_t numElements;
  size_t numNodesPerElement;
  ARRAY_REAL_VIEW termx, termy, termz;
  ARRAY_REAL_VIEW termx2, termy2, termz2;
  VECTOR_INT_VIEW element, element2;
  ARRAY_REAL_VIEW weights, weights2;
};

TEST_F(RhsElasticTest, Constructor)
{
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
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(rhs.m_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(rhs.m_termx(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(rhs.m_termy(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(rhs.m_termz(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(rhs.m_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsElasticTest, CopyConstructor)
{
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
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(copy.m_element(i), original.m_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(copy.m_termx(i, j), original.m_termx(i, j));
      EXPECT_FLOAT_EQ(copy.m_termy(i, j), original.m_termy(i, j));
      EXPECT_FLOAT_EQ(copy.m_termz(i, j), original.m_termz(i, j));
      EXPECT_FLOAT_EQ(copy.m_weights(i, j), original.m_weights(i, j));
    }
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  original.m_termx(0, 0) = 999.0f;
  EXPECT_FLOAT_EQ(copy.m_termx(0, 0), 999.0f);
  original.m_termy(1, 1) = 888.0f;
  EXPECT_FLOAT_EQ(copy.m_termy(1, 1), 888.0f);
  original.m_termz(2, 2) = 777.0f;
  EXPECT_FLOAT_EQ(copy.m_termz(2, 2), 777.0f);
  original.m_element(3) = 666;
  EXPECT_EQ(copy.m_element(3), 666);
  original.m_weights(4, 5) = 555.0f;
  EXPECT_FLOAT_EQ(copy.m_weights(4, 5), 555.0f);
}

TEST_F(RhsElasticTest, CopyAssignmentOperator)
{
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
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(rhs2.m_element(i), rhs1.m_element(i));
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(rhs2.m_termx(i, j), rhs1.m_termx(i, j));
      EXPECT_FLOAT_EQ(rhs2.m_termy(i, j), rhs1.m_termy(i, j));
      EXPECT_FLOAT_EQ(rhs2.m_termz(i, j), rhs1.m_termz(i, j));
      EXPECT_FLOAT_EQ(rhs2.m_weights(i, j), rhs1.m_weights(i, j));
    }
  }

  // Modify original data to verify shallow copy behavior (Kokkos views are
  // shallow by default)
  rhs1.m_termx(0, 0) = 777.0f;
  EXPECT_FLOAT_EQ(rhs2.m_termx(0, 0), 777.0f);
  rhs1.m_termy(1, 1) = 666.0f;
  EXPECT_FLOAT_EQ(rhs2.m_termy(1, 1), 666.0f);
  rhs1.m_element(2) = 555;
  EXPECT_EQ(rhs2.m_element(2), 555);
}

TEST_F(RhsElasticTest, CopyAssignmentSelfAssignment)
{
  RhsElastic rhs(termx, termy, termz, element, weights);

  // Self-assignment should not cause issues
  rhs = rhs;

  // Verify data is unchanged
  EXPECT_EQ(rhs.m_termx.extent(0), numElements);
  EXPECT_EQ(rhs.m_termx.extent(1), numNodesPerElement);
  EXPECT_EQ(rhs.m_termy.extent(0), numElements);
  EXPECT_EQ(rhs.m_termz.extent(0), numElements);
  EXPECT_EQ(rhs.m_element.extent(0), numElements);

  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(rhs.m_element(i), i * 10);
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(rhs.m_termx(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(rhs.m_termy(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(rhs.m_termz(i, j), i * 300 + j);
      EXPECT_FLOAT_EQ(rhs.m_weights(i, j), i + j * 0.1);
    }
  }
}

TEST_F(RhsElasticTest, GetTerm)
{
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

  for (size_t i = 0; i < numElements; ++i)
  {
    for (size_t j = 0; j < numNodesPerElement; ++j)
    {
      EXPECT_FLOAT_EQ(tx(i, j), i * 100 + j);
      EXPECT_FLOAT_EQ(ty(i, j), i * 200 + j);
      EXPECT_FLOAT_EQ(tz(i, j), i * 300 + j);
    }
  }
}

TEST_F(RhsElasticTest, GetElement)
{
  RhsElastic rhs(termx, termy, termz, element, weights);

  auto elem = rhs.getElement();

  EXPECT_EQ(elem.extent(0), numElements);
  for (size_t i = 0; i < numElements; ++i)
  {
    EXPECT_EQ(elem(i), i * 10);
  }
}

TEST_F(RhsElasticTest, GetWeights)
{
  RhsElastic rhs(termx, termy, termz, element, weights);

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

TEST_F(RhsElasticTest, EmptyFields)
{
  auto emptyTermx = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "emptyTermx");
  auto emptyTermy = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "emptyTermy");
  auto emptyTermz = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "emptyTermz");
  auto emptyElement = allocateVector<VECTOR_INT_VIEW>(0, "emptyElement");
  auto emptyWeights = allocateArray2D<ARRAY_REAL_VIEW>(0, 0, "emptyWeights");

  RhsElastic rhs(emptyTermx, emptyTermy, emptyTermz, emptyElement,
                 emptyWeights);

  EXPECT_EQ(rhs.m_termx.extent(0), 0);
  EXPECT_EQ(rhs.m_termy.extent(0), 0);
  EXPECT_EQ(rhs.m_termz.extent(0), 0);
  EXPECT_EQ(rhs.m_element.extent(0), 0);
  EXPECT_EQ(rhs.m_weights.extent(0), 0);
}

TEST_F(RhsElasticTest, ModifyAfterConstruction)
{
  RhsElastic rhs(termx, termy, termz, element, weights);

  // Modify data through the RHS object
  rhs.m_termx(5, 3) = 123.456f;
  rhs.m_termy(6, 2) = 234.567f;
  rhs.m_termz(7, 1) = 345.678f;
  rhs.m_element(8) = 789;
  rhs.m_weights(2, 1) = 0.999f;

  // Verify modifications
  EXPECT_FLOAT_EQ(rhs.m_termx(5, 3), 123.456f);
  EXPECT_FLOAT_EQ(rhs.m_termy(6, 2), 234.567f);
  EXPECT_FLOAT_EQ(rhs.m_termz(7, 1), 345.678f);
  EXPECT_EQ(rhs.m_element(8), 789);
  EXPECT_FLOAT_EQ(rhs.m_weights(2, 1), 0.999f);

  // Original views should reflect changes (shallow copy)
  EXPECT_FLOAT_EQ(termx(5, 3), 123.456f);
  EXPECT_FLOAT_EQ(termy(6, 2), 234.567f);
  EXPECT_FLOAT_EQ(termz(7, 1), 345.678f);
  EXPECT_EQ(element(8), 789);
  EXPECT_FLOAT_EQ(weights(2, 1), 0.999f);
}

TEST_F(RhsElasticTest, CopyInContainerClass)
{
  // Create a simple container class that stores rhs by copy
  struct RhsContainer
  {
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
  original.m_termx(0, 0) = 999.0f;
  original.m_termy(1, 1) = 888.0f;
  original.m_termz(2, 2) = 777.0f;
  original.m_element(3) = 666;

  // Container should reflect changes (shallow copy via Kokkos views)
  EXPECT_FLOAT_EQ(container.rhs.m_termx(0, 0), 999.0f);
  EXPECT_FLOAT_EQ(container.rhs.m_termy(1, 1), 888.0f);
  EXPECT_FLOAT_EQ(container.rhs.m_termz(2, 2), 777.0f);
  EXPECT_EQ(container.rhs.m_element(3), 666);

  // Modify through container
  container.rhs.m_weights(4, 5) = 0.555f;

  // Original should reflect change
  EXPECT_FLOAT_EQ(original.m_weights(4, 5), 0.555f);
}

}  // namespace test
}  // namespace fe
}  // namespace solver
