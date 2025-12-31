#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

#include "cartesian_params.h"
#include "cartesian_partitioner.h"

namespace model
{
namespace test
{

/// @brief Test suite for CartesianXPartitioner.
class CartesianPartitionerTest : public ::testing::Test
{
 protected:
  using FloatType = float;
  using ScalarType = int;

  /// @brief Helper to create global mesh parameters.
  CartesianParams<FloatType, ScalarType> createGlobalParams(int ex, float lx)
  {
    CartesianParams<FloatType, ScalarType> p;
    p.ex = ex;
    p.ey = 10;
    p.ez = 10;
    p.lx = lx;
    p.ly = 10.0f;
    p.lz = 10.0f;
    p.origin_x = 0.0f;
    p.origin_y = 0.0f;
    p.origin_z = 0.0f;
    p.order = 1;
    p.isModelOnNodes = true;
    p.isElastic = false;
    return p;
  }
};

TEST_F(CartesianPartitionerTest, SingleRankIsIdentity)
{
  auto global = createGlobalParams(100, 1000.0f);
  CartesianXPartitioner<FloatType, ScalarType> partitioner;

  auto local = partitioner.partition(global, 0, 1);

  EXPECT_EQ(local.ex, 100);
  EXPECT_FLOAT_EQ(local.lx, 1000.0f);
  EXPECT_FLOAT_EQ(local.origin_x, 0.0f);
  EXPECT_FLOAT_EQ(local.global_lx, 1000.0f);
}

TEST_F(CartesianPartitionerTest, EvenSplit)
{
  const int globalElements = 100;
  const int numRanks = 4;
  const float dx = 1.0f;

  auto global = createGlobalParams(globalElements, globalElements * dx);
  CartesianXPartitioner<FloatType, ScalarType> partitioner;

  int elementsPerRank = globalElements / numRanks;  // 25

  for (int rank = 0; rank < numRanks; ++rank)
  {
    auto local = partitioner.partition(global, rank, numRanks);

    EXPECT_EQ(local.ex, elementsPerRank) << "Rank " << rank;
    EXPECT_FLOAT_EQ(local.lx, elementsPerRank * dx) << "Rank " << rank;
    EXPECT_FLOAT_EQ(local.origin_x, rank * elementsPerRank * dx)
        << "Rank " << rank;
  }
}

TEST_F(CartesianPartitionerTest, UnevenSplitWithRemainder)
{
  // 10 elements, 3 ranks: 10 = 3*3 + 1
  // Remainder 1 goes to first rank: distribution [4, 3, 3]
  const int globalElements = 10;
  const int numRanks = 3;
  const float dx = 1.0f;

  auto global = createGlobalParams(globalElements, globalElements * dx);
  CartesianXPartitioner<FloatType, ScalarType> partitioner;

  std::vector<int> expectedEx = {4, 3, 3};
  std::vector<float> expectedOrigin = {0.0f, 4.0f, 7.0f};

  for (int rank = 0; rank < numRanks; ++rank)
  {
    auto local = partitioner.partition(global, rank, numRanks);

    EXPECT_EQ(local.ex, expectedEx[rank]) << "Rank " << rank;
    EXPECT_FLOAT_EQ(local.origin_x, expectedOrigin[rank]) << "Rank " << rank;
  }
}

TEST_F(CartesianPartitionerTest, PreservesNonPartitionedAxes)
{
  auto global = createGlobalParams(10, 10.0f);
  global.ey = 50;
  global.ly = 500.0f;
  global.ez = 25;
  global.lz = 250.0f;
  global.origin_y = 5.0f;
  global.origin_z = -5.0f;

  CartesianXPartitioner<FloatType, ScalarType> partitioner;
  auto local = partitioner.partition(global, 1, 2);

  // Y and Z should be unchanged
  EXPECT_EQ(local.ey, global.ey);
  EXPECT_EQ(local.ez, global.ez);
  EXPECT_FLOAT_EQ(local.ly, global.ly);
  EXPECT_FLOAT_EQ(local.lz, global.lz);
  EXPECT_FLOAT_EQ(local.origin_y, global.origin_y);
  EXPECT_FLOAT_EQ(local.origin_z, global.origin_z);
}

TEST_F(CartesianPartitionerTest, GlobalOriginOffset)
{
  auto global = createGlobalParams(10, 10.0f);  // dx = 1.0
  global.origin_x = 100.0f;                     // Global mesh starts at x=100

  CartesianXPartitioner<FloatType, ScalarType> partitioner;

  // Split into 2: Rank 0 [100, 105), Rank 1 [105, 110)
  auto r0 = partitioner.partition(global, 0, 2);
  auto r1 = partitioner.partition(global, 1, 2);

  EXPECT_FLOAT_EQ(r0.origin_x, 100.0f);
  EXPECT_FLOAT_EQ(r1.origin_x, 105.0f);

  EXPECT_FLOAT_EQ(r0.global_lx, 10.0f);
  EXPECT_FLOAT_EQ(r1.global_lx, 10.0f);
}

TEST_F(CartesianPartitionerTest, PartitionsAreContinuousAndNonOverlapping)
{
  const int globalElements = 17;  // Prime number
  auto global =
      createGlobalParams(globalElements, static_cast<float>(globalElements));
  CartesianXPartitioner<FloatType, ScalarType> partitioner;
  const int numRanks = 5;

  float minX = std::numeric_limits<float>::max();
  float maxX = std::numeric_limits<float>::lowest();
  int totalElements = 0;

  for (int rank = 0; rank < numRanks; ++rank)
  {
    auto local = partitioner.partition(global, rank, numRanks);

    float rankMinX = local.origin_x;
    float rankMaxX = local.origin_x + local.lx;

    minX = std::min(minX, rankMinX);
    maxX = std::max(maxX, rankMaxX);
    totalElements += local.ex;

    // Verify contiguousness with next rank (except last)
    if (rank < numRanks - 1)
    {
      auto next = partitioner.partition(global, rank + 1, numRanks);
      // Current Right Boundary should equal Next Left Boundary
      EXPECT_FLOAT_EQ(rankMaxX, next.origin_x);
    }
  }

  EXPECT_FLOAT_EQ(minX, global.origin_x);
  EXPECT_FLOAT_EQ(maxX, global.origin_x + global.lx);
  EXPECT_EQ(totalElements, global.ex);
}

TEST_F(CartesianPartitionerTest, FloatingPointConsistency)
{
  // Use tricky floating-point case: 13 elements at lx=13.3
  auto global = createGlobalParams(13, 13.3f);
  CartesianXPartitioner<FloatType, ScalarType> partitioner;
  const int numRanks = 5;

  float reconstructedLx = 0.0f;
  float expectedMax = global.origin_x + global.lx;

  for (int rank = 0; rank < numRanks; ++rank)
  {
    auto local = partitioner.partition(global, rank, numRanks);
    reconstructedLx += local.lx;

    if (rank == numRanks - 1)
    {
      float rankMaxX = local.origin_x + local.lx;
      EXPECT_NEAR(rankMaxX, expectedMax, 1e-5f * std::abs(expectedMax));
    }
  }

  EXPECT_NEAR(reconstructedLx, global.lx, 1e-6f * global.lx);
}

TEST_F(CartesianPartitionerTest, ThrowsOnInvalidRank)
{
  auto global = createGlobalParams(10, 10.0f);
  CartesianXPartitioner<FloatType, ScalarType> partitioner;

  EXPECT_THROW(partitioner.partition(global, -1, 4), std::invalid_argument);
  EXPECT_THROW(partitioner.partition(global, 4, 4), std::invalid_argument);
  EXPECT_THROW(partitioner.partition(global, 0, 0), std::invalid_argument);
  EXPECT_THROW(partitioner.partition(global, 0, -1), std::invalid_argument);
}

TEST_F(CartesianPartitionerTest, StressTestManyRanks)
{
  auto global = createGlobalParams(100, 100.0f);
  CartesianXPartitioner<FloatType, ScalarType> partitioner;
  const int numRanks = 1024;  // More ranks than elements

  // Verify behavior for first, middle, and last ranks
  std::vector<int> ranksToCheck = {0, 50, 99, 100, 500, 1023};

  for (int rank : ranksToCheck)
  {
    EXPECT_NO_THROW({
      auto local = partitioner.partition(global, rank, numRanks);
      EXPECT_GE(local.ex, 0);  // Should be 0 or 1
      EXPECT_GE(local.lx, 0.0f);

      // Verify origin is still correct
      // For ranks >= 100, origin should be at the end of the domain
      if (rank >= 100)
      {
        EXPECT_FLOAT_EQ(local.origin_x, 100.0f);
        EXPECT_EQ(local.ex, 0);
      }
    });
  }
}

}  // namespace test
}  // namespace model
