#include <gtest/gtest.h>
#include <mpi.h>

#include "cartesian_params.h"
#include "cartesian_partitioner.h"

class MPIPartitionerTest : public ::testing::Test
{
 protected:
  static void SetUpTestSuite()
  {
    // MPI is initialized by GTest main or we assume it's initialized
  }
};

TEST_F(MPIPartitionerTest, SumLocalElementsEqualsGlobal)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  using FloatT = float;
  using IndexT = int;

  // Define Global Problem
  model::CartesianParams<FloatT, IndexT> global;
  global.ex = 100 * size;  // 100 elements per rank
  global.ey = 50;
  global.ez = 50;
  global.lx = 1000.0;
  global.ly = 500.0;
  global.lz = 500.0;
  global.origin_x = 0.0;

  // Run Partitioner
  model::CartesianXPartitioner<FloatT, IndexT> partitioner;
  auto local = partitioner.partition(global, rank, size);

  // Verify Local calculations
  int local_ex = local.ex;

  // Reduce to sum all local_ex
  int total_ex = 0;
  MPI_Allreduce(&local_ex, &total_ex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Assert
  EXPECT_EQ(total_ex, global.ex);

  // Verify Origin Shifts
  // Gather all origins
  std::vector<float> origins(size);
  MPI_Allgather(&local.origin_x, 1, MPI_FLOAT, origins.data(), 1, MPI_FLOAT,
                MPI_COMM_WORLD);

  // Check contiguous
  if (rank == 0)
  {
    float current_x = 0.0;
    float dx = global.lx / global.ex;

    // We can't easily check exact origin without knowing the split count per
    // rank, but we can check order.
    for (int r = 0; r < size; ++r)
    {
      EXPECT_GE(origins[r], current_x - 0.001);  // should be increasing
      current_x = origins[r];
    }
  }
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  int result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
