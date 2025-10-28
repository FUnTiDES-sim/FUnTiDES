#include <gtest/gtest.h>

#include "common_config.h"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

/*
 * Main entry point for the Google Test framework with Kokkos initialization.
 */
int main(int argc, char **argv)
{
#ifdef USE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
#ifdef USE_KOKKOS
  Kokkos::finalize();
#endif
  return result;
}