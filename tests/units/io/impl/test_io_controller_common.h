#ifndef FUNTIDES_TESTS_UNITS_IO_IMPL_TEST_IO_CONTROLLER_COMMON_H_
#define FUNTIDES_TESTS_UNITS_IO_IMPL_TEST_IO_CONTROLLER_COMMON_H_

#include <gtest/gtest.h>

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "io_controller_base.h"

namespace funtides {
namespace io {
namespace test {

/// Each test gets its own directory, wiped on teardown, so a file left over
/// from a previous run can never make a test pass or fail spuriously.
class IOControllerFixture : public ::testing::Test {
 protected:
  void SetUp() override {
    const ::testing::TestInfo* info = ::testing::UnitTest::GetInstance()->current_test_info();
    dir_ = std::filesystem::temp_directory_path() / ("funtides_io_" + std::string(info->name()));
    std::filesystem::remove_all(dir_);
    std::filesystem::create_directories(dir_);
  }

  void TearDown() override { std::filesystem::remove_all(dir_); }

  IOConfig makeConfig(std::vector<std::size_t> dims = {4, 5, 6}, std::size_t nb_receiver = 3,
                      std::size_t nb_iter = 10) const {
    IOConfig cfg;
    cfg.output_dir = dir_.string();
    cfg.prefix = "test";
    cfg.local_dims = dims;
    cfg.global_dims = dims;
    cfg.start_offsets.assign(dims.size(), 0);
    cfg.nb_receiver = nb_receiver;
    cfg.nb_iter = nb_iter;
    cfg.async_snapshots = false;
    return cfg;
  }

  std::filesystem::path dir_;
};

/// Deterministic, non-trivial values: an all-zero buffer would pass even if the
/// payload were never written at all.
static HostVectorReal makeField(std::size_t n, float seed = 0.0f) {
  HostVectorReal v("field", n);
  for (std::size_t i = 0; i < n; ++i) v(i) = seed + static_cast<float>(i) * 0.25f - 3.0f;
  return v;
}

static HostArrayReal makeArray(std::size_t n0, std::size_t n1, float seed = 0.0f) {
  HostArrayReal a("array", n0, n1);
  for (std::size_t i = 0; i < n0; ++i) {
    for (std::size_t j = 0; j < n1; ++j) {
      a(i, j) = seed + static_cast<float>(i) * 100.0f + static_cast<float>(j);
    }
  }
  return a;
}

/// EXPECT_EQ and not EXPECT_FLOAT_EQ: a binary round trip must be bit exact.
static void expectEqual(const HostVectorReal& a, const HostVectorReal& b) {
  ASSERT_EQ(a.size(), b.size());
  for (std::size_t i = 0; i < a.size(); ++i) EXPECT_EQ(a(i), b(i)) << "at index " << i;
}

static void expectEqual(const HostArrayReal& a, const HostArrayReal& b) {
  ASSERT_EQ(a.extent(0), b.extent(0));
  ASSERT_EQ(a.extent(1), b.extent(1));
  for (std::size_t i = 0; i < a.extent(0); ++i) {
    for (std::size_t j = 0; j < a.extent(1); ++j) {
      EXPECT_EQ(a(i, j), b(i, j)) << "at (" << i << ", " << j << ")";
    }
  }
}

}  // namespace test
}  // namespace io
}  // namespace funtides

#endif  // FUNTIDES_TESTS_UNITS_IO_IMPL_TEST_IO_CONTROLLER_COMMON_H_
