#include <gtest/gtest.h>
#include <unistd.h>

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include "cartesian_model_file_reader.h"

namespace {

static std::string writeTempFile(const std::string& content) {
  char path[] = "/tmp/model_reader_test_XXXXXX";
  int fd = mkstemp(path);
  close(fd);
  std::ofstream f(path);
  f << content;
  return std::string(path);
}

// ── happy paths ──────────────────────────────────────────────────────────────

TEST(CartesianModelFileReaderTest, SinglePropertyElement) {
  auto path = writeTempFile(
      "Model Vp element\n"
      "3\n"
      "1500.0\n"
      "2000.0\n"
      "2500.0\n");
  model::CartesianModelFileReader r(path);
  std::remove(path.c_str());

  EXPECT_FALSE(r.onNodes());
  EXPECT_EQ(r.count(), 3u);
  ASSERT_TRUE(r.has("Vp"));
  EXPECT_DOUBLE_EQ(r.get("Vp")[0], 1500.0);
  EXPECT_DOUBLE_EQ(r.get("Vp")[1], 2000.0);
  EXPECT_DOUBLE_EQ(r.get("Vp")[2], 2500.0);
}

TEST(CartesianModelFileReaderTest, MultipleProperties) {
  auto path = writeTempFile(
      "Model Vp element\n"
      "2\n"
      "1500.0\n"
      "3000.0\n"
      "\n"
      "Model Rho element\n"
      "2\n"
      "1000.0\n"
      "2000.0\n"
      "\n"
      "Model Vs element\n"
      "2\n"
      "0.0\n"
      "1500.0\n");
  model::CartesianModelFileReader r(path);
  std::remove(path.c_str());

  EXPECT_EQ(r.count(), 2u);
  ASSERT_TRUE(r.has("Vp"));
  ASSERT_TRUE(r.has("Rho"));
  ASSERT_TRUE(r.has("Vs"));
  EXPECT_DOUBLE_EQ(r.get("Vp")[0], 1500.0);
  EXPECT_DOUBLE_EQ(r.get("Rho")[1], 2000.0);
  EXPECT_DOUBLE_EQ(r.get("Vs")[0], 0.0);
}

TEST(CartesianModelFileReaderTest, NodeSupport) {
  auto path = writeTempFile(
      "Model Vs node\n"
      "1\n"
      "800.0\n");
  model::CartesianModelFileReader r(path);
  std::remove(path.c_str());

  EXPECT_TRUE(r.onNodes());
  EXPECT_EQ(r.count(), 1u);
  EXPECT_DOUBLE_EQ(r.get("Vs")[0], 800.0);
}

TEST(CartesianModelFileReaderTest, EmptyLinesSkipped) {
  auto path = writeTempFile(
      "\n"
      "\n"
      "Model Vp element\n"
      "\n"
      "2\n"
      "\n"
      "1.0\n"
      "2.0\n");
  model::CartesianModelFileReader r(path);
  std::remove(path.c_str());

  EXPECT_EQ(r.count(), 2u);
  EXPECT_DOUBLE_EQ(r.get("Vp")[0], 1.0);
  EXPECT_DOUBLE_EQ(r.get("Vp")[1], 2.0);
}

TEST(CartesianModelFileReaderTest, HasReturnsFalseForAbsentProperty) {
  auto path = writeTempFile(
      "Model Vp element\n"
      "1\n"
      "1500.0\n");
  model::CartesianModelFileReader r(path);
  std::remove(path.c_str());

  EXPECT_FALSE(r.has("Vs"));
  EXPECT_FALSE(r.has("Rho"));
}

TEST(CartesianModelFileReaderTest, GetThrowsForAbsentProperty) {
  auto path = writeTempFile(
      "Model Vp element\n"
      "1\n"
      "1500.0\n");
  model::CartesianModelFileReader r(path);
  std::remove(path.c_str());

  EXPECT_THROW(r.get("Vs"), std::runtime_error);
}

// ── error paths ──────────────────────────────────────────────────────────────

TEST(CartesianModelFileReaderTest, FileNotFound) {
  EXPECT_THROW({ model::CartesianModelFileReader r("/tmp/__nonexistent_funtides_model__.txt"); }, std::runtime_error);
}

TEST(CartesianModelFileReaderTest, MalformedHeaderMissingSupport) {
  auto path = writeTempFile("Model Vp\n2\n1.0\n2.0\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

TEST(CartesianModelFileReaderTest, UnknownSupport) {
  auto path = writeTempFile("Model Vp blocks\n2\n1.0\n2.0\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

TEST(CartesianModelFileReaderTest, MixedSupports) {
  auto path = writeTempFile(
      "Model Vp element\n"
      "2\n"
      "1.0\n"
      "2.0\n"
      "Model Rho node\n"
      "2\n"
      "1.0\n"
      "2.0\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

TEST(CartesianModelFileReaderTest, CountMismatch) {
  auto path = writeTempFile(
      "Model Vp element\n"
      "2\n"
      "1.0\n"
      "2.0\n"
      "Model Rho element\n"
      "3\n"
      "1.0\n"
      "2.0\n"
      "3.0\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

TEST(CartesianModelFileReaderTest, NonNumericValue) {
  auto path = writeTempFile("Model Vp element\n2\n1.0\nnotanumber\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

TEST(CartesianModelFileReaderTest, UnexpectedEof) {
  auto path = writeTempFile("Model Vp element\n3\n1.0\n2.0\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

TEST(CartesianModelFileReaderTest, DuplicateProperty) {
  auto path = writeTempFile(
      "Model Vp element\n"
      "2\n"
      "1.0\n"
      "2.0\n"
      "Model Vp element\n"
      "2\n"
      "3.0\n"
      "4.0\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

TEST(CartesianModelFileReaderTest, NoValidSection) {
  auto path = writeTempFile("This is not a model file\nsome garbage data\n");
  EXPECT_THROW({ model::CartesianModelFileReader r(path); }, std::runtime_error);
  std::remove(path.c_str());
}

}  // namespace
