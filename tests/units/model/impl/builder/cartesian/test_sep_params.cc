/**
 * @file test_sep_params.cpp
 * @brief Unit tests for SepParams SEP header file parser
 */

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

#include "cartesian_struct_builder.h"
#include "sep_params.h"

namespace model
{
namespace
{

using SepParamsFloat = SepParams<float, int>;
using SepParamsDouble = SepParams<double, int>;
using SepParamsSizeT = SepParams<float, size_t>;
using SepParamsDoubleSizeT = SepParams<double, size_t>;

const std::string kTestDataDir = "tests/units/data/sep/";

// =============================================================================
// Static File Tests - Velocity Model
// =============================================================================

TEST(SepParamsVelocityModel, ParsesDimensions)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getN1(), 201);
  EXPECT_EQ(params.getN2(), 401);
  EXPECT_EQ(params.getN3(), 101);
}

TEST(SepParamsVelocityModel, ParsesSpacing)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_FLOAT_EQ(params.getD1(), 10.0f);
  EXPECT_FLOAT_EQ(params.getD2(), 10.0f);
  EXPECT_FLOAT_EQ(params.getD3(), 10.0f);
}

TEST(SepParamsVelocityModel, ParsesOrigin)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_FLOAT_EQ(params.getO1(), 0.0f);
  EXPECT_FLOAT_EQ(params.getO2(), 0.0f);
  EXPECT_FLOAT_EQ(params.getO3(), 0.0f);
}

TEST(SepParamsVelocityModel, ParsesDataFormat)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getEsize(), 4);
  EXPECT_EQ(params.getDataFormat(), "native_float");
}

TEST(SepParamsVelocityModel, ParsesDataLabel)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getDataLabel(), "Velocity (m/s)");
}

TEST(SepParamsVelocityModel, ParsesDataFile)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  std::string expected = kTestDataDir + "velocity_3d.bin";
  EXPECT_EQ(params.getDataFile(), expected);
}

TEST(SepParamsVelocityModel, ParsesOrder)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getOrder(), 3);
}

TEST(SepParamsVelocityModel, ParsesModelOnNode)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_TRUE(params.isModelOnNode());
}

TEST(SepParamsVelocityModel, ParsesElastic)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_FALSE(params.isElastic());
}

TEST(SepParamsVelocityModel, CalculatesTotalElements)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getTotalElements(), 201 * 401 * 101);
}

TEST(SepParamsVelocityModel, CalculatesTotalBytes)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  size_t expected = static_cast<size_t>(201) * 401 * 101 * 4;
  EXPECT_EQ(params.getTotalBytes(), expected);
}

// =============================================================================
// Static File Tests - Density Model
// =============================================================================

TEST(SepParamsDensityModel, ParsesDimensions)
{
  SepParamsFloat params(kTestDataDir + "density_model.H");

  EXPECT_EQ(params.getN1(), 201);
  EXPECT_EQ(params.getN2(), 401);
  EXPECT_EQ(params.getN3(), 101);
}

TEST(SepParamsDensityModel, ParsesDataLabel)
{
  SepParamsFloat params(kTestDataDir + "density_model.H");

  EXPECT_EQ(params.getDataLabel(), "Density (kg/m3)");
}

TEST(SepParamsDensityModel, ParsesDataFile)
{
  SepParamsFloat params(kTestDataDir + "density_model.H");

  std::string expected = kTestDataDir + "density_3d.bin";
  EXPECT_EQ(params.getDataFile(), expected);
}

// =============================================================================
// Static File Tests - Elastic Model
// =============================================================================

TEST(SepParamsElasticModel, ParsesDimensions)
{
  SepParamsFloat params(kTestDataDir + "elastic_model.H");

  EXPECT_EQ(params.getN1(), 128);
  EXPECT_EQ(params.getN2(), 256);
  EXPECT_EQ(params.getN3(), 64);
}

TEST(SepParamsElasticModel, ParsesSpacing)
{
  SepParamsFloat params(kTestDataDir + "elastic_model.H");

  EXPECT_FLOAT_EQ(params.getD1(), 5.0f);
  EXPECT_FLOAT_EQ(params.getD2(), 5.0f);
  EXPECT_FLOAT_EQ(params.getD3(), 5.0f);
}

TEST(SepParamsElasticModel, ParsesNegativeOrigin)
{
  SepParamsFloat params(kTestDataDir + "elastic_model.H");

  EXPECT_FLOAT_EQ(params.getO1(), -320.0f);
  EXPECT_FLOAT_EQ(params.getO2(), -640.0f);
  EXPECT_FLOAT_EQ(params.getO3(), 0.0f);
}

TEST(SepParamsElasticModel, ParsesHigherOrder)
{
  SepParamsFloat params(kTestDataDir + "elastic_model.H");

  EXPECT_EQ(params.getOrder(), 6);
}

TEST(SepParamsElasticModel, ParsesModelOnNodeFalse)
{
  SepParamsFloat params(kTestDataDir + "elastic_model.H");

  EXPECT_FALSE(params.isModelOnNode());
}

TEST(SepParamsElasticModel, ParsesElasticTrue)
{
  SepParamsFloat params(kTestDataDir + "elastic_model.H");

  EXPECT_TRUE(params.isElastic());
}

TEST(SepParamsElasticModel, ParsesDataFiletype)
{
  SepParamsFloat params(kTestDataDir + "elastic_model.H");

  EXPECT_EQ(params.getDataFiletype(), "native_float");
}

// =============================================================================
// Static File Tests - XDR Velocity Model
// =============================================================================

TEST(SepParamsXdrVelocity, ParsesDimensions)
{
  SepParamsFloat params(kTestDataDir + "xdr_velocity.H");

  EXPECT_EQ(params.getN1(), 512);
  EXPECT_EQ(params.getN2(), 512);
  EXPECT_EQ(params.getN3(), 256);
}

TEST(SepParamsXdrVelocity, ParsesSpacing)
{
  SepParamsFloat params(kTestDataDir + "xdr_velocity.H");

  EXPECT_FLOAT_EQ(params.getD1(), 12.5f);
  EXPECT_FLOAT_EQ(params.getD2(), 12.5f);
  EXPECT_FLOAT_EQ(params.getD3(), 12.5f);
}

TEST(SepParamsXdrVelocity, ParsesXdrFormat)
{
  SepParamsFloat params(kTestDataDir + "xdr_velocity.H");

  EXPECT_EQ(params.getDataFormat(), "xdr_float");
  EXPECT_EQ(params.getDataFiletype(), "xdr_float");
}

TEST(SepParamsXdrVelocity, ParsesDataLabel)
{
  SepParamsFloat params(kTestDataDir + "xdr_velocity.H");

  EXPECT_EQ(params.getDataLabel(), "P-wave velocity (m/s)");
}

// =============================================================================
// Static File Tests - Minimal Header
// =============================================================================

TEST(SepParamsMinimal, ParsesDimensions)
{
  SepParamsFloat params(kTestDataDir + "minimal.H");

  EXPECT_EQ(params.getN1(), 100);
  EXPECT_EQ(params.getN2(), 100);
  EXPECT_EQ(params.getN3(), 50);
}

TEST(SepParamsMinimal, ParsesSpacing)
{
  SepParamsFloat params(kTestDataDir + "minimal.H");

  EXPECT_FLOAT_EQ(params.getD1(), 1.0f);
  EXPECT_FLOAT_EQ(params.getD2(), 1.0f);
  EXPECT_FLOAT_EQ(params.getD3(), 1.0f);
}

TEST(SepParamsMinimal, ParsesOrigin)
{
  SepParamsFloat params(kTestDataDir + "minimal.H");

  EXPECT_FLOAT_EQ(params.getO1(), 0.0f);
  EXPECT_FLOAT_EQ(params.getO2(), 0.0f);
  EXPECT_FLOAT_EQ(params.getO3(), 0.0f);
}

TEST(SepParamsMinimal, ParsesEsize)
{
  SepParamsFloat params(kTestDataDir + "minimal.H");

  EXPECT_EQ(params.getEsize(), 4);
}

TEST(SepParamsMinimal, DefaultModelOnNodeIsTrue)
{
  SepParamsFloat params(kTestDataDir + "minimal.H");

  EXPECT_TRUE(params.isModelOnNode());
}

TEST(SepParamsMinimal, DefaultElasticIsFalse)
{
  SepParamsFloat params(kTestDataDir + "minimal.H");

  EXPECT_FALSE(params.isElastic());
}

TEST(SepParamsMinimal, CalculatesTotalElements)
{
  SepParamsFloat params(kTestDataDir + "minimal.H");

  EXPECT_EQ(params.getTotalElements(), 100 * 100 * 50);
}

// =============================================================================
// Template Instantiation Tests
// =============================================================================

TEST(SepParamsTemplates, WorksWithFloatInt)
{
  SepParamsFloat params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getN1(), 201);
  EXPECT_FLOAT_EQ(params.getD1(), 10.0f);
}

TEST(SepParamsTemplates, WorksWithDoubleInt)
{
  SepParamsDouble params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getN1(), 201);
  EXPECT_DOUBLE_EQ(params.getD1(), 10.0);
}

TEST(SepParamsTemplates, WorksWithFloatSizeT)
{
  SepParamsSizeT params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getN1(), 201UL);
  EXPECT_FLOAT_EQ(params.getD1(), 10.0f);
}

TEST(SepParamsTemplates, WorksWithDoubleSizeT)
{
  SepParamsDoubleSizeT params(kTestDataDir + "velocity_model.H");

  EXPECT_EQ(params.getN1(), 201UL);
  EXPECT_DOUBLE_EQ(params.getD1(), 10.0);
}

// =============================================================================
// Error Handling Tests
// =============================================================================

TEST(SepParamsErrors, ThrowsOnMissingFile)
{
  EXPECT_THROW(SepParamsFloat params("/nonexistent/path/to/file.H"),
               std::runtime_error);
}

TEST(SepParamsErrors, ThrowsOnEmptyPath)
{
  EXPECT_THROW(SepParamsFloat params(""), std::runtime_error);
}

// =============================================================================
// Test Fixture for Dynamic File Tests (edge cases)
// =============================================================================

class SepParamsDynamicTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    temp_prefix_ = "/tmp/sep_test_" + std::to_string(getpid()) + "_";
  }

  void TearDown() override
  {
    for (const auto& file : created_files_)
    {
      std::remove(file.c_str());
    }
  }

  std::string CreateTempSepFile(const std::string& content)
  {
    std::string filename =
        temp_prefix_ + std::to_string(file_counter_++) + ".H";
    std::ofstream file(filename);
    file << content;
    file.close();
    created_files_.push_back(filename);
    return filename;
  }

  std::string temp_prefix_;
  std::vector<std::string> created_files_;
  int file_counter_ = 0;
};

TEST_F(SepParamsDynamicTest, ThrowsOnMissingDimensions)
{
  std::string content = R"(
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  EXPECT_THROW(SepParamsFloat params(filename), std::runtime_error);
}

TEST_F(SepParamsDynamicTest, ThrowsOnZeroDimensionN1)
{
  std::string content = R"(
n1=0
n2=100
n3=100
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  EXPECT_THROW(SepParamsFloat params(filename), std::runtime_error);
}

TEST_F(SepParamsDynamicTest, ThrowsOnZeroDimensionN2)
{
  std::string content = R"(
n1=100
n2=0
n3=100
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  EXPECT_THROW(SepParamsFloat params(filename), std::runtime_error);
}

TEST_F(SepParamsDynamicTest, ThrowsOnZeroDimensionN3)
{
  std::string content = R"(
n1=100
n2=100
n3=0
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  EXPECT_THROW(SepParamsFloat params(filename), std::runtime_error);
}

TEST_F(SepParamsDynamicTest, ThrowsOnPartialDimensions)
{
  std::string content = R"(
n1=100
n2=100
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  EXPECT_THROW(SepParamsFloat params(filename), std::runtime_error);
}

TEST_F(SepParamsDynamicTest, HandlesLeadingWhitespace)
{
  std::string content = R"(
    n1=100
	n2=200
  n3=50
d1=10.0
d2=10.0
d3=5.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_EQ(params.getN1(), 100);
  EXPECT_EQ(params.getN2(), 200);
  EXPECT_EQ(params.getN3(), 50);
}

TEST_F(SepParamsDynamicTest, HandlesTrailingWhitespace)
{
  std::string content =
      "n1=100   \n"
      "n2=200\t\n"
      "n3=50  \t  \n"
      "d1=10.0\n"
      "d2=10.0\n"
      "d3=5.0\n"
      "o1=0.0\n"
      "o2=0.0\n"
      "o3=0.0\n"
      "esize=4\n"
      "in=/data/test.bin\n";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_EQ(params.getN1(), 100);
  EXPECT_EQ(params.getN2(), 200);
  EXPECT_EQ(params.getN3(), 50);
}

TEST_F(SepParamsDynamicTest, SkipsCommentLines)
{
  std::string content = R"(
# This is a comment
n1=100
# Another comment
n2=200
n3=50
d1=10.0
d2=10.0
d3=5.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_EQ(params.getN1(), 100);
  EXPECT_EQ(params.getN2(), 200);
  EXPECT_EQ(params.getN3(), 50);
}

TEST_F(SepParamsDynamicTest, SkipsEmptyLines)
{
  std::string content = R"(

n1=100

n2=200

n3=50
d1=10.0
d2=10.0
d3=5.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_EQ(params.getN1(), 100);
  EXPECT_EQ(params.getN2(), 200);
  EXPECT_EQ(params.getN3(), 50);
}

TEST_F(SepParamsDynamicTest, SkipsLinesWithoutEquals)
{
  std::string content = R"(
n1=100
this line has no equals sign
n2=200
n3=50
d1=10.0
d2=10.0
d3=5.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_EQ(params.getN1(), 100);
  EXPECT_EQ(params.getN2(), 200);
}

TEST_F(SepParamsDynamicTest, ParsesModelOnNodeNumericOne)
{
  std::string content = R"(
n1=10
n2=10
n3=10
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
model_on_node=1
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_TRUE(params.isModelOnNode());
}

TEST_F(SepParamsDynamicTest, ParsesModelOnNodeNumericZero)
{
  std::string content = R"(
n1=10
n2=10
n3=10
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
model_on_node=0
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_FALSE(params.isModelOnNode());
}

TEST_F(SepParamsDynamicTest, ParsesElasticNumericOne)
{
  std::string content = R"(
n1=10
n2=10
n3=10
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
elastic=1
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_TRUE(params.isElastic());
}

TEST_F(SepParamsDynamicTest, PreservesAbsoluteDataPath)
{
  std::string content = R"(
n1=10
n2=10
n3=10
d1=1.0
d2=1.0
d3=1.0
o1=0.0
o2=0.0
o3=0.0
esize=4
in=/absolute/path/to/data.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsFloat params(filename);

  EXPECT_EQ(params.getDataFile(), "/absolute/path/to/data.bin");
}

TEST_F(SepParamsDynamicTest, ParsesDoublePrecisionValues)
{
  std::string content = R"(
n1=10
n2=20
n3=30
d1=0.123456789012345
d2=1.987654321098765
d3=3.141592653589793
o1=-100.5
o2=200.25
o3=0.001
esize=8
in=/data/test.bin
)";

  std::string filename = CreateTempSepFile(content);
  SepParamsDouble params(filename);

  EXPECT_DOUBLE_EQ(params.getD1(), 0.123456789012345);
  EXPECT_DOUBLE_EQ(params.getD2(), 1.987654321098765);
  EXPECT_DOUBLE_EQ(params.getD3(), 3.141592653589793);
  EXPECT_DOUBLE_EQ(params.getO1(), -100.5);
  EXPECT_DOUBLE_EQ(params.getO2(), 200.25);
  EXPECT_DOUBLE_EQ(params.getO3(), 0.001);
}

}  // namespace
}  // namespace model
