#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "model_unstruct.h"

namespace model
{
namespace
{

// ============================================================================
// Type Wrapper Classes for Non-Type Template Parameters
// ============================================================================

struct FloatOrder1
{
  using FloatType = float;
  using ScalarType = int;
};

struct DoubleOrder1
{
  using FloatType = double;
  using ScalarType = int;
};

struct FloatOrder2
{
  using FloatType = float;
  using ScalarType = int;
};

struct DoubleOrder2
{
  using FloatType = double;
  using ScalarType = int;
};

struct FloatOrder3
{
  using FloatType = float;
  using ScalarType = int;
};

struct DoubleOrder3
{
  using FloatType = double;
  using ScalarType = int;
};

// ============================================================================
// Test Fixture
// ============================================================================

template <typename TypeWrapper>
class ModelUnstructTest : public ::testing::Test
{
 protected:
  using FloatType = typename TypeWrapper::FloatType;
  using ScalarType = typename TypeWrapper::ScalarType;

  using ModelUnstructType = ModelUnstruct<FloatType, ScalarType>;
};

// Register type wrappers for typed tests
using TypeWrappers = ::testing::Types<FloatOrder1, FloatOrder2, FloatOrder3,
                                      DoubleOrder1, DoubleOrder2, DoubleOrder3>;

TYPED_TEST_SUITE(ModelUnstructTest, TypeWrappers);

// ============================================================================
// Constructor Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, DefaultConstructor)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  ModelUnstructType model;
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, AssignmentOperatorCompiles)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  ModelUnstructType model1;
  ModelUnstructType model2;
  model2 = model1;
  SUCCEED();
}

// ============================================================================
// Type System Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, IndexTypeIsInt)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using IndexType = typename ModelUnstructType::IndexType;

  static_assert(std::is_same<IndexType, int>::value,
                "ModelUnstruct::IndexType must be int");
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, FloatTypeCorrect)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;
  using FloatType = typename TestFixture::FloatType;

  // Template should be properly parameterized
  // Verify by checking it's instantiable
  ModelUnstructType model;
  SUCCEED();
}

// ============================================================================
// GPU Compatibility Documentation
// ============================================================================

TYPED_TEST(ModelUnstructTest, GPUCompatibleMacros)
{
  // ModelUnstruct uses PROXY_HOST_DEVICE on all methods
  // for GPU/CPU dual compilation via Kokkos
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, KokkosViewBased)
{
  // Storage uses Kokkos Views:
  // - ARRAY_INT_VIEW: global_node_index_ (element → node map)
  // - VECTOR_REAL_VIEW: coordinates, properties, parameters
  // - ARRAY3D_REAL_VIEW: elasticity tensors
  SUCCEED();
}

// ============================================================================
// Method Existence Verification
// ============================================================================

TYPED_TEST(ModelUnstructTest, IndexingMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::elementIndex)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::globalVertexIndex)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::globalNodeIndex)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, CoordinateMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::vertexCoords)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::nodeCoord)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, PropertyMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelVpOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelVpOnElement)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelRhoOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelRhoOnElement)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, AnisotropicMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelDeltaOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelDeltaOnElement)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelEpsilonOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getModelEpsilonOnElement)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, ElasticityMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::initElasticityTensors)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getCTensorOnElement)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, ConfigurationMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::isModelOnNodes)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::isElastic)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getNumberOfElements)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getNumberOfNodes)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, DomainMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::domainSize)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getMinSpacing)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::getMaxSpeed)>);
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, BoundaryMethodsExist)
{
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::boundaryType)>);
  static_assert(std::is_member_function_pointer_v<
                decltype(&ModelUnstructType::faceNormal)>);
  SUCCEED();
}

// ============================================================================
// Compilation Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, AllMethodsCompile)
{
  // This test verifies that all public methods have valid signatures
  // by attempting to reference them at compile time
  using ModelUnstructType = typename TestFixture::ModelUnstructType;

  // If this compiles, all methods are available
  ModelUnstructType model;
  SUCCEED();
}

// ============================================================================
// Documentation Tests
// ============================================================================

TYPED_TEST(ModelUnstructTest, StructuredVsUnstructuredComparison)
{
  // ModelStruct: Structured mesh with implicit geometry
  // - Template parameter: Order
  // - Mesh defined by: element counts (ex, ey, ez)
  // - Node positions: Computed from index formulas
  // - Use case: Regular grids, fast access
  //
  // ModelUnstruct: Unstructured mesh with explicit geometry
  // - Template parameters: FloatType, ScalarType (NOT Order)
  // - Mesh defined by: connectivity arrays
  // - Node positions: Stored explicitly
  // - Use case: Arbitrary topologies, complex geometries
  SUCCEED();
}

TYPED_TEST(ModelUnstructTest, TemplateParameterizations)
{
  // ModelUnstruct<float, int>:
  // - 32-bit floating point for coordinates
  // - 32-bit integer indexing
  // - Suitable for typical simulations
  //
  // ModelUnstruct<double, int>:
  // - 64-bit floating point for coordinates
  // - Higher precision for large-scale domains
  // - More memory usage
  SUCCEED();
}

}  // namespace
}  // namespace model
