#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "gradients_acoustic.h"

namespace solver {
namespace fe {
namespace test {

class GradientsAcousticTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create test vectors for gradients
    num_elements = 10;
    num_nodes = 100;

    // Allocate gradient vectors
    grad_kappa_elem = allocateVector<vectorReal>(num_elements, "grad_kappa_elem");
    grad_buoyancy_elem = allocateVector<vectorReal>(num_elements, "grad_buoyancy_elem");
    grad_kappa_node = allocateVector<vectorReal>(num_nodes, "grad_kappa_node");
    grad_buoyancy_node = allocateVector<vectorReal>(num_nodes, "grad_buoyancy_node");

    // Initialize with test values
    for (size_t i = 0; i < num_elements; ++i) {
      grad_kappa_elem(i) = i * 1.5f;
      grad_buoyancy_elem(i) = i * 2.5f;
    }

    for (size_t i = 0; i < num_nodes; ++i) {
      grad_kappa_node(i) = i * 0.5f;
      grad_buoyancy_node(i) = i * 1.0f;
    }
  }

  size_t num_elements;
  size_t num_nodes;
  vectorReal grad_kappa_elem;
  vectorReal grad_buoyancy_elem;
  vectorReal grad_kappa_node;
  vectorReal grad_buoyancy_node;
};

TEST_F(GradientsAcousticTest, Constructor) {
  GradientsAcoustic gradients(grad_kappa_elem, grad_buoyancy_elem);

  EXPECT_EQ(gradients.m_gradKappa.extent(0), num_elements);
  EXPECT_EQ(gradients.m_gradBuoyancy.extent(0), num_elements);
}

TEST_F(GradientsAcousticTest, GetNumGrads) {
  GradientsAcoustic gradients(grad_kappa_elem, grad_buoyancy_elem);

  EXPECT_EQ(gradients.getNumGrads(), 2);
}

TEST_F(GradientsAcousticTest, GetGradsNames) {
  GradientsAcoustic gradients(grad_kappa_elem, grad_buoyancy_elem);

  const char* const* names = gradients.getGradsNames();
  EXPECT_STREQ(names[0], "gradKappa");
  EXPECT_STREQ(names[1], "gradBuoyancy");
}

TEST_F(GradientsAcousticTest, GetCurrentGrads) {
  GradientsAcoustic gradients(grad_kappa_elem, grad_buoyancy_elem);

  auto grad_kappa = gradients.getCurrentGrads(0);
  auto grad_buoyancy = gradients.getCurrentGrads(1);

  // Verify data access
  for (size_t i = 0; i < num_elements; ++i) {
    EXPECT_FLOAT_EQ(grad_kappa(i), i * 1.5f);
    EXPECT_FLOAT_EQ(grad_buoyancy(i), i * 2.5f);
  }
}

TEST_F(GradientsAcousticTest, NodeBasedGradients) {
  GradientsAcoustic gradients(grad_kappa_node, grad_buoyancy_node);

  auto grad_kappa = gradients.getCurrentGrads(0);
  auto grad_buoyancy = gradients.getCurrentGrads(1);

  // Verify node-based access
  for (size_t i = 0; i < num_nodes; ++i) {
    EXPECT_FLOAT_EQ(grad_kappa(i), i * 0.5f);
    EXPECT_FLOAT_EQ(grad_buoyancy(i), i * 1.0f);
  }
}

TEST_F(GradientsAcousticTest, OutOfBoundsAccess) {
  GradientsAcoustic gradients(grad_kappa_elem, grad_buoyancy_elem);

  // Accessing out-of-bounds gradient index should return first gradient
  auto fallback = gradients.getCurrentGrads(999);
  EXPECT_EQ(fallback.extent(0), num_elements);
}

TEST_F(GradientsAcousticTest, AccumulationPattern) {
  // Test that gradients can be accumulated/modified
  GradientsAcoustic gradients(grad_kappa_elem, grad_buoyancy_elem);

  auto grad_kappa = gradients.getCurrentGrads(0);

  // Accumulate a value
  float original_val = grad_kappa(0);
  float add_val = 10.5f;

  // Simulate atomic add
  grad_kappa(0) = original_val + add_val;

  EXPECT_FLOAT_EQ(grad_kappa(0), original_val + add_val);
}

}  // namespace test
}  // namespace fe
}  // namespace solver
