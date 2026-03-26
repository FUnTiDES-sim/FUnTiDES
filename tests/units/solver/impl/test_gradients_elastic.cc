#include <gtest/gtest.h>

#include "common_macros.h"
#include "data_type.h"
#include "gradients_elastic.h"

namespace solver
{
namespace fe
{
namespace test
{

class GradientsElasticTest : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    // Create test vectors for gradients
    num_elements = 10;
    num_nodes = 100;

    // Allocate element-based gradient vectors
    grad_rho_elem =
        allocateVector<VECTOR_REAL_VIEW>(num_elements, "grad_rho_elem");
    grad_lambda_elem =
        allocateVector<VECTOR_REAL_VIEW>(num_elements, "grad_lambda_elem");
    grad_mu_elem =
        allocateVector<VECTOR_REAL_VIEW>(num_elements, "grad_mu_elem");

    // Allocate node-based gradient vectors
    grad_rho_node =
        allocateVector<VECTOR_REAL_VIEW>(num_nodes, "grad_rho_node");
    grad_lambda_node =
        allocateVector<VECTOR_REAL_VIEW>(num_nodes, "grad_lambda_node");
    grad_mu_node = allocateVector<VECTOR_REAL_VIEW>(num_nodes, "grad_mu_node");

    // Initialize element-based gradients
    for (size_t i = 0; i < num_elements; ++i)
    {
      grad_rho_elem(i) = i * 1.0f;
      grad_lambda_elem(i) = i * 2.0f;
      grad_mu_elem(i) = i * 3.0f;
    }

    // Initialize node-based gradients
    for (size_t i = 0; i < num_nodes; ++i)
    {
      grad_rho_node(i) = i * 0.1f;
      grad_lambda_node(i) = i * 0.2f;
      grad_mu_node(i) = i * 0.3f;
    }
  }

  size_t num_elements;
  size_t num_nodes;
  VECTOR_REAL_VIEW grad_rho_elem;
  VECTOR_REAL_VIEW grad_lambda_elem;
  VECTOR_REAL_VIEW grad_mu_elem;
  VECTOR_REAL_VIEW grad_rho_node;
  VECTOR_REAL_VIEW grad_lambda_node;
  VECTOR_REAL_VIEW grad_mu_node;
};

TEST_F(GradientsElasticTest, Constructor)
{
  GradientsElastic gradients(grad_rho_elem, grad_lambda_elem, grad_mu_elem);

  EXPECT_EQ(gradients.m_gradRho.extent(0), num_elements);
  EXPECT_EQ(gradients.m_gradLambda.extent(0), num_elements);
  EXPECT_EQ(gradients.m_gradMu.extent(0), num_elements);
}

TEST_F(GradientsElasticTest, GetNumGrads)
{
  GradientsElastic gradients(grad_rho_elem, grad_lambda_elem, grad_mu_elem);

  EXPECT_EQ(gradients.getNumGrads(), 3);
}

TEST_F(GradientsElasticTest, GetGradsNames)
{
  GradientsElastic gradients(grad_rho_elem, grad_lambda_elem, grad_mu_elem);

  const char* const* names = gradients.getGradsNames();
  EXPECT_STREQ(names[0], "gradRho");
  EXPECT_STREQ(names[1], "gradLambda");
  EXPECT_STREQ(names[2], "gradMu");
}

TEST_F(GradientsElasticTest, GetCurrentGradsElement)
{
  GradientsElastic gradients(grad_rho_elem, grad_lambda_elem, grad_mu_elem);

  auto grad_rho = gradients.getCurrentGrads(0);
  auto grad_lambda = gradients.getCurrentGrads(1);
  auto grad_mu = gradients.getCurrentGrads(2);

  // Verify element-based data access
  for (size_t i = 0; i < num_elements; ++i)
  {
    EXPECT_FLOAT_EQ(grad_rho(i), i * 1.0f);
    EXPECT_FLOAT_EQ(grad_lambda(i), i * 2.0f);
    EXPECT_FLOAT_EQ(grad_mu(i), i * 3.0f);
  }
}

TEST_F(GradientsElasticTest, GetCurrentGradsNode)
{
  GradientsElastic gradients(grad_rho_node, grad_lambda_node, grad_mu_node);

  auto grad_rho = gradients.getCurrentGrads(0);
  auto grad_lambda = gradients.getCurrentGrads(1);
  auto grad_mu = gradients.getCurrentGrads(2);

  // Verify node-based data access
  for (size_t i = 0; i < num_nodes; ++i)
  {
    EXPECT_FLOAT_EQ(grad_rho(i), i * 0.1f);
    EXPECT_FLOAT_EQ(grad_lambda(i), i * 0.2f);
    EXPECT_FLOAT_EQ(grad_mu(i), i * 0.3f);
  }
}

TEST_F(GradientsElasticTest, MixedAccess)
{
  GradientsElastic gradients(grad_rho_elem, grad_lambda_elem, grad_mu_elem);

  auto grad_rho = gradients.getCurrentGrads(0);
  auto grad_lambda = gradients.getCurrentGrads(1);

  // Verify we can access different gradients independently
  EXPECT_EQ(grad_rho.extent(0), num_elements);
  EXPECT_EQ(grad_lambda.extent(0), num_elements);
  EXPECT_NE(grad_rho(0), grad_lambda(0));
  EXPECT_FLOAT_EQ(grad_rho(5), 5.0f);
  EXPECT_FLOAT_EQ(grad_lambda(5), 10.0f);
}

TEST_F(GradientsElasticTest, OutOfBoundsAccess)
{
  GradientsElastic gradients(grad_rho_elem, grad_lambda_elem, grad_mu_elem);

  // Accessing out-of-bounds gradient index should return first gradient
  auto fallback = gradients.getCurrentGrads(999);
  EXPECT_EQ(fallback.extent(0), num_elements);
}

TEST_F(GradientsElasticTest, Print)
{
  GradientsElastic gradients(grad_rho_elem, grad_lambda_elem, grad_mu_elem);

  // Just verify print() doesn't crash
  EXPECT_NO_THROW(gradients.print());
}

}  // namespace test
}  // namespace fe
}  // namespace solver
