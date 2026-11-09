#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>

#include "data_type.h"
#include "gllpoints.h"
#include "gmsh_2d_reader.h"
#include "quad_observation.h"

namespace {
const std::string mesh = R"($MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
7
10 0 0 0
20 1 0 0
30 2 0 0
40 0 0 1
50 1 0 1
60 2 0 1
99 8 0 8
$EndNodes
$Elements
8
1 1 2 2 1 10 20
2 1 2 2 1 20 30
3 1 2 2 1 30 60
4 1 2 1 1 60 50
5 1 2 1 1 50 40
6 1 2 2 1 40 10
7 3 2 10 1 10 20 50 40
8 3 2 20 2 20 50 60 30
$EndElements
)";
class Gmsh2DReader : public ::testing::Test {
 protected:
  std::filesystem::path path;
  void SetUp() override {
    path = std::filesystem::temp_directory_path() /
           ("funtides_gmsh_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".msh");
  }
  void TearDown() override { std::filesystem::remove(path); }
  void write(std::string text = mesh) { std::ofstream(path) << text; }
  auto read(int order = 2) {
    return model::readGmsh2D(path, order, {{10, {1500, 1}}, {20, {3000, 2}}},
                             {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}});
  }
};

TEST_F(Gmsh2DReader, SparseTagsClockwiseQuadAndSharedGLLNodes) {
  write();
  for (int order : {1, 2, 4, 9}) {
    auto d = read(order);
    EXPECT_EQ(d.n_element_, 2);
    EXPECT_EQ(d.n_node_, (2 * order + 1) * (order + 1));
    EXPECT_EQ(d.dim_, 2);
    EXPECT_FLOAT_EQ(d.lx_, 2);
    EXPECT_FLOAT_EQ(d.lz_, 1);
    EXPECT_FLOAT_EQ(d.model_vp_element_.view_host()(0), 1500);
    EXPECT_FLOAT_EQ(d.model_vp_element_.view_host()(1), 3000);
    std::set<int> a, b;
    for (int q = 0; q < (order + 1) * (order + 1); ++q) {
      a.insert(d.global_node_index_.view_host()(0, q));
      b.insert(d.global_node_index_.view_host()(1, q));
    }
    int shared = 0;
    for (int n : a) shared += b.count(n);
    EXPECT_EQ(shared, order + 1);
    for (int j = 0; j <= order; ++j) {
      int n = d.global_node_index_.view_host()(0, order + j * (order + 1));
      EXPECT_NEAR(d.nodes_coords_x_.view_host()(n), 1, 1e-6);
      EXPECT_NEAR(d.nodes_coords_z_.view_host()(n), (GLLPoints::get(order, j) + 1) / 2, 1e-6);
    }
    for (int n = 0; n < d.n_node_; ++n)
      if (d.nodes_coords_z_.view_host()(n) == 1)
        EXPECT_EQ(d.boundaries_t_.view_host()(n), static_cast<int>(model::BoundaryFlag::Surface));
    // Construct the actual mesh too: detects inconsistent dimensions/connectivity.
    model::ModelUnstruct<float, int> actual(d);
    EXPECT_EQ(actual.getNumberOfNodes(), d.n_node_);
  }
}

TEST_F(Gmsh2DReader, RejectsUnsupportedFormatAndElementType) {
  auto text = mesh;
  text.replace(text.find("2.2 0 8"), 7, "4.1 0 8");
  write(text);
  EXPECT_THROW(read(), std::runtime_error);
  text = mesh;
  text.replace(text.find("7 3 2 10"), 8, "7 2 2 10");
  write(text);
  EXPECT_THROW(read(), std::runtime_error);
  write();
  EXPECT_THROW(read(0), std::runtime_error);
}

TEST_F(Gmsh2DReader, InterpolatesGLLFieldsWithoutSnappingReceivers) {
  write();
  model::ModelUnstruct<float, int> actual(read());
  std::vector<float> field(actual.getNumberOfNodes());
  for (int i = 0; i < actual.getNumberOfNodes(); ++i) {
    double x = actual.nodeCoord(i, 0), z = actual.nodeCoord(i, 2);
    field[i] = 1 + x * x + 2 * z * z + x * z;
  }
  for (auto point : {std::pair<double, double>{.7, .4}, {1., .33}, {1.6, .2}}) {
    auto observation = model::locateQuadObservation(actual, point.first, point.second);
    double x = point.first, z = point.second;
    EXPECT_NEAR(model::observeQuad(observation, field), 1 + x * x + 2 * z * z + x * z, 1e-6);
    double sum = 0;
    for (double w : observation.weights) sum += w;
    EXPECT_NEAR(sum, 1, 1e-12);
  }
  EXPECT_THROW(model::locateQuadObservation(actual, 2.1, .5), std::runtime_error);
}

TEST_F(Gmsh2DReader, RejectsUnknownTagsAndNonplanarGeometry) {
  auto text = mesh;
  text.replace(text.find("7 3 2 10"), 8, "7 3 2 77");
  write(text);
  EXPECT_THROW(read(), std::runtime_error);
  text = mesh;
  text.replace(text.find("10 0 0 0"), 8, "10 0 1 0");
  write(text);
  EXPECT_THROW(read(), std::runtime_error);
  text = mesh;
  text.replace(text.find("10 20 50 40"), 11, "10 20 999 40");
  write(text);
  EXPECT_THROW(read(), std::runtime_error);
}

TEST_F(Gmsh2DReader, RejectsFoldedElementsAndUnmarkedBoundary) {
  auto text = mesh;
  text.replace(text.find("10 20 50 40"), 11, "10 50 20 40");
  write(text);
  EXPECT_THROW(read(), std::runtime_error);
  text = mesh;
  text.replace(text.find("1 1 2 2 1 10 20"), 15, "1 15 2 2 1 10");
  write(text);
  EXPECT_THROW(read(), std::runtime_error);
}
}  // namespace
