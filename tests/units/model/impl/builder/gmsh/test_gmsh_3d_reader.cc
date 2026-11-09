#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

#include "gmsh_3d_reader.h"
#ifdef COMPILE_SEM
#include <cmath>

#include "Integrals.h"
#include "cartesian_unstruct_builder.h"
#include "rhs_acoustic.h"
#include "sem_solver.h"
#include "sem_solver_impl.h"
#include "wavefield_acoustic.h"
#endif
namespace {
// Two cubes, sparse tags, second cell rotated about its local x axis.
std::string meshText() {
  std::ostringstream s;
  s << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n12\n";
  auto tag = [](int x, int y, int z) { return 10 * (1 + x + 3 * y + 6 * z); };
  for (int z = 0; z < 2; ++z)
    for (int y = 0; y < 2; ++y)
      for (int x = 0; x < 3; ++x) s << tag(x, y, z) << " " << x << " " << y << " " << z << "\n";
  s << "$EndNodes\n$Elements\n12\n";
  int id = 0;
  auto quad = [&](int a, int b, int c, int d, int physical) {
    s << ++id << " 3 2 " << physical << " 1 " << a << " " << b << " " << c << " " << d << "\n";
  };
  for (int x = 0; x < 2; ++x) {
    quad(tag(x, 0, 0), tag(x + 1, 0, 0), tag(x + 1, 1, 0), tag(x, 1, 0), 2);
    quad(tag(x, 0, 1), tag(x + 1, 0, 1), tag(x + 1, 1, 1), tag(x, 1, 1), 1);
    quad(tag(x, 0, 0), tag(x + 1, 0, 0), tag(x + 1, 0, 1), tag(x, 0, 1), 2);
    quad(tag(x, 1, 0), tag(x + 1, 1, 0), tag(x + 1, 1, 1), tag(x, 1, 1), 2);
  }
  for (int x : {0, 2}) quad(tag(x, 0, 0), tag(x, 1, 0), tag(x, 1, 1), tag(x, 0, 1), 2);
  s << "11 5 2 10 1 10 20 50 40 70 80 110 100\n";
  s << "12 5 2 20 2 80 90 30 20 110 120 60 50\n$EndElements\n";
  return s.str();
}
class Gmsh3DReader : public ::testing::Test {
 protected:
  std::filesystem::path path;
  void SetUp() override {
    path = std::filesystem::temp_directory_path() /
           ("gmsh3d_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".msh");
    write(meshText());
  }
  void TearDown() override { std::filesystem::remove(path); }
  void write(const std::string& s) { std::ofstream(path) << s; }
  auto read(int p = 2) {
    return model::readGmsh3D(path, p, {{10, {2, 1}}, {20, {3, 2}}},
                             {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}});
  }
};
TEST_F(Gmsh3DReader, SharedRotatedFaceAndBoundaryPriority) {
  for (int p : {1, 2, 3, 4, 9}) {
    auto d = read(p);
    EXPECT_EQ(d.dim_, 3);
    EXPECT_EQ(d.n_node_, (2 * p + 1) * (p + 1) * (p + 1));
    EXPECT_FLOAT_EQ(d.lx_, 2);
    EXPECT_FLOAT_EQ(d.ly_, 1);
    EXPECT_FLOAT_EQ(d.lz_, 1);
    EXPECT_FLOAT_EQ(d.model_vp_element_.view_host()(1), 3);
    std::set<int> a, b;
    for (int q = 0; q < (p + 1) * (p + 1) * (p + 1); ++q) {
      a.insert(d.global_node_index_.view_host()(0, q));
      b.insert(d.global_node_index_.view_host()(1, q));
    }
    int common = 0;
    for (int id : a) common += b.count(id);
    EXPECT_EQ(common, (p + 1) * (p + 1));
    for (int id = 0; id < d.n_node_; ++id) {
      float z = d.nodes_coords_z_.view_host()(id);
      if (z == 1) EXPECT_EQ(d.boundaries_t_.view_host()(id), int(model::BoundaryFlag::Surface));
      if (d.nodes_coords_x_.view_host()(id) == 1 && z > 0 && z < 1 && d.nodes_coords_y_.view_host()(id) > 0 &&
          d.nodes_coords_y_.view_host()(id) < 1)
        EXPECT_EQ(d.boundaries_t_.view_host()(id), int(model::BoundaryFlag::InteriorNode));
    }
    model::ModelUnstruct<float, int> m(d);
    m.buildFaceConnectivity();
    EXPECT_EQ(m.getNumberOfFaces(), 11);
  }
}
TEST_F(Gmsh3DReader, TrilinearWarpedHexGeometry) {
  auto text = meshText();
  auto pos = text.find("120 2 1 1");
  ASSERT_NE(pos, std::string::npos);
  text.replace(pos, 9, "120 2 1 1.25");
  write(text);
  auto d = read(2);
  int center = d.global_node_index_.view_host()(1, 13);
  EXPECT_FLOAT_EQ(d.nodes_coords_x_.view_host()(center), 1.5f);
  EXPECT_FLOAT_EQ(d.nodes_coords_y_.view_host()(center), 0.5f);
  EXPECT_FLOAT_EQ(d.nodes_coords_z_.view_host()(center), 0.53125f);
  EXPECT_EQ(d.n_node_, 45);
}
TEST_F(Gmsh3DReader, RejectsBadInput) {
  for (auto pair : {std::pair<std::string, std::string>{"2.2 0 8", "4.1 0 8"},
                    {"11 5 2 10 1 10 20 50 40 70 80 110 100", "11 5 2 10 1 20 10 40 50 80 70 100 110"},
                    {"11 5 2 10", "11 4 2 10"},
                    {"11 5 2 10", "11 5 2 99"},
                    {"1 3 2 2 1 10 20 50 40", "1 3 2 99 1 10 20 50 40"},
                    {"11 5 2 10 1 10 20", "11 5 2 10 1 10 10"}}) {
    auto s = meshText();
    auto pos = s.find(pair.first);
    ASSERT_NE(pos, std::string::npos);
    s.replace(pos, pair.first.size(), pair.second);
    write(s);
    EXPECT_THROW(read(), std::runtime_error);
  }
  write(meshText());
  EXPECT_THROW(read(0), std::runtime_error);
  auto s = meshText();
  auto start = s.find("1 3 2");
  auto end = s.find('\n', start);
  s.erase(start, end - start + 1);
  s.replace(s.find("$Elements\n12"), 12, "$Elements\n11");
  write(s);
  EXPECT_THROW(read(), std::runtime_error);
}
#ifdef COMPILE_SEM
TEST_F(Gmsh3DReader, NonconstantFieldMatchesCartesianSEM3D) {
  using Integral = typename IntegralTypeSelector<2, IntegralType::MAKUTU>::type;
  using Mesh = model::ModelUnstruct<float, int>;
  using Solver = solver::SEMsolver<2, Integral, Mesh, false, utils::enums::physicType::kAcoustic>;
  auto d = model::readGmsh3D(path, 2, {{10, {1500, 1}}, {20, {1500, 1}}},
                             {{1, model::BoundaryFlag::Damping}, {2, model::BoundaryFlag::Damping}});
  auto imported = std::make_shared<Mesh>(d);
  imported->buildFaceConnectivity();
  model::CartesianParams<float, int> params(2, 2, 1, 1, 2, 1, 1, false, false);
  auto reference = model::CartesianUnstructBuilder<float, int>(params).getModel(false);
  std::vector<std::array<float, 5>> results[2];
  int which = 0;
  for (auto mesh : {std::shared_ptr<model::ModelApi<float, int>>(imported), reference}) {
    Solver sem;
    sem.computeFEInit(*mesh, {0, 0, 0}, false, 0);
    int count = mesh->getNumberOfNodes();
    auto prev = allocateVector<vectorReal>(count, "prev");
    auto curr = allocateVector<vectorReal>(count, "curr");
    auto host = Kokkos::create_mirror_view(curr);
    for (int n = 0; n < count; ++n) {
      float xyz[3];
      for (int c = 0; c < 3; ++c) xyz[c] = mesh->nodeCoord(n, c);
      host(n) = std::sin(xyz[0] + 0.3f * xyz[1] + 0.7f * xyz[2]);
    }
    Kokkos::deep_copy(curr, host);
    Kokkos::deep_copy(prev, host);
    auto term = allocateArray2D<arrayReal>(1, 1, "rhs");
    Kokkos::deep_copy(term, 0.f);
    auto elem = allocateVector<vectorInt>(1, "element");
    Kokkos::deep_copy(elem, 0);
    auto weights = allocateArray2D<arrayReal>(1, 27, "weights");
    Kokkos::deep_copy(weights, 0.f);
    solver::WavefieldAcoustic wf(prev, curr);
    solver::RhsAcoustic rhs(term, elem, weights);
    Solver::DataType data(wf, rhs);
    sem.resetGlobalVectors(count);
    sem.computeElementContributions_Acoustic_Flat(data);
    auto force = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sem.getForceVector(0));
    auto mass = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sem.getMassMatrixAcoustic());
    for (int n = 0; n < count; ++n) {
      float xyz[3];
      for (int c = 0; c < 3; ++c) xyz[c] = mesh->nodeCoord(n, c);
      results[which].push_back({xyz[0], xyz[1], xyz[2], force(n), mass(n)});
    }
    std::sort(results[which].begin(), results[which].end());
    ++which;
  }
  ASSERT_EQ(results[0].size(), results[1].size());
  for (size_t n = 0; n < results[0].size(); ++n) {
    for (int c = 0; c < 3; ++c) EXPECT_FLOAT_EQ(results[0][n][c], results[1][n][c]);
    EXPECT_NEAR(results[0][n][3], results[1][n][3], 2e-6f);
    EXPECT_NEAR(results[0][n][4], results[1][n][4], 1e-6f * std::abs(results[1][n][4]));
  }
}
TEST_F(Gmsh3DReader, AcousticSEM3DInitializesAndAdvancesConstantField) {
  using Integral = typename IntegralTypeSelector<2, IntegralType::MAKUTU>::type;
  using Mesh = model::ModelUnstruct<float, int>;
  using Solver = solver::SEMsolver<2, Integral, Mesh, false, utils::enums::physicType::kAcoustic>;
  // Interior flags give natural boundaries for the constant-field patch test.
  auto d = model::readGmsh3D(path, 2, {{10, {2, 1}}, {20, {3, 2}}},
                             {{1, model::BoundaryFlag::InteriorNode}, {2, model::BoundaryFlag::InteriorNode}});
  Mesh mesh(d);
  mesh.buildFaceConnectivity();
  Solver sem;
  sem.computeFEInit(mesh, {0, 0, 0}, false, 0);
  auto mass = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sem.getMassMatrixAcoustic());
  for (int i = 0; i < d.n_node_; ++i) EXPECT_GT(mass(i), 0);
  auto prev = allocateVector<vectorReal>(d.n_node_, "prev");
  auto curr = allocateVector<vectorReal>(d.n_node_, "curr");
  Kokkos::deep_copy(prev, 1.f);
  Kokkos::deep_copy(curr, 1.f);
  auto term = allocateArray2D<arrayReal>(1, 1, "rhs");
  Kokkos::deep_copy(term, 0.f);
  auto elem = allocateVector<vectorInt>(1, "element");
  Kokkos::deep_copy(elem, 0);
  auto weights = allocateArray2D<arrayReal>(1, 27, "weights");
  Kokkos::deep_copy(weights, 0.f);
  solver::WavefieldAcoustic wf(prev, curr);
  solver::RhsAcoustic rhs(term, elem, weights);
  Solver::DataType data(wf, rhs);
  sem.resetGlobalVectors(d.n_node_);
  sem.computeElementContributions_Acoustic_Flat(data);
  auto force = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sem.getForceVector(0));
  for (int i = 0; i < d.n_node_; ++i) EXPECT_NEAR(force(i), 0, 1e-6);
  sem.updateSolutionForward(0.001f, data);
  auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), prev);
  for (int i = 0; i < d.n_node_; ++i) EXPECT_NEAR(result(i), 1, 1e-6);
}
#endif
}  // namespace
