/**
 * @file test_sponge_taper.cc
 * @brief Unit tests for SEMsolver::initSpongeValues.
 *
 * The taper is defined on a length normalised by the sponge thickness, so the
 * profile must span the whole layer whatever the units of the mesh, and must
 * be measured from the actual domain box rather than from the origin of the
 * coordinate system.
 */
#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "Integrals.h"
#include "cartesian_struct_builder.h"
#include "common_macros.h"
#include "data_type.h"
#include "model_struct.h"
#include "sem_solver_impl.h"

namespace solver {
namespace fe {
namespace test {

namespace {

constexpr int kOrder = 1;
// Eight elements put a node exactly half way through the sponge layer, which is
// where the tests probe the decay.
constexpr int kNElem = 8;
constexpr float kDomain = 2000.0f;
constexpr float kSponge = 500.0f;

using MeshType = model::ModelStruct<float, int, kOrder>;
using IntType = typename IntegralTypeSelector<kOrder, IntegralType::MAKUTU>::type;
using SolverT = SEMsolver<kOrder, IntType, MeshType, false, physicType::kAcoustic>;

/// Build a kNElem^3 Cartesian mesh whose lower corner sits at (ox, oy, oz).
/// The global box is the local one, so all six faces are external.
std::shared_ptr<model::ModelApi<float, int>> makeMesh(float ox, float oy, float oz, bool free_surface = false) {
  model::CartesianStructBuilder<float, int, kOrder> builder(kNElem, kDomain, kNElem, kDomain, kNElem, kDomain, false,
                                                            false, ox, oy, oz, kDomain, kDomain, kDomain, ox, oy, oz);
  return builder.getModel(free_surface);
}

}  // namespace

TEST(SpongeTaperTest, NoSpongeLeavesEveryCoefficientAtOne) {
  auto mesh = makeMesh(0.0f, 0.0f, 0.0f);
  SolverT solver;
  solver.computeFEInit(*mesh, {0.0f, 0.0f, 0.0f}, false, 0.333f);

  auto const& taper = solver.getSpongeTaperCoeff();
  for (int n = 0; n < mesh->getNumberOfNodes(); ++n) EXPECT_FLOAT_EQ(taper(n), 1.0f);
}

TEST(SpongeTaperTest, TaperAttenuatesAcrossTheWholeLayer) {
  auto mesh = makeMesh(0.0f, 0.0f, 0.0f);
  SolverT solver;
  solver.computeFEInit(*mesh, {kSponge, kSponge, kSponge}, false, 0.333f);

  auto const& taper = solver.getSpongeTaperCoeff();

  // A node on the boundary must be damped, and the damping must decrease
  // monotonically with the distance to the boundary. The regression this
  // guards against made every coefficient 1.0 because the decay length was
  // read as an absolute distance of a few centimetres.
  float boundary = 1.0f, mid_layer = 1.0f, interior = 1.0f;
  for (int n = 0; n < mesh->getNumberOfNodes(); ++n) {
    const float x = mesh->nodeCoord(n, 0);
    const float y = mesh->nodeCoord(n, 1);
    const float z = mesh->nodeCoord(n, 2);
    const float d = std::min({x, kDomain - x, y, kDomain - y, z, kDomain - z});
    if (d == 0.0f) boundary = std::min(boundary, taper(n));
    if (std::abs(d - kSponge * 0.5f) < 1.0f) mid_layer = std::min(mid_layer, taper(n));
    if (d >= kSponge) interior = std::min(interior, taper(n));
  }

  EXPECT_LT(boundary, 0.9f);          // 1 / (1 + sigma_max) with sigma_max = 0.15
  EXPECT_GT(mid_layer, boundary);     // damping decays inwards
  EXPECT_LT(mid_layer, 1.0f);         // but is still active in the middle of the layer
  EXPECT_FLOAT_EQ(interior, 1.0f);    // and stops at the inner edge
}

TEST(SpongeTaperTest, TaperFollowsAShiftedDomainOrigin) {
  const float shift = 1.0e5f;  // a georeferenced or partitioned mesh
  auto mesh = makeMesh(shift, shift, shift);
  SolverT solver;
  solver.computeFEInit(*mesh, {kSponge, kSponge, kSponge}, false, 0.333f);

  auto const& taper = solver.getSpongeTaperCoeff();
  for (int n = 0; n < mesh->getNumberOfNodes(); ++n) {
    const float x = mesh->nodeCoord(n, 0) - shift;
    const float y = mesh->nodeCoord(n, 1) - shift;
    const float z = mesh->nodeCoord(n, 2) - shift;
    const float d = std::min({x, kDomain - x, y, kDomain - y, z, kDomain - z});
    if (d >= kSponge) {
      EXPECT_FLOAT_EQ(taper(n), 1.0f) << "interior node at distance " << d << " was damped";
    } else {
      EXPECT_LT(taper(n), 1.0f) << "sponge node at distance " << d << " was not damped";
    }
  }
}

TEST(SpongeTaperTest, FreeSurfaceIsNotDampedUnlessAsked) {
  // The Cartesian builder puts the free surface on the z = origin + lz face.
  auto mesh = makeMesh(0.0f, 0.0f, 0.0f, true);
  SolverT kept, damped;
  kept.computeFEInit(*mesh, {0.0f, 0.0f, kSponge}, false, 0.333f);
  damped.computeFEInit(*mesh, {0.0f, 0.0f, kSponge}, true, 0.333f);

  auto const& taper_kept = kept.getSpongeTaperCoeff();
  auto const& taper_damped = damped.getSpongeTaperCoeff();
  bool saw_surface = false;
  for (int n = 0; n < mesh->getNumberOfNodes(); ++n) {
    if (mesh->nodeCoord(n, 2) < kDomain - 1.0f) continue;
    saw_surface = true;
    EXPECT_FLOAT_EQ(taper_kept(n), 1.0f);
    EXPECT_LT(taper_damped(n), 1.0f);
  }
  EXPECT_TRUE(saw_surface);
}

TEST(SpongeTaperTest, InternalPartitionFacesAreNotDamped) {
  // Middle slab of a three-rank decomposition along x. Both x faces are shared
  // with a neighbour, so a sponge there would eat the wavefield in the middle
  // of the global domain.
  model::CartesianStructBuilder<float, int, kOrder> builder(kNElem, kDomain, kNElem, kDomain, kNElem, kDomain, false,
                                                            false, kDomain, 0.0f, 0.0f, 3.0f * kDomain, kDomain,
                                                            kDomain, 0.0f, 0.0f, 0.0f);
  auto mesh = builder.getModel(false);

  SolverT along_x;
  along_x.computeFEInit(*mesh, {kSponge, 0.0f, 0.0f}, false, 0.333f);
  auto const& taper_x = along_x.getSpongeTaperCoeff();
  for (int n = 0; n < mesh->getNumberOfNodes(); ++n)
    EXPECT_FLOAT_EQ(taper_x(n), 1.0f) << "node at x = " << mesh->nodeCoord(n, 0) << " was damped by an internal face";

  // The y faces of the same slab are global, so they must still get a layer.
  SolverT along_y;
  along_y.computeFEInit(*mesh, {0.0f, kSponge, 0.0f}, false, 0.333f);
  auto const& taper_y = along_y.getSpongeTaperCoeff();
  bool damped = false;
  for (int n = 0; n < mesh->getNumberOfNodes(); ++n)
    if (mesh->nodeCoord(n, 1) < 1.0f) damped = damped || taper_y(n) < 1.0f;
  EXPECT_TRUE(damped) << "the global y face of the slab was left undamped";
}

}  // namespace test
}  // namespace fe
}  // namespace solver
