#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>
#include <vector>

#include "solver.h"
#include "solver_factory.h"

#ifdef COMPILE_SEM

namespace solver::solver_factory {

struct SolverTestParams {
  utils::enums::physicType physic;
  utils::enums::meshType mesh;
  utils::enums::modelLocationType location;
  int order;
};

class SolverFactoryMatrixTest : public ::testing::TestWithParam<SolverTestParams> {};

/* Validates all physics, mesh types, and model locations (OnNodes vs OnElements) */
TEST_P(SolverFactoryMatrixTest, ProducesNonNullSolver) {
  auto const& p = GetParam();

  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, p.mesh, p.location, p.physic,
                     p.order);
  });

  ASSERT_NE(s, nullptr);
}

namespace {

/* Both mesh types and both model locations, for one physics at one order */
void AppendPhysics(std::vector<SolverTestParams>& params, utils::enums::physicType physic, int order) {
  for (auto mesh : {utils::enums::meshType::kStruct, utils::enums::meshType::kUnstruct})
    for (auto location : {utils::enums::modelLocationType::kOnNodes, utils::enums::modelLocationType::kOnElements})
      params.push_back(SolverTestParams{physic, mesh, location, order});
}

/* Only the physics compiled in this build */
std::vector<SolverTestParams> CompiledPhysicsMatrix() {
  std::vector<SolverTestParams> params;
#ifdef COMPILE_ACOUSTIC
  AppendPhysics(params, utils::enums::physicType::kAcoustic, 1);
#endif
#ifdef COMPILE_ELASTIC
  AppendPhysics(params, utils::enums::physicType::kElastic, 2);
#endif
#ifdef COMPILE_ACOUSTOELASTIC
  AppendPhysics(params, utils::enums::physicType::kAcoustoElastic, 1);
#endif
  return params;
}

}  // namespace

INSTANTIATE_TEST_SUITE_P(FullMatrix, SolverFactoryMatrixTest, ::testing::ValuesIn(CompiledPhysicsMatrix()));
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(SolverFactoryMatrixTest);

class SolverFactoryErrorTest : public ::testing::Test {};

/* Covers the 'throw std::runtime_error("Unknown physics type")' line */
TEST_F(SolverFactoryErrorTest, ThrowsOnUnknownPhysics) {
  auto unknownPhysic = static_cast<utils::enums::physicType>(99);
  EXPECT_THROW(
      createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                   utils::enums::modelLocationType::kOnNodes, unknownPhysic, 1),
      std::runtime_error);
}

/* Covers orderDispatch failure paths */
TEST_F(SolverFactoryErrorTest, ThrowsOnInvalidPolynomialOrder) {
  EXPECT_THROW(
      createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                   utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 999),
      std::runtime_error);
}

/* Covers unsupported implementation fallback */
TEST_F(SolverFactoryErrorTest, ThrowsOnUnsupportedMethodOrImpl) {
  auto unsupportedMethod = static_cast<utils::enums::methodType>(99);
  EXPECT_THROW(createSolver(unsupportedMethod, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1),
               std::runtime_error);
}

/* Covers "Unknown implementation type" inside the kSem switch */
TEST_F(SolverFactoryErrorTest, ThrowsOnUnknownImplTypeForSem) {
  auto unknownImpl = static_cast<utils::enums::implemType>(99);
  EXPECT_THROW(createSolver(utils::enums::methodType::kSem, unknownImpl, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1),
               std::runtime_error);
}

#ifdef COMPILE_ACOUSTIC
/* Covers the no-op default bodies of the Solver base API, i.e. the hooks that
   only some physics/methods override (SEM acoustic overrides none of them) */
TEST_F(SolverFactoryErrorTest, BaseClassOptionalHooksAreNoOps) {
  auto s =
      createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                   utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1);
  ASSERT_NE(s, nullptr);

  EXPECT_NO_THROW(s->setZBoundary(0.f));
  EXPECT_NO_THROW(s->setInterfacePropertyConvention(utils::enums::interfacePropertyConvention::kFluidOnInterfaceNodes));
  EXPECT_EQ(s->getInterfaceCouplingCoeff(0).size(), 0u);
}
#endif  // COMPILE_ACOUSTIC

#if !defined(COMPILE_ACOUSTIC) || !defined(COMPILE_ELASTIC) || !defined(COMPILE_ACOUSTOELASTIC)
/* A physics left out of the build is reported, not silently replaced */
TEST_F(SolverFactoryErrorTest, ThrowsOnPhysicsNotCompiled) {
#ifndef COMPILE_ACOUSTIC
  auto const missing = utils::enums::physicType::kAcoustic;
#elif !defined(COMPILE_ELASTIC)
  auto const missing = utils::enums::physicType::kElastic;
#else
  auto const missing = utils::enums::physicType::kAcoustoElastic;
#endif
  EXPECT_THROW(
      createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                   utils::enums::modelLocationType::kOnNodes, missing, 1),
      std::runtime_error);
}
#endif

}  // namespace solver::solver_factory

#endif  // COMPILE_SEM
