#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

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
    s = createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, p.mesh, p.location, p.physic, p.order);
  });

  ASSERT_NE(s, nullptr);
}

INSTANTIATE_TEST_SUITE_P(FullMatrix, SolverFactoryMatrixTest,
                         ::testing::Values(
                             /* Structured Mesh: Covering both OnNodes and OnElements for all physics */
                             SolverTestParams{utils::enums::physicType::kAcoustic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnNodes, 1},
                             SolverTestParams{utils::enums::physicType::kAcoustic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnElements, 1},
                             SolverTestParams{utils::enums::physicType::kElastic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnNodes, 2},
                             SolverTestParams{utils::enums::physicType::kElastic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnElements, 2},
                             SolverTestParams{utils::enums::physicType::kAcoustoElastic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnNodes, 1},
                             SolverTestParams{utils::enums::physicType::kAcoustoElastic, utils::enums::meshType::kStruct,
                                              utils::enums::modelLocationType::kOnElements, 1},

                             /* Unstructured Mesh: Covering both OnNodes and OnElements for all physics */
                             SolverTestParams{utils::enums::physicType::kAcoustic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnNodes, 1},
                             SolverTestParams{utils::enums::physicType::kAcoustic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnElements, 1},
                             SolverTestParams{utils::enums::physicType::kElastic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnNodes, 2},
                             SolverTestParams{utils::enums::physicType::kElastic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnElements, 2},
                             SolverTestParams{utils::enums::physicType::kAcoustoElastic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnNodes, 1},
                             SolverTestParams{utils::enums::physicType::kAcoustoElastic, utils::enums::meshType::kUnstruct,
                                              utils::enums::modelLocationType::kOnElements, 1}));

class SolverFactoryErrorTest : public ::testing::Test {};

/* Covers the 'throw std::runtime_error("Unknown physics type")' line */
TEST_F(SolverFactoryErrorTest, ThrowsOnUnknownPhysics) {
  auto unknownPhysic = static_cast<utils::enums::physicType>(99);
  EXPECT_THROW(createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnNodes, unknownPhysic, 1),
               std::runtime_error);
}

/* Covers orderDispatch failure paths */
TEST_F(SolverFactoryErrorTest, ThrowsOnInvalidPolynomialOrder) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
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

/* Covers the no-op default bodies of the Solver base API, i.e. the hooks that
   only some physics/methods override (SEM acoustic overrides none of them) */
TEST_F(SolverFactoryErrorTest, BaseClassOptionalHooksAreNoOps) {
  auto s = createSolver(utils::enums::methodType::kSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                        utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1);
  ASSERT_NE(s, nullptr);

  EXPECT_NO_THROW(s->setZBoundary(0.f));
  EXPECT_NO_THROW(s->setInterfacePropertyConvention(utils::enums::interfacePropertyConvention::kFluidOnInterfaceNodes));
  EXPECT_EQ(s->getInterfaceCouplingCoeff(0).size(), 0u);
}

}  // namespace solver::solver_factory

#endif  // COMPILE_SEM
