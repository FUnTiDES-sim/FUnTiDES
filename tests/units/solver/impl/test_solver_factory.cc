#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "solver.h"
#include "solver_factory.h"

namespace solver::fe::solver_factory {

namespace feenum = utils::enums;

struct SolverTestParams {
  feenum::physicType physic;
  feenum::meshType mesh;
  feenum::modelLocationType location;
  int order;
};

class SolverFactoryMatrixTest : public ::testing::TestWithParam<SolverTestParams> {};

/* Validates all physics, mesh types, and model locations (OnNodes vs OnElements) */
TEST_P(SolverFactoryMatrixTest, ProducesNonNullSolver) {
  auto const& p = GetParam();

  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, p.mesh, p.location, p.physic, p.order);
  });

  ASSERT_NE(s, nullptr);
}

INSTANTIATE_TEST_SUITE_P(FullMatrix, SolverFactoryMatrixTest,
                         ::testing::Values(
                             /* Structured Mesh: Covering both OnNodes and OnElements for all physics */
                             SolverTestParams{feenum::physicType::kAcoustic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnNodes, 1},
                             SolverTestParams{feenum::physicType::kAcoustic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnElements, 1},
                             SolverTestParams{feenum::physicType::kElastic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnNodes, 2},
                             SolverTestParams{feenum::physicType::kElastic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnElements, 2},
                             SolverTestParams{feenum::physicType::kAcoustoElastic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnNodes, 1},
                             SolverTestParams{feenum::physicType::kAcoustoElastic, feenum::meshType::kStruct,
                                              feenum::modelLocationType::kOnElements, 1},

                             /* Unstructured Mesh: Covering both OnNodes and OnElements for all physics */
                             SolverTestParams{feenum::physicType::kAcoustic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnNodes, 1},
                             SolverTestParams{feenum::physicType::kAcoustic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnElements, 1},
                             SolverTestParams{feenum::physicType::kElastic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnNodes, 2},
                             SolverTestParams{feenum::physicType::kElastic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnElements, 2},
                             SolverTestParams{feenum::physicType::kAcoustoElastic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnNodes, 1},
                             SolverTestParams{feenum::physicType::kAcoustoElastic, feenum::meshType::kUnstruct,
                                              feenum::modelLocationType::kOnElements, 1}));

class SolverFactoryErrorTest : public ::testing::Test {};

/* Covers the 'throw std::runtime_error("Unknown physics type")' line */
TEST_F(SolverFactoryErrorTest, ThrowsOnUnknownPhysics) {
  auto unknownPhysic = static_cast<feenum::physicType>(99);
  EXPECT_THROW(createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnNodes, unknownPhysic, 1),
               std::runtime_error);
}

/* Covers orderDispatch failure paths */
TEST_F(SolverFactoryErrorTest, ThrowsOnInvalidPolynomialOrder) {
  EXPECT_THROW(createSolver(feenum::methodType::kSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 999),
               std::runtime_error);
}

/* Covers unsupported implementation fallback */
TEST_F(SolverFactoryErrorTest, ThrowsOnUnsupportedMethodOrImpl) {
  auto unsupportedMethod = static_cast<feenum::methodType>(99);
  EXPECT_THROW(createSolver(unsupportedMethod, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 1),
               std::runtime_error);
}

}  // namespace solver::fe::solver_factory
