#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "solver.h"
#include "solver_factory.h"

#ifdef COMPILE_DG_PADAPTIVE

namespace solver::solver_factory {


class DgPAdaptiveSolverFactoryTest : public ::testing::Test {};

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticStructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                     utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticStructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                     utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticUnstructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kUnstruct,
                     utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticUnstructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kUnstruct,
                     utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, ElasticThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kElastic, 2, 1),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcoustoElasticThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustoElastic, 2, 1),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, UnknownImplTypeThrows) {
  auto unknownImpl = static_cast<utils::enums::implemType>(99);
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgPAdaptive, unknownImpl, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 2, 1),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, OrderMinZeroThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 2, 0),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, OrderMinGreaterOrEqualOrderMaxThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgPAdaptive, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 2, 2),
               std::runtime_error);
}

}  // namespace solver::solver_factory

#endif  // COMPILE_DG_PADAPTIVE
