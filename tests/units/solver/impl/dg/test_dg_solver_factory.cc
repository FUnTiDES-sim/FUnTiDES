#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "solver.h"
#include "solver_factory.h"

#ifdef COMPILE_DG

namespace solver::solver_factory {


class DgSolverFactoryTest : public ::testing::Test {};

TEST_F(DgSolverFactoryTest, AcousticStructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDg, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                     utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, AcousticStructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDg, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                     utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, AcousticUnstructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDg, utils::enums::implemType::kMakutu, utils::enums::meshType::kUnstruct,
                     utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, AcousticUnstructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDg, utils::enums::implemType::kMakutu, utils::enums::meshType::kUnstruct,
                     utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, ElasticThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDg, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kElastic, 1),
               std::runtime_error);
}

TEST_F(DgSolverFactoryTest, AcoustoElasticThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDg, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustoElastic, 1),
               std::runtime_error);
}

TEST_F(DgSolverFactoryTest, UnknownImplTypeThrows) {
  auto unknownImpl = static_cast<utils::enums::implemType>(99);
  EXPECT_THROW(createSolver(utils::enums::methodType::kDg, unknownImpl, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 1),
               std::runtime_error);
}

}  // namespace solver::solver_factory

#endif  // COMPILE_DG
