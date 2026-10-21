#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "solver.h"
#include "solver_factory.h"

#ifdef COMPILE_DG_SEM

namespace solver::solver_factory {


class DgSemSolverFactoryTest : public ::testing::Test {};

TEST_F(DgSemSolverFactoryTest, AcousticStructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                     utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, AcousticStructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                     utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, AcousticUnstructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kUnstruct,
                     utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, AcousticUnstructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(utils::enums::methodType::kDgSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kUnstruct,
                     utils::enums::modelLocationType::kOnNodes, utils::enums::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, ElasticThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kElastic, 1),
               std::runtime_error);
}

TEST_F(DgSemSolverFactoryTest, AcoustoElasticThrows) {
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgSem, utils::enums::implemType::kMakutu, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustoElastic, 1),
               std::runtime_error);
}

TEST_F(DgSemSolverFactoryTest, UnknownImplTypeThrows) {
  auto unknownImpl = static_cast<utils::enums::implemType>(99);
  EXPECT_THROW(createSolver(utils::enums::methodType::kDgSem, unknownImpl, utils::enums::meshType::kStruct,
                            utils::enums::modelLocationType::kOnElements, utils::enums::physicType::kAcoustic, 1),
               std::runtime_error);
}

}  // namespace solver::solver_factory

#endif  // COMPILE_DG_SEM
