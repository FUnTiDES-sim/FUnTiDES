#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "solver.h"
#include "solver_factory.h"

#ifdef COMPILE_DG_SEM

namespace solver::fe::solver_factory {

namespace feenum = utils::enums;

class DgSemSolverFactoryTest : public ::testing::Test {};

TEST_F(DgSemSolverFactoryTest, AcousticStructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, AcousticStructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                     feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, AcousticUnstructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgSem, feenum::implemType::kMakutu, feenum::meshType::kUnstruct,
                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, AcousticUnstructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgSem, feenum::implemType::kMakutu, feenum::meshType::kUnstruct,
                     feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSemSolverFactoryTest, ElasticThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDgSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1),
               std::runtime_error);
}

TEST_F(DgSemSolverFactoryTest, AcoustoElasticThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDgSem, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustoElastic, 1),
               std::runtime_error);
}

TEST_F(DgSemSolverFactoryTest, UnknownImplTypeThrows) {
  auto unknownImpl = static_cast<feenum::implemType>(99);
  EXPECT_THROW(createSolver(feenum::methodType::kDgSem, unknownImpl, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1),
               std::runtime_error);
}

}  // namespace solver::fe::solver_factory

#endif  // COMPILE_DG_SEM
