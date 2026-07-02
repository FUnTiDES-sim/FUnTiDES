#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "solver.h"
#include "solver_factory.h"

#ifdef COMPILE_DG

namespace solver::fe::solver_factory {

namespace feenum = utils::enums;

class DgSolverFactoryTest : public ::testing::Test {};

TEST_F(DgSolverFactoryTest, AcousticStructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDg, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, AcousticStructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDg, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                     feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, AcousticUnstructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDg, feenum::implemType::kMakutu, feenum::meshType::kUnstruct,
                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, AcousticUnstructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDg, feenum::implemType::kMakutu, feenum::meshType::kUnstruct,
                     feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgSolverFactoryTest, ElasticThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDg, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 1),
               std::runtime_error);
}

TEST_F(DgSolverFactoryTest, AcoustoElasticThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDg, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustoElastic, 1),
               std::runtime_error);
}

TEST_F(DgSolverFactoryTest, UnknownImplTypeThrows) {
  auto unknownImpl = static_cast<feenum::implemType>(99);
  EXPECT_THROW(createSolver(feenum::methodType::kDg, unknownImpl, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 1),
               std::runtime_error);
}

}  // namespace solver::fe::solver_factory

#endif  // COMPILE_DG
