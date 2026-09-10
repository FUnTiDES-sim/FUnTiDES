#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "solver.h"
#include "solver_factory.h"

#ifdef COMPILE_DG_PADAPTIVE

namespace solver::fe::solver_factory {

namespace feenum = utils::enums;

class DgPAdaptiveSolverFactoryTest : public ::testing::Test {};

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticStructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticStructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                     feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticUnstructOnElements_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kUnstruct,
                     feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcousticUnstructOnNodes_NonNull) {
  std::unique_ptr<Solver> s;
  EXPECT_NO_THROW({
    s = createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kUnstruct,
                     feenum::modelLocationType::kOnNodes, feenum::physicType::kAcoustic, 2, 1);
  });
  EXPECT_NE(s, nullptr);
}

TEST_F(DgPAdaptiveSolverFactoryTest, ElasticThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kElastic, 2, 1),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, AcoustoElasticThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustoElastic, 2, 1),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, UnknownImplTypeThrows) {
  auto unknownImpl = static_cast<feenum::implemType>(99);
  EXPECT_THROW(createSolver(feenum::methodType::kDgPAdaptive, unknownImpl, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 2, 1),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, OrderMinZeroThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 2, 0),
               std::runtime_error);
}

TEST_F(DgPAdaptiveSolverFactoryTest, OrderMinGreaterOrEqualOrderMaxThrows) {
  EXPECT_THROW(createSolver(feenum::methodType::kDgPAdaptive, feenum::implemType::kMakutu, feenum::meshType::kStruct,
                            feenum::modelLocationType::kOnElements, feenum::physicType::kAcoustic, 2, 2),
               std::runtime_error);
}

}  // namespace solver::fe::solver_factory

#endif  // COMPILE_DG_PADAPTIVE
