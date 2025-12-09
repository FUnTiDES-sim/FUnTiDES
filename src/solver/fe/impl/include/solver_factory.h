#pragma once

#include <fe/Integrals.hpp>
#include <memory>

#include "sem_enums.h"
#include "sem_solver_base.h"

using namespace solver::fe;

namespace SolverFactory
{

std::unique_ptr<SEMSolverBase> createSolver(
    methodType const methodType, implemType const implemType,
    meshType const meshType, modelLocationType const modelLocation,
    physicType const physicType, int const order);

}  // namespace SolverFactory
