#pragma once

#include <fe/Integrals.hpp>
#include <memory>

#include "solver_base.h"

namespace SolverFactory
{
enum methodType
{
  SEM,
  DG
};
enum implemType
{
  MAKUTU,
  SHIVA
};
enum meshType
{
  Struct,
  Unstruct
};

std::unique_ptr<SolverBase> createSolver(methodType const methodType,
                                         implemType const implemType,
                                         meshType const meshType,
                                         int const order);
}  // namespace SolverFactory
