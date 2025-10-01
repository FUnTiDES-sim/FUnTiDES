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
  GEOS,
  SHIVA
};
enum meshType
{
  Struct,
  Unstruct
};

inline std::string to_string(methodType m)
{
  switch (m)
  {
    case SEM:
      return "SEM";
    case DG:
      return "DG";
    default:
      return "Unknown";
  }
}

inline std::string to_string(implemType i)
{
  switch (i)
  {
    case GEOS:
      return "GEOS";
    case SHIVA:
      return "SHIVA";
    default:
      return "Unknown";
  }
}

inline std::string to_string(meshType m)
{
  switch (m)
  {
    case Struct:
      return "Struct";
    case Unstruct:
      return "Unstruct";
    default:
      return "Unknown";
  }
}

std::unique_ptr<SolverBase> createSolver(methodType const methodType,
                                         implemType const implemType,
                                         meshType const meshType,
                                         int const order);
}  // namespace SolverFactory
