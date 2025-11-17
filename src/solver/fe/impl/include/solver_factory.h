#pragma once

#include <fe/Integrals.hpp>
#include <memory>

#include "sem_solver_base.h"

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
enum modelLocationType
{
  OnNodes,
  OnElements
};
enum physicType
{
  Acoustic,
  Elastic
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
    case MAKUTU:
      return "MAKUTU";
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

inline std::string to_string(modelLocationType loc)
{
  switch (loc)
  {
    case OnNodes:
      return "OnNodes";
    case OnElements:
      return "OnElements";
    default:
      return "Unknown";
  }
}

inline std::string to_string(physicType p)
{
  switch (p)
  {
    case Acoustic:
      return "Acoustic";
    case Elastic:
      return "Elastic";
    default:
      return "Unknown";
  }
}

std::unique_ptr<SEMSolverBase> createSolver(
    methodType const methodType, implemType const implemType,
    meshType const meshType, modelLocationType const modelLocation,
    physicType const physicType, int const order);
}  // namespace SolverFactory
