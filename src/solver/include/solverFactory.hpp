#pragma once

#include "SolverBase.hpp"
#include "fe/Integrals.hpp"
#include <memory>

namespace SolverFactory
{
    enum methodType { SEM, DG };
    enum implemType { CLASSIC, GEOS, OPTIM, SHIVA };
    enum meshType   { Struct, Unstruct};

    std::unique_ptr<SolverBase> createSolver( methodType const methodType,
                                              implemType const implemType,
                                              meshType const meshType,
                                              int const order );
}

