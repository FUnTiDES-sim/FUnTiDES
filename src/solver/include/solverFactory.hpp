#ifndef SOLVERFACTORY_H_
#define SOLVERFACTORY_H_

#pragma once

#include "SolverBase.hpp"
#include "fe/Integrals.hpp"
#include <memory>

namespace SolverFactory
{
    enum methodType { SEM, DG };
    enum implemType { CLASSIC, GEOS, OPTIM, SHIVA };

    std::unique_ptr<SolverBase> createSolver( methodType const methodType,
                                              implemType const implemType,
                                              int const order );
}


#endif // SOLVERFACTORY_H_
