#pragma once

#include "SolverBase.hpp"
#include "fe/Integrals.hpp"

#include <memory>

std::unique_ptr<SolverBase> createSolver( int const physicsType,
                                          int const methodType,
                                          int const order );

