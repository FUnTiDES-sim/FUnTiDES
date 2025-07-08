#ifndef SOLVER_BASE_HPP_
#define SOLVER_BASE_HPP_

#include "dataType.hpp"
#include "SEMdata.hpp"
#include "model.hpp"

class SolverBase
{
public:
  // Constructor
  SolverBase() = default;

  // Destructor
  virtual ~SolverBase() = default;

  // Pure virtual function to compute one step of the solver
  virtual void computeOneStep(const int & timeSample,
                              const int & nPointsPerElement,
                              const int & i1,
                              const int & i2,
                              SEMinfo const & myInfo,
                              const ARRAY_REAL_VIEW & rhsTerm,
                              const ARRAY_REAL_VIEW & pnGlobal,
                              const VECTOR_INT_VIEW & rhsElement ) = 0;

  // Pure virtual function to initialize finite element components
  virtual void computeFEInit( SEMinfo const & myInfo, 
                              Mesh const & mesh) = 0;

  virtual void outputPnValues( Mesh mesh, 
                               const int &indexTimeStep, 
                               int &i1,
                               int &myElementSource, 
                               const ARRAY_REAL_VIEW &pnGlobal) = 0;

};


#endif // SOLVER_BASE_HPP_