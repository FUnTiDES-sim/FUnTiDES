#ifndef SOLVER_BASE_HPP_
#define SOLVER_BASE_HPP_

#include "dataType.hpp"
#include "model.hpp"

class SolverBase
{
public:
  // Constructor
  SolverBase() = default;

  // Destructor
  virtual ~SolverBase() = default;

  struct DataStruct
  {
    // Base structure for data passed to the solver
    virtual ~DataStruct() = default;
  };

  // Pure virtual function to compute one step of the solver
  virtual void computeOneStep( const float & dt,
                               const int & timeSample,
                               DataStruct & data ) = 0;

  // Pure virtual function to initialize finite element components
  virtual void computeFEInit( model::ModelApi< float, int > & mesh ) = 0;

  virtual void outputPnValues( const int &indexTimeStep,
                               int &i1,
                               int &myElementSource,
                               const ARRAY_REAL_VIEW &pnGlobal ) = 0;

};


#endif // SOLVER_BASE_HPP_
