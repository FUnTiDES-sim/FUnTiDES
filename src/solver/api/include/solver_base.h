#ifndef SOLVER_BASE_HPP_
#define SOLVER_BASE_HPP_

#include <data_type.h>
#include <model.h>

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

    virtual void print() const = 0;
  };

  // Pure virtual function to compute one step of the solver
  virtual void computeOneStep(const float &dt, const int &timeSample,
                              DataStruct &data) = 0;
};
#endif  // SOLVER_BASE_HPP_
