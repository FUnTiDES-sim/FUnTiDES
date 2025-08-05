//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include <SEMsolver.hpp>
#include <argsparse.hpp>
#include <utils.hpp>
#ifdef USE_CALIPER
#include <caliper/cali.h>
#endif

/**
 * @class SEMproxy
 */

class SEMproxy {
public:
  /**
   * @brief Constructor of the SEMproxy class
   */
  SEMproxy(int argc, char *argv[]);

  SEMproxy(int ex, int ey, int ez, float lx)
      : myMesh(ex, ey, ez, lx, lx, lx, order), mySolver(myMesh) {}

  /**
   * @brief Destructor of the SEMproxy class
   */
  ~SEMproxy(){};

  /**
   * @brief Initialize the simulation.
   * @post run()
   */
  void initFiniteElem() {
    // allocate arrays and vectors
    init_arrays();
    // initialize source and RHS
    init_source();
  };

  // void saveCtrlSlice(int iteration, int i);

  /**
   * @brief Run the simulation.
   * @pre This must be called after init()
   * @post Optional printout performance resutls
   */
  void run();

  /**
  * Extract an XY slice from a Kokkos 1D View representing a 3D cubic array
  *
  * For Kokkos:
  * Extract XY slice using subview (more efficient, zero-copy)
  * Returns a subview that shares memory with the original array
  *
  * @param array_1d: Input Kokkos 1D View containing 3D data (stored in X-Y-Z order)
  * @param size: Size of each dimension (assuming cubic array: size x size x size)
  * @param z: Z-level to extract (0 to size-1)
  * @return: Kokkos 1D View containing the XY slice
  */
  VECTOR_REAL_VIEW extractXYSlice(const VECTOR_REAL_VIEW& array, int size, int z);

  /**
  * Save slice in gnuplot matrix format (default - best for gnuplot)
  * Format: space-separated matrix with blank lines between rows for 3D plotting
  */
  void saveSlice(const VECTOR_REAL_VIEW& host_slice,
                int size, const std::string& filepath);

private:
  int i1 = 0;
  int i2 = 1;

  const int myNumberOfRHS = 1;
  const float myTimeStep = 0.001;
  const float f0 = 10.;
  const float myTimeMax = 1.5;
  const int sourceOrder = 1;
  const static int order = 3;
  int myNumSamples = myTimeMax / myTimeStep;
  int myElementSource = 0;

  Mesh myMesh;

  SEMsolver mySolver;
  SolverUtils myUtils;

  // arrays
  arrayReal myRHSTerm;
  arrayReal pnGlobal;
  vectorInt rhsElement;
  arrayReal rhsWeights;

  // initialize source and RHS
  void init_source();

  // allocate arrays and vectors
  void init_arrays();
};

#endif /* SEMPROXY_HPP_ */
