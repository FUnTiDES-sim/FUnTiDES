//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.hpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#ifndef SEM_SOLVER_HPP_
#define SEM_SOLVER_HPP_

//#include "SEMQkGL.hpp"
#include <BasisFunctions.hpp>
#include <Integrals.hpp>
#include <model.hpp>
#include "commonMacros.hpp"
#include "SEMdata.hpp"
#include "SEMQkGLBasisFunctions.hpp"

#ifdef USE_CALIPER
#include <caliper/cali.h>
#endif

class SEMsolver
{
public:

  PROXY_HOST_DEVICE SEMsolver(){};
  PROXY_HOST_DEVICE ~SEMsolver(){};

  /**
   * @brief computeFEInit function:
   * init all FE components for computing mass and stiffness matrices
   */
  void computeFEInit ( SEMinfo & myInfo,
                       Mesh mesh );

  /**
   * @brief computeOneStep function:
   * init all FE components for computing mass and stiffness matrices
   */

  void computeOneStep ( const int & timeSample,
                        const int & order,
                        const int & nPointsPerElement,
                        const int & i1,
                        const int & i2,
                        SEMinfo & myInfo,
                        const arrayReal & myRHSTerm,
                        arrayReal const & myPnGlobal,
                        const vectorInt & myRhsElement );

  void outputPnValues ( Mesh mesh,
                        const int & indexTimeStep,
                        int & i1,
                        int & myElementSource,
                        const arrayReal & pnGlobal );

  void initFEarrays( SEMinfo & myInfo, Mesh mesh );

  void allocateFEarrays( SEMinfo & myInfo );

private:

  int order;
  SEMQkGLIntegrals myQkIntegrals;

  //shared arrays
  arrayInt globalNodesList;
  arrayReal globalNodesCoordsX;
  arrayReal globalNodesCoordsY;
  arrayReal globalNodesCoordsZ;
  vectorInt listOfInteriorNodes;

  // get model
  vectorReal model;

  // get quadrature points and weights
  vectorDouble quadraturePoints;
  vectorDouble weights;

  // get basis function and corresponding derivatives
  arrayDouble derivativeBasisFunction1D;

  //shared arrays
  vectorReal massMatrixGlobal;
  vectorReal yGlobal;

};
#endif //SEM_SOLVER_HPP_
