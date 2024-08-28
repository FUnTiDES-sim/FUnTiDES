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

#include "SEMQkGL.hpp"
#include "SEMmesh.hpp"

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
                       SEMmesh mesh );

  /**
   * @brief computeOneStep function:
   * init all FE components for computing mass and stiffness matrices
   */

  void computeOneStep ( const int & timeSample,
                        const int & order_a,
			const int & order_b,
			const int & order_c,
                        const int & nPointsPerElement,
                        const int & i1,
                        const int & i2,
                        SEMinfo & myInfo,
                        const arrayReal & myRHSTerm,
                        arrayReal const & myPnGlobal,
                        const vectorInt & myRhsElement );

  void outputPnValues ( SEMmesh mesh,
                        const int & indexTimeStep,
                        int & i1,
                        int & myElementSource,
                        const arrayReal & pnGlobal );

  void initFEarrays( SEMinfo & myInfo, SEMmesh mesh );

  void allocateFEarrays( SEMinfo & myInfo );

private:

  int order_a;
  int order_b;
  int order_c;
  SEMQkGL myQk;

  //shared arrays
  arrayInt globalNodesList;
  arrayReal globalNodesCoords;
  arrayReal globalNodesCoordsX;
  arrayReal globalNodesCoordsY;
  arrayReal globalNodesCoordsZ;
  vectorInt listOfInteriorNodes;
  vectorInt listOfBoundaryNodes;
  arrayInt faceInfos;
  arrayInt localFaceNodeToGlobalFaceNode;

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
  vectorReal ShGlobal;

  arrayInt listOfElementsByColor;
};
#endif //SEM_SOLVER_HPP_
