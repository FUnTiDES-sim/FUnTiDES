//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include "SEMsolver.hpp"

void SEMsolver::computeFEInit( SEMinfo & myInfo, SEMmesh mesh )
{
  order = myInfo.myOrderNumber;
  allocateFEarrays( myInfo );
  initFEarrays( myInfo, mesh );
}


// compute one step of the time dynamic wave equation solver
void SEMsolver::computeOneStep( const int & timeSample,
                                const int & order,
                                const int & nPointsPerElement,
                                const int & i1,
                                const int & i2,
                                SEMinfo & myInfo,
                                const arrayReal & RHS_Term,
                                arrayReal const & PN_Global,
                                const vectorInt & RHS_Element )
{
  CREATEVIEWS

  LOOPHEAD( myInfo.numberOfNodes, i )
  massMatrixGlobal[i]=0;
  yGlobal[i]=0;
  LOOPEND

  // update pnGLobal with right hade side
  LOOPHEAD( myInfo.myNumberOfRHS, i )
  int nodeRHS=globalNodesList( rhsElement[i], 0 );
  pnGlobal( nodeRHS, i2 )+=myInfo.myTimeStep*myInfo.myTimeStep*model[rhsElement[i]]
                          *model[rhsElement[i]]*rhsTerm( i, timeSample );
  LOOPEND

  // start main parallel section
  MAINLOOPHEAD( myInfo.numberOfElements, elementNumber )

  float massMatrixLocal[ROW];
  float pnLocal[ROW];
  float Y[ROW];

  // get pnGlobal to pnLocal
  for( int i=0; i<nPointsPerElement; i++ )
  {
    int localToGlobal=globalNodesList( elementNumber, i );
    pnLocal[i]=pnGlobal( localToGlobal, i2 );
  }

#ifdef USE_SEMCLASSIC
  myQkIntegrals.computeMassMatrixAndStiffnessVector(elementNumber,order,nPointsPerElement,
                                           globalNodesCoordsX,globalNodesCoordsY,globalNodesCoordsZ,
                                           weights,derivativeBasisFunction1D,
                                           massMatrixLocal,pnLocal,Y);
#endif

#ifdef USE_SEMOPTIM
  myQkIntegrals.computeMassMatrixAndStiffnessVector(elementNumber, order, nPointsPerElement,
                                           globalNodesCoordsX,globalNodesCoordsY,globalNodesCoordsZ,
                                           massMatrixLocal, pnLocal, Y);
#endif

#ifdef USE_SHIVA
  constexpr int ORDER=SEMinfo::myOrderNumber;
  myQkIntegrals.computeMassMatrixAndStiffnessVector<ORDER>(elementNumber, nPointsPerElement,
                                           globalNodesCoordsX,globalNodesCoordsY,globalNodesCoordsZ,
                                           massMatrixLocal, pnLocal, Y);
#endif


  //compute global mass Matrix and global stiffness vector
  for( int i=0; i<nPointsPerElement; i++ )
  {
    int gIndex=globalNodesList( elementNumber, i );
    massMatrixLocal[i]/=(model[elementNumber]*model[elementNumber]);
    ATOMICADD( massMatrixGlobal[gIndex], massMatrixLocal[i] );
    ATOMICADD( yGlobal[gIndex], Y[i] );
  }

  MAINLOOPEND

  // update pressure
  LOOPHEAD( myInfo.numberOfInteriorNodes, i )
  int I=listOfInteriorNodes[i];
  pnGlobal( I, i1 )=2*pnGlobal( I, i2 )-pnGlobal( I, i1 )
                   -myInfo.myTimeStep*myInfo.myTimeStep*yGlobal[I]/massMatrixGlobal[I];
  LOOPEND

  FENCE
}

void SEMsolver::outputPnValues(  SEMmesh mesh,
		                 const int & indexTimeStep,
                                 int & i1,
                                 int & myElementSource, 
                                 const arrayReal & pnGlobal)
{
    //writes debugging ascii file.
    if( indexTimeStep%50==0 )
    {   
      cout<<"TimeStep="<<indexTimeStep<<";  pnGlobal @ elementSource location "<<myElementSource
          <<" after computeOneStep = "<< pnGlobal(globalNodesList(myElementSource,0),i1)<<endl;
      #ifdef SEM_SAVE_SNAPSHOTS
      mesh.saveSnapShot( indexTimeStep, i1, pnGlobal );
      #endif
    }  
}

void SEMsolver::initFEarrays( SEMinfo & myInfo, SEMmesh mesh )
{
  //interior elements
  mesh.globalNodesList( myInfo.numberOfElements, globalNodesList );
  mesh.getListOfInteriorNodes( myInfo.numberOfInteriorNodes, listOfInteriorNodes );
  // mesh coordinates
  mesh.nodesCoordinates( myInfo.numberOfNodes, globalNodesCoords );
  mesh.nodesCoordinates( globalNodesCoordsX,globalNodesCoordsZ,globalNodesCoordsY);
  // boundary elements
  mesh.getListOfBoundaryNodes( myInfo.numberOfBoundaryNodes, listOfBoundaryNodes );
  mesh.getBoundaryFacesInfos( faceInfos );
  mesh.getLocalFaceNodeToGlobalFaceNode( localFaceNodeToGlobalFaceNode );
  // get model
  mesh.getModel( myInfo.numberOfElements, model );
  // get quadrature points
  myQkBasis.gaussLobattoQuadraturePoints( order, quadraturePoints );
  // get gauss-lobatto weights
  myQkBasis.gaussLobattoQuadratureWeights( order, weights );
  // get basis function and corresponding derivatives
  myQkBasis.getDerivativeBasisFunction1D( order, quadraturePoints, derivativeBasisFunction1D );
 
  // sort element by color
  #ifdef SEM_MESHCOLOR
  mesh.sortElementsByColor(myInfo.numberOfElementsByColor,listOfElementsByColor);
  printf("number of elements color red %d\n", myInfo.numberOfElementsByColor[0]);
  printf("number of elements color green %d\n", myInfo.numberOfElementsByColor[1]);
  printf("number of elements color blue %d\n", myInfo.numberOfElementsByColor[2]);
  printf("number of elements color yellow %d\n", myInfo.numberOfElementsByColor[3]);
  #endif

}

void SEMsolver::allocateFEarrays( SEMinfo & myInfo )
{
  int nbQuadraturePoints=(order+1)*(order+1)*(order+1);
  //interior elements
  cout<<"Allocate host memory for arrays in the solver ..."<<endl;
  globalNodesList=allocateArray2D< arrayInt >( myInfo.numberOfElements, myInfo.numberOfPointsPerElement, "globalNodesList" );
  listOfInteriorNodes=allocateVector< vectorInt >( myInfo.numberOfInteriorNodes, "listOfInteriorNodes" );
  
  // global coordinates
  globalNodesCoordsX=allocateArray2D< arrayReal >( myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsX");
  globalNodesCoordsY=allocateArray2D< arrayReal >( myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsY");
  globalNodesCoordsZ=allocateArray2D< arrayReal >( myInfo.numberOfElements, nbQuadraturePoints, "globalNodesCoordsZ");
  
  listOfBoundaryNodes=allocateVector< vectorInt >( myInfo.numberOfBoundaryNodes, "listOfBoundaryNodes" );

  faceInfos=allocateArray2D< arrayInt >( myInfo.numberOfBoundaryFaces, 2+(order+1), "faceInfos" );
  localFaceNodeToGlobalFaceNode=allocateArray2D< arrayInt >( myInfo.numberOfBoundaryFaces, order+1, "localFaceNodeToGlobalFaceNode" );

  model=allocateVector< vectorReal >( myInfo.numberOfElements, "model" );

  quadraturePoints=allocateVector< vectorDouble >( order+1, "quadraturePoints" );

  weights=allocateVector< vectorDouble >( order+1, "weights" );

  derivativeBasisFunction1D=allocateArray2D< arrayDouble >( order+1, order+1, "derivativeBasisFunction1D" );

  //shared arrays
  massMatrixGlobal=allocateVector< vectorReal >( myInfo.numberOfNodes, "massMatrixGlobal" );
  yGlobal=allocateVector< vectorReal >( myInfo.numberOfNodes, "yGlobal" );
  ShGlobal=allocateVector< vectorReal >( myInfo.numberOfBoundaryNodes, "ShGlobal" );

  #ifdef SEM_MESHCOLOR
  //allocate list of elements by color
  listOfElementsByColor=allocateArray2D<arrayInt>(myInfo.numberOfColors, myInfo.numberMaxOfElementsByColor, "listOfElemByColor");
  #endif
}
