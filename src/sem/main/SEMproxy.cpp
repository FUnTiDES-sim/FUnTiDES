//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "SEMproxy.hpp"

SEMproxy::SEMproxy(int argc, char *argv[])
{
  int ex = (argc > 1)? std::stoi(argv[1]) : 100; 
  int ey = ex;
  int ez = ex;

  float lx = (argc > 2)? std::stof(argv[2]) : 2000; 
  float ly = lx;
  float lz = lx;

  SEMmesh simpleMesh {ex, ey, ez, lx, ly, lz, myInfo.myOrderNumber_a, myInfo.myOrderNumber_b, myInfo.myOrderNumber_c};
  myMesh = simpleMesh;
}

// Initialize the simulation.
void SEMproxy::initFiniteElem()
{
  // get information from mesh
  getMeshInfo();

  // allocate arrays and vectors
  init_arrays();

  // initialize source and RHS
  init_source();

  mySolver.computeFEInit( myInfo, myMesh );

}

// Run the simulation.
void SEMproxy::run()
{
  time_point< system_clock > startComputeTime, startOutputTime, totalComputeTime, totalOutputTime;

  for( int indexTimeSample=0; indexTimeSample<myInfo.myNumSamples; indexTimeSample++ )
  {
    startComputeTime = system_clock::now();
    mySolver.computeOneStep( indexTimeSample, myInfo.myOrderNumber_a, myInfo.myOrderNumber_b, myInfo.myOrderNumber_c, myInfo.nPointsPerElement, i1, i2, myInfo, myRHSTerm, pnGlobal, rhsElement );
    totalComputeTime += system_clock::now() - startComputeTime;

    startOutputTime = system_clock::now();
    mySolver.outputPnValues( myMesh, indexTimeSample, i1, myInfo.myElementSource, pnGlobal );
    swap( i1, i2 );
    totalOutputTime += system_clock::now() - startOutputTime;
  }
  float kerneltime_ms = time_point_cast< microseconds >( totalComputeTime ).time_since_epoch().count();
  float outputtime_ms = time_point_cast< microseconds >( totalOutputTime ).time_since_epoch().count();

  cout << "------------------------------------------------ "<< endl;
  cout << "\n---- Elapsed Kernel Time : "<< kerneltime_ms/1E6<<" seconds."<< endl;
  cout << "---- Elapsed Output Time : "<< outputtime_ms/1E6<<" seconds."<< endl;
  cout << "------------------------------------------------ "<< endl;

}

// Initialize arrays
void SEMproxy::getMeshInfo()
{
  // get information from mesh
  myInfo.numberOfNodes=myMesh.getNumberOfNodes();
  myInfo.numberOfElements=myMesh.getNumberOfElements();
  myInfo.numberOfPointsPerElement=myMesh.getNumberOfPointsPerElement();
  myInfo.numberOfInteriorNodes=myMesh.getNumberOfInteriorNodes();
  myInfo.numberOfPointsPerElement=myMesh.getNumberOfPointsPerElement();
  myInfo.numberOfBoundaryNodes=myMesh.getNumberOfBoundaryNodes();
  myInfo.numberOfBoundaryFaces=myMesh.getNumberOfBoundaryFaces();
  #ifdef SEM_MESHCOLOR
  myInfo.numberMaxOfElementsByColor=myMesh.getNumberOfElementsByColor();
  #endif
  
}

// Initialize arrays
void SEMproxy::init_arrays()
{
  cout<<"Allocate host memory for source and pressure values ..."<< endl;
  myRHSTerm=allocateArray2D< arrayReal >( myInfo.myNumberOfRHS, myInfo.myNumSamples, "RHSTerm" );
  rhsElement=allocateVector< vectorInt >( myInfo.myNumberOfRHS, "rhsElement" );
  pnGlobal=allocateArray2D< arrayReal >( myInfo.numberOfNodes, 2, "pnGlobal" );
}

// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation=allocateArray2D< arrayReal >( myInfo.myNumberOfRHS, 3, "RHSLocation" );
  // set number of rhs and location
  myRHSLocation( 0, 0 )=1001;
  myRHSLocation( 0, 1 )=1001;
  myRHSLocation( 0, 2 )=1001;
  cout << "\nSource location: "<<myRHSLocation( 0, 0 )<<", "<<myRHSLocation( 0, 1 )<<", "<<myRHSLocation( 0, 2 )<< endl;
  for( int i=0; i<myInfo.myNumberOfRHS; i++ )
  {
    //extract element number for current rhs
    rhsElement[i]=myMesh.getElementNumberFromPoints( myRHSLocation( i, 0 ), myRHSLocation( i, 1 ), myRHSLocation( i, 2 ));
    //printf("Element number for the source %d is: %d\n", i, rhsElement[i]);
  }

  // initialize source term
  vector< float > sourceTerm=myUtils.computeSourceTerm( myInfo.myNumSamples, myInfo.myTimeStep, myInfo.f0, myInfo.sourceOrder );
  for( int j=0; j<myInfo.myNumSamples; j++ )
  {
    myRHSTerm( 0, j )=sourceTerm[j];
    if( j%100==0 )
      cout<<"Sample "<<j<<"\t: sourceTerm = "<<sourceTerm[j]<< endl;
  }
  // get element number of source term
  myInfo.myElementSource=rhsElement[0];
  cout <<"Element number for the source location: "<<myInfo.myElementSource<< endl<< endl;
}
