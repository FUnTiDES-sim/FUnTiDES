//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <chrono>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>
#include <omp.h>

//#define MemSpace Kokkos::CudaUVMSpace
#define MemSpace Kokkos::SharedSpace
//#define MemSpace Kokkos::HostSpace


using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy< ExecSpace >;

/*
   typedef Kokkos::View<int*, Kokkos::LayoutRight, MemSpace>   vectorInt;
   typedef Kokkos::View<int**, Kokkos::LayoutRight, MemSpace>   arrayInt;
   typedef Kokkos::View<float**, Kokkos::LayoutRight, MemSpace>   arrayReal;
   typedef Kokkos::View<double**, Kokkos::LayoutRight, MemSpace>   arrayDouble;
 */
typedef Kokkos::View< int *, MemSpace >   vectorInt;
typedef Kokkos::View< int * *, MemSpace >   arrayInt;
typedef Kokkos::View< float * *, MemSpace >   arrayReal;
typedef Kokkos::View< double * *, MemSpace >   arrayDouble;
typedef Kokkos::View< double * * *, MemSpace >   array3DDouble;

template< class T >
T allocateVector( int n1 )
{
  std::cout<<"allocate vector of size "<<n1<<std::endl;
  T vect( "v", n1 );
  return vect;
}
template< class T >
T allocateArray2D( int n1, int n2 )
{
  std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
  T array( "a", n1, n2 );
  return array;
}


void globalNodesList( const int & order, const int & ex, const int & ey, const int & nx, arrayInt const & nodesList )
{
  for( int j=0; j<ey; j++ )
  {
    for( int i=0; i<ex; i++ )
    {
      int n0=i+j*ex;
      int offset=i*order+j*order*nx;
      for( int k=0; k<order+1; k++ )
      {
        for( int l=0; l<order+1; l++ )
        {
          int dofLocal=l+k*(order+1);
          int dofGlobal=offset+l+k*nx;
          nodesList( n0, dofLocal )=dofGlobal;
        }
      }
    }
  }
}

void nodesCoordinates( const int & nx, const int & ny, const int & order, const int & ex, const int & ey,
                       const float & hx, const float & hy, const int & numberOfNodes, arrayReal const & nodeCoords )
{
  //arrayReal nodeCoords( numberOfNodes, 2 );
  std::vector< float > coordX( nx );
  std::vector< float > coordY( ny );
  std::vector< float > xi( order+1 );

  switch( order )
  {
    case 1:
      xi[0]=-1.;
      xi[1]=1.;
      break;
    case 2:
      xi[0]=-1.;
      xi[1]=0.;
      xi[2]=1.;
      break;
    case 3:
      static constexpr double sqrt5 = 2.2360679774997897;
      xi[0] = -1.0;
      xi[1] = -1./sqrt5;
      xi[2] = 1./sqrt5;
      xi[3] = 1.;
      break;
    case 4:
      static constexpr double sqrt3_7 = 0.6546536707079771;
      xi[0] = -1.0;
      xi[1] = -sqrt3_7;
      xi[2] = 0.0;
      xi[3] = sqrt3_7;
      xi[4] = 1.0;
      break;
    case 5:
      static constexpr double sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
      static constexpr double sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
      static constexpr double sqrt_inv21 = 0.218217890235992381;
      xi[0] = -1.0;
      xi[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
      xi[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
      xi[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
      xi[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
      xi[5] = 1.0;
      break;
    default:
      break;
  }
  for( int i=0; i<ex; i++ )
  {
    float x0=i*hx;
    float x1=(i+1)*hx;
    float b=(x1+x0)/2.;
    float a=b-x0;

    for( int j=0; j<order+1; j++ )
    {
      coordX[j+i*order]=a*xi[j]+b;
    }
  }
  for( int i=0; i<ey; i++ )
  {
    float y0=i*hy;
    float y1=(i+1)*hy;
    float b=(y1+y0)/2.;
    float a=b-y0;
    for( int j=0; j<order+1; j++ )
    {
      coordY[j+i*order]=a*xi[j]+b;
    }
  }
  for( int j=0; j<ny; j++ )
  {
    for( int i=0; i<nx; i++ )
    {
      nodeCoords( i+nx*j, 0 )=coordX[i];
      nodeCoords( i+nx*j, 1 )=coordY[j];
    }
  }
}

KOKKOS_FUNCTION int localToGlobalNodes( const int & elementNumber,
                                        const int & nPointsPerElement,
                                        arrayInt const & nodesList,
                                        int localToGlobal[] )
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    localToGlobal[i]=nodesList( elementNumber, i );
  }
  return 0;
}

KOKKOS_FUNCTION int getXi( const int & numberOfPointsPerElement,
                           arrayReal const & globalNodesCoords,
                           int const localToGlobal[],
                           double Xi[][2] )
{
  for( int i=0; i<numberOfPointsPerElement; i++ )
  {
    Xi[i][0]=globalNodesCoords( localToGlobal[i], 0 );
    Xi[i][1]=globalNodesCoords( localToGlobal[i], 1 );
  }
  return 0;
}

void sortElementsByColor( const int & ex, const int & ey, arrayInt const & listOfElementsByColor )
{
  // red
  int k=0;
  for( int j=0; j<ey; j=j+2 )
  {
    for( int i=0; i<ex; i=i+2 )
    {
      k=k+1;
      listOfElementsByColor( 0, k )=i+j*ex;
    }
  }
  listOfElementsByColor( 0, 0 )=k;
  // green
  k=0;
  for( int j=0; j<ey; j=j+2 )
  {
    for( int i=1; i<ex; i=i+2 )
    {
      k=k+1;
      listOfElementsByColor( 1, k )=i+j*ex;
    }
  }
  listOfElementsByColor( 1, 0 )=k;
  // blue
  k=0;
  for( int j=1; j<ey; j=j+2 )
  {
    for( int i=0; i<ex; i=i+2 )
    {
      k=k+1;
      listOfElementsByColor( 2, k )=i+j*ex;
    }
  }
  listOfElementsByColor( 2, 0 )=k;
  // yellow
  k=0;
  for( int j=1; j<ey; j=j+2 )
  {
    for( int i=1; i<ex; i=i+2 )
    {
      k=k+1;
      listOfElementsByColor( 3, k )=i+j*ex;
    }
  }
  listOfElementsByColor( 3, 0 )=k;
}


int main( int argc, char *argv[] )
{
  Kokkos::initialize( argc, argv );
  {
    int order=1;
    int ex=100;
    int ey=100;
    int nx=ex*order+1;
    int ny=ey*order+1;
    float hx=10;
    float hy=10;
    int numberOfElements=ex*ey;
    int numberOfPointsPerElement=(order+1)*(order+1);
    int numberOfNodes=(ex*order+1)*(ey*order+1);

    int numberOfElementsPerColor=(ex/2+ex%2)*(ey/2+ey%2);
    int nbColors=4;
    arrayInt listOfElementsByColor( "elementByColor", nbColors, numberOfElementsPerColor );
    printf ( "maximum number of mesh per color %d\n", numberOfElementsPerColor );
    sortElementsByColor( ex, ey, listOfElementsByColor );
    for( int e=0; e<10; e++ )
    {
      printf( "color red element    %d global element %d\n", e, listOfElementsByColor( 0, e+1 ));
      printf( "color green element  %d global element %d\n", e, listOfElementsByColor( 1, e+1 ));
      printf( "color blue element   %d global element %d\n", e, listOfElementsByColor( 2, e+1 ));
      printf( "color yellow element %d global element %d\n", e, listOfElementsByColor( 3, e+1 ));
    }
    printf( "ici \n" );
    int nthreads=omp_get_max_threads();
    printf( "nthreads %d\n", nthreads );
    arrayInt nodesList( "nodeList", numberOfElements, numberOfPointsPerElement );
    arrayReal nodeCoords( "nodeCoords", numberOfNodes, 2 );
    globalNodesList( order, ex, ey, nx, nodesList );
    nodesCoordinates( nx, ny, order, ex, ey, hx, hy, numberOfNodes, nodeCoords );


// we expect :
// 0 0 0 10 0 0 10 10 10
// 1 10 0 20 0 10 10 20 10
// 2 20 0 30 0 20 10 30 10
// 3 30 0 40 0 30 10 40 10
// 4 40 0 50 0 40 10 50 10
// 5 50 0 60 0 50 10 60 10
// 6 60 0 70 0 60 10 70 10
// 7 70 0 80 0 70 10 80 10
// 8 80 0 90 0 80 10 90 10
// 9 90 0 100 0 90 10 100 10
//
// in this case results are correct.
    for( int color=0; color<nbColors; color++ )
    {
      numberOfElements=listOfElementsByColor( color, 0 );
      printf( " color %d numberOf Elements %d\n", color, numberOfElements );
      Kokkos::parallel_for( numberOfElements, KOKKOS_LAMBDA ( int eColor )
      {
        int localToGlobal[36];
        double Xi[36][2];
        int e=listOfElementsByColor( color, eColor+1 );
        //printf("eColor %d e %d\n",eColor,e);
        int i=localToGlobalNodes( e, numberOfPointsPerElement, nodesList, localToGlobal );
        //get global coordinates Xi of element e
        int j=getXi( numberOfPointsPerElement, nodeCoords, localToGlobal, Xi );

        if( e<10 )
          printf( " element %d (%f %f) (%f %f) (%f %f) (%f %f)   \n ", e, Xi[0][0], Xi[0][1], Xi[1][0], Xi[1][1], Xi[2][0], Xi[2][1], Xi[3][0], Xi[3][1] );
        if( e>100 && e<111 )
          printf( " element %d (%f %f) (%f %f) (%f %f) (%f %f)   \n ", e, Xi[0][0], Xi[0][1], Xi[1][0], Xi[1][1], Xi[2][0], Xi[2][1], Xi[3][0], Xi[3][1] );
      } );
      Kokkos::fence();
    }
  }
  Kokkos::finalize();
  return (0);
}
