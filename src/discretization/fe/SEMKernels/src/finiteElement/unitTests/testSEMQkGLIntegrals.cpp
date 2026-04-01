
#include "../SEMQkGLIntegralsClassic.hpp"
#include "../SEMQkGLIntegralsOptim.hpp"
#include "../SEMQkGLIntegralsShiva.hpp"

#include <gtest/gtest.h>

using namespace shiva;
using namespace shiva::functions;
using namespace shiva::geometry;
using namespace shiva::discretizations::finiteElementMethod;

void setX( double (&X)[8][3] )
{
  double x0 = 0, y0 = 0, z0 = 0;
  double x1 = 20, y1 = 20, z1 = 20;

  double X0[8][3] = { {x0, y0, z0},
    {x1, y0, z0},
    {x0, y1, z0},
    {x1, y1, z0},
    {x0, y0, z1},
    {x1, y0, z1},
    {x0, y1, z1},
    {x1, y1, z1} };

  for ( int a = 0; a < 8; ++a )
  {
    for ( int i = 0; i < 3; ++i )
    {
      X[a][i] = X0[a][i] ;//* ( 0.9 + 0.1 * (rand() % 20) );
    }
  }
}

template< typename T >
void setXYZ( T & X, T & Y, T & Z )
{
  double x0 = -1.0, y0 = -1.0, z0 = -1.0;
  double x1 =  1.0, y1 =  1.0, z1 =  1.0;

  // double x0 = -1.1, y0 = -1.2, z0 = -1.2;
  // double x1 =  1.2, y1 =  1.1, z1 =  1.4;

  X( 0, 0 ) = x0; Y( 0, 0 ) = y0; Z( 0, 0 ) = z0;
  X( 0, 1 ) = x1; Y( 0, 1 ) = y0; Z( 0, 1 ) = z0;
  X( 0, 2 ) = x0; Y( 0, 2 ) = y1; Z( 0, 2 ) = z0;
  X( 0, 3 ) = x1; Y( 0, 3 ) = y1; Z( 0, 3 ) = z0;
  X( 0, 4 ) = x0; Y( 0, 4 ) = y0; Z( 0, 4 ) = z1;
  X( 0, 5 ) = x1; Y( 0, 5 ) = y0; Z( 0, 5 ) = z1;
  X( 0, 6 ) = x0; Y( 0, 6 ) = y1; Z( 0, 6 ) = z1;
  X( 0, 7 ) = x1; Y( 0, 7 ) = y1; Z( 0, 7 ) = z1;
  
  srand(123);
  for ( int a = 0; a < 8; ++a )
  {
    X( 0, a ) = X( 0, a ) * ( 0.9 + 0.01 * (rand() % 21) );
    Y( 0, a ) = Y( 0, a ) * ( 0.9 + 0.01 * (rand() % 21) );
    Z( 0, a ) = Z( 0, a ) * ( 0.9 + 0.01 * (rand() % 21) );
  }
}



template< int ORDER >
struct MassAndStiffnessSolutions;

template<>
struct MassAndStiffnessSolutions<1>
{ 
  constexpr static int length = 8;
  static constexpr float mass[length] = { 9.74300026893616e-01, 9.87788021564484e-01, 1.00481104850769e+00, 1.03112494945526e+00, 9.57599997520447e-01, 1.00162827968597e+00, 9.57065880298615e-01, 1.01409268379211e+00 };
  static constexpr float pN[length] = { 1.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00 };
  static constexpr float Y[length] = { 1.63464689254761e+00, 5.55009663105011e-01, 4.78901594877243e-01, -8.43870174139738e-03, 4.89260762929916e-01, -3.51113788783550e-02, -6.79248124361038e-02, 0.00000000000000e+00 };
};


template< typename INTEGRALS >
void computeMassMatrixAndStiffnessVectorTester()
{
  using Integrals = INTEGRALS;
  static constexpr int order = INTEGRALS::order;


  CArrayNd<double, 1, 8> Xcoords;
  CArrayNd<double, 1, 8> Ycoords;
  CArrayNd<double, 1, 8> Zcoords;
  setXYZ( Xcoords, Ycoords, Zcoords );

  constexpr int massOffset = 0;
  constexpr int pOffset = (order+1)*(order+1)*(order+1);
  constexpr int YOffset = 2 * (order+1)*(order+1)*(order+1);
  constexpr int length = (order+1)*(order+1)*(order+1);
  
  constexpr int dataSize = 3 * (order+1)*(order+1)*(order+1);
  float * hostData = new float[ dataSize ];

  hostData[pOffset] = 1.0;

  pmpl::genericKernelWrapper( dataSize, 
                              hostData, 
                              [ Xcoords,
                                Ycoords, 
                                Zcoords,
                                massOffset,
                                pOffset,
                                YOffset,
                                length ] SHIVA_DEVICE ( auto * device_data )
  {
    float * const massMatrixLocal = device_data + massOffset;
    float * const pnLocal = device_data + pOffset;
    float * const Y = device_data + YOffset;

    Integrals integralInterface;

    integralInterface.init();
    integralInterface.computeMassMatrixAndStiffnessVector( 0,
                                                    length,
                                                    Xcoords,
                                                    Ycoords,
                                                    Zcoords,
                                                    massMatrixLocal,
                                                    pnLocal,
                                                    Y );  
  } );


  float maxMassSoln = 0.0;
  float maxYSoln = 0.0;
  for( int i = 0; i < length; ++i )
  {
    maxMassSoln = std::max( maxMassSoln, std::abs( MassAndStiffnessSolutions<order>::mass[i] ) );
    maxYSoln = std::max( maxYSoln, std::abs( MassAndStiffnessSolutions<order>::Y[i] ) );
  }


  printf( "YSoln = { " );
  for( int i=0; i < length; ++i )
  {
    printf( "%18.14e, ", hostData[ YOffset + i ] );
  }
  printf( "};\n" );

  for( int i = 0; i < length; ++i )
  {
     EXPECT_NEAR( MassAndStiffnessSolutions<order>::mass[i], hostData[massOffset + i], maxMassSoln * 1e-4 );
  //   EXPECT_NEAR( MassAndStiffnessSolutions<order>::Y[i], hostData[YOffset + i], maxYSoln * 1e-4 );
  }
}

TEST( testSEMQkGLIntegralsShiva, computeMassMatrixAndStiffnessVector )
{
  using TransformType =
    LinearTransform< double,
                     InterpolatedShape< double,
                                        Cube< double >,
                                        LagrangeBasis< double, 1, EqualSpacing >,
                                        LagrangeBasis< double, 1, EqualSpacing >,
                                        LagrangeBasis< double, 1, EqualSpacing > > >;

  constexpr int order = 1;
  using ParentElementType =
    ParentElement< double,
                   Cube< double >,
                   LagrangeBasis< double, order, EqualSpacing >,
                   LagrangeBasis< double, order, EqualSpacing >,
                   LagrangeBasis< double, order, EqualSpacing > >;

  using Integrals = SEMQkGLIntegralsShiva< double, order, TransformType, ParentElementType >;

  computeMassMatrixAndStiffnessVectorTester< Integrals >();
}


TEST( testSEMQkGLIntegralsOptim, computeMassMatrixAndStiffnessVector )
{
  constexpr int order = 1;
  using Integrals = SEMQkGLIntegralsOptim<order>;

  computeMassMatrixAndStiffnessVectorTester< Integrals >();
}

TEST( testSEMQkGLIntegralsClassic, computeMassMatrixAndStiffnessVector )
{
  constexpr int order = 1;
  using Integrals = SEMQkGLIntegralsClassic<order>;

  computeMassMatrixAndStiffnessVectorTester< Integrals >();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
