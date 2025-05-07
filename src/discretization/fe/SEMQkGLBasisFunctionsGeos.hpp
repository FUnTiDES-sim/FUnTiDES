#ifndef SEMQKGLBASISFUNCTIONSGEOS_HPP_
#define SEMQKGLBASISFUNCTIONSGEOS_HPP_
#include <dataType.hpp>
#include <SEMmacros.hpp>
#include <SEMdata.hpp>

#include "LagrangeBasis1.hpp"
#include "LagrangeBasis2.hpp"
#include "LagrangeBasis3GL.hpp"
#include "LagrangeBasis4GL.hpp"
#include "LagrangeBasis5GL.hpp"

using namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
class SEMQkGLBasisFunctionsGeos
{
private:
  int order;
  struct SEMinfo infos;

  ////////////////////////////////////////////////////////////////////////////////////
  //  from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
  constexpr static double sqrt5 = 2.2360679774997897;
  // order of polynomial approximation
  constexpr static  int r=SEMinfo::myOrderNumber;
  // number of support/quadrature/nodes points in one direction
  constexpr static int numSupport1dPoints=r+1;
  constexpr static int num1dNodes=numSupport1dPoints;
  // Half the number of support points, rounded down. Precomputed for efficiency
  constexpr static int halfNodes = ( numSupport1dPoints - 1 )/ 2;
  // the number of nodes/support points per element
  constexpr static int numSupportPoints=(r+1)*(r+1)*(r+1);

public:
  PROXY_HOST_DEVICE SEMQkGLBasisFunctionsGeos(){};
  PROXY_HOST_DEVICE ~SEMQkGLBasisFunctionsGeos(){};

  
  template<typename SEMinfo>
  PROXY_HOST_DEVICE
  constexpr static double parentSupportCoord(  const int supportPointIndex ) 
  {
      double result=0.0;
      switch( SEMinfo::myOrderNumber )
      {
        case 1:
          return -1.0 + 2.0 * (supportPointIndex & 1);
        case 2:
          switch( supportPointIndex )
          {
            case 0:
              return -1.0;
              break;
            case 2:
              return 1.0;
            case 1:
            default:
              return 0.0;
          }
        case 3:
          switch( supportPointIndex )
          {
            case 0:
               result = -1.0;
               break;
            case 1:
              result = -1.0/sqrt5;
              break;
            case 2:
              result = 1.0/sqrt5;
              break;
            case 3:
              result = 1.0;
              break;
            default:
              break;
          }
        default:
           return 0;
      }
      return result;
   }

  /**
   * @brief The gradient of the basis function for a support point evaluated at
   *  a given support point. By symmetry, p is assumed to be in 0, ..., (N-1)/2
   * @param q The index of the basis function
   * @param p The index of the support point
   * @return The gradient of basis function.
  */
  template<typename SEMinfo>
  PROXY_HOST_DEVICE  
  constexpr static double gradientAt( const int q, const int p ) 
  {
      switch( SEMinfo::myOrderNumber )
      {
        case 1:
          return q == 0 ? -0.5 : 0.5;
        case 2:
          switch( q )
          {
            case 0:
              return p == 0 ? -1.5 : -0.5;
            case 1:
              return p == 0 ? 2.0 : 0.0;
            case 2:
              return p == 0 ? -0.5 : 0.5;
            default:
              return 0;
          }
        case 3:
          switch( q )
          {
            case 0:
              return p == 0 ? -3.0 : -0.80901699437494742410;
            case 1:
              return p == 0 ? 4.0450849718747371205 : 0.0;
            case 2:
              return p == 0 ? -1.5450849718747371205 : 1.1180339887498948482;
            case 3:
              return p == 0 ? 0.5 : -0.30901699437494742410;
            default:
              return 0;
          }
        default:
           return 0;
      }
   }

  /*
   * @brief Compute the 1st derivative of the q-th 1D basis function at quadrature point p
   * @param q the index of the 1D basis funcion
   * @param p the index of the 1D quadrature point
   * @return The derivative value
  */
  PROXY_HOST_DEVICE 
  constexpr static double basisGradientAt(  const int q, const int p )
  {
     if( p <= halfNodes )
     {
       return gradientAt<SEMinfo>( q, p );
     }
     else
     {
       return -gradientAt<SEMinfo>( numSupport1dPoints - 1 - q, numSupport1dPoints - 1 - p );
     }
  }

  /**
   * @brief The value of the weight for the given support point
   * @param q The index of the support point
   * @return The value of the weight
  */
  template<typename SEMinfo>
  PROXY_HOST_DEVICE
  constexpr static double weight( const int q )
  {
      switch(SEMinfo::myOrderNumber)
      {
        case 1:
           return 1;
        case 2:
           switch( q )
           {
             case 0:
             case 2:
               return 1.0/3.0;
             default:
               return 4.0/3.0;
           }
        case 3:
           switch( q )
           {
             case 1:
             case 2:
               return 5.0/6.0;
             default:
              return 1.0/6.0;
           }
       default:
           return 0;
      }
   }

};
  
#endif //SEMQKGLBASISFUNCTIONS_HPP_
