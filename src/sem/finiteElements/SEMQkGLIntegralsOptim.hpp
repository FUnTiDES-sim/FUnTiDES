#ifndef SEMQKGLINTEGRALSOPTIM_HPP_
#define SEMQKGLINTEGRALSOPTIM_HPP_

#include "dataType.hpp"
#include "SEMmacros.hpp"
#include "SEMdata.hpp"
#include "SEMQkGLBasisFunctions.hpp"
using namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
class SEMQkGLIntegralsOptim
{
private:
  int order;
  struct SEMinfo infos;
  SEMQkGLBasisFunctions GLBasis;

  ////////////////////////////////////////////////////////////////////////////////////
  //  from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
  constexpr static double sqrt5 = 2.2360679774997897;
  // order of polynomial approximation
  //constexpr static  int r=2;
  constexpr static  int r=SEMinfo::myOrderNumber;
  // number of support/quadrature/nodes points in one direction
  constexpr static int numSupport1dPoints=r+1;
  constexpr static int num1dNodes=numSupport1dPoints;
  // Half the number of support points, rounded down. Precomputed for efficiency
  constexpr static int halfNodes = ( numSupport1dPoints - 1 )/ 2;
  // the number of nodes/support points per element
  constexpr static int numSupportPoints=(r+1)*(r+1)*(r+1);

public:
  PROXY_HOST_DEVICE SEMQkGLIntegralsOptim(){};
  PROXY_HOST_DEVICE ~SEMQkGLIntegralsOptim(){};
  
  /////////////////////////////////////////////////////////////////////////////////////
  //  from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Calculates the linear index for support/quadrature points from ijk
   *   coordinates.
   * @param r order of polynomial approximation
   * @param i The index in the xi0 direction (0,r)
   * @param j The index in the xi1 direction (0,r)
   * @param k The index in the xi2 direction (0,r)
   * @return The linear index of the support/quadrature point (0-(r+1)^3)
  */
  PROXY_HOST_DEVICE 
  constexpr static int linearIndex( const int r,
                                    const int i,
                                    const int j,
                                    const int k ) 
  {
           return i + (r+1) * j + (r+1)*(r+1) * k;
  }

  /**
   * @brief Calculate the Cartesian/TensorProduct index given the linear index
   *   of a support point.
   * @param linearIndex The linear index of support point
   * @param r order of polynomial approximation
   * @param i0 The Cartesian index of the support point in the xi0 direction.
   * @param i1 The Cartesian index of the support point in the xi1 direction.
   * @param i2 The Cartesian index of the support point in the xi2 direction.
  */
  PROXY_HOST_DEVICE 
  constexpr static void multiIndex( int const r,int const linearIndex, int & i0, int & i1, int & i2 ) 
  {
        i2 = linearIndex/((r+1)*(r+1));
        i1 = (linearIndex%((r+1)*(r+1)))/(r+1);
        i0 = (linearIndex%((r+1)*(r+1)))%(r+1);
  }

  /**
   * @brief Compute the interpolation coefficients of the q-th quadrature point in a given direction
   * @param q the index of the quadrature point in 1D
   * @param k the index of the interval endpoint (0 or 1)
   * @return The interpolation coefficient
  */
  PROXY_HOST_DEVICE 
  constexpr static double interpolationCoord( const int order, const int q, const int k ) 
  {
      const double alpha = (SEMQkGLBasisFunctions::parentSupportCoord<SEMinfo>( q ) + 1.0 ) / 2.0;
      return k == 0 ? ( 1.0 - alpha ) : alpha;
   }

  /**
   * @brief Compute the 1D factor of the coefficient of the jacobian on the q-th quadrature point,
   * with respect to the k-th interval endpoint (0 or 1). The computation depends on the position
   * in the basis tensor product of this term (i, equal to 0, 1 or 2) and on the direction in which
   * the gradient is being computed (dir, from 0 to 2)
   * @param q The index of the quadrature point in 1D
   * @param i The index of the position in the tensor product
   * @param k The index of the interval endpoint (0 or 1)
   * @param dir The direction in which the derivatives are being computed
   * @return The value of the jacobian factor
  */
  PROXY_HOST_DEVICE 
  constexpr static double jacobianCoefficient1D( const int order, const int q, const int i, const int k, const int dir )
  {
      if( i == dir )
      {
        return k== 0 ? -1.0/2.0 : 1.0/2.0;
      }
      else
      {
        return interpolationCoord( order, q, k );
      }
  }

  
  PROXY_HOST_DEVICE
  double determinant(double  m[3][3]) const
  {
     return abs(m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
               -m[0][1]*(m[1][0]*m[2][2]-m[2][0]*m[1][2])
               +m[0][2]*(m[1][0]*m[2][1]-m[2][0]*m[1][1]));
  }
   
  /**
     * @brief Calculates the isoparametric "Jacobian" transformation
     *  matrix/mapping from the parent space to the physical space.
     * @param qa The 1d quadrature point index in xi0 direction (0,1)
     * @param qb The 1d quadrature point index in xi1 direction (0,1)
     * @param qc The 1d quadrature point index in xi2 direction (0,1)
     * @param X Array containing the coordinates of the mesh support points.
     * @param J Array to store the Jacobian transformation.
  */
   PROXY_HOST_DEVICE
   void jacobianTransformation( int e, int const r,
                                int const qa,
                                int const qb,
                                int const qc,
                                double const (&X)[8][3],
                                double ( & J )[3][3] ) const
   {
      for( int k = 0; k < 8; k++ )
      {
        const int ka = k % 2;
        const int kb = ( k % 4 ) / 2;
        const int kc = k / 4;
        for( int j = 0; j < 3; j++ )
        {
          double jacCoeff = jacobianCoefficient1D(r, qa, 0, ka, j ) *
                            jacobianCoefficient1D(r, qb, 1, kb, j ) *
                            jacobianCoefficient1D(r, qc, 2, kc, j );
          for( int i = 0; i < 3; i++ )
          {
            J[i][j] +=  jacCoeff * X[k][i];
          }
        }
      }
   }

   /**
    * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
    *   mass matrix M, i.e., the superposition matrix of the shape functions.
    * @param q The quadrature point index
    * @param X Array containing the coordinates of the mesh support points.
    * @return The diagonal mass term associated to q
   */
   PROXY_HOST_DEVICE
   double computeMassTerm( int e, int const r, int const q, double const (&X)[8][3] ) const
   {
      int qa, qb, qc;
      multiIndex( r,q, qa, qb, qc );
      const double w3D = SEMQkGLBasisFunctions::weight<SEMinfo>( qa )*SEMQkGLBasisFunctions::weight<SEMinfo>( qb )*SEMQkGLBasisFunctions::weight<SEMinfo>( qc );
      double J[3][3] = {{0}};
      jacobianTransformation(e, r, qa, qb, qc, X, J );
      return determinant( J )*w3D;
   }

  /**
   * @brief Invert the symmetric matrix @p srcSymMatrix and store the result in @p dstSymMatrix.
   * @param dstSymMatrix The 3x3 symmetric matrix to write the inverse to.
   * @param srcSymMatrix The 3x3 symmetric matrix to take the inverse of.
   * @return The determinant.
   * @note @p srcSymMatrix can contain integers but @p dstMatrix must contain floating point values.
  */
  PROXY_HOST_DEVICE 
  void symInvert( double  dstSymMatrix[6], double  srcSymMatrix[6]) const
  {
   
     using FloatingPoint = std::decay_t< decltype( dstSymMatrix[ 0 ] ) >;
   
     dstSymMatrix[ 0 ] = srcSymMatrix[ 1 ] * srcSymMatrix[ 2 ] - srcSymMatrix[ 3 ] * srcSymMatrix[ 3 ];
     dstSymMatrix[ 5 ] = srcSymMatrix[ 4 ] * srcSymMatrix[ 3 ] - srcSymMatrix[ 5 ] * srcSymMatrix[ 2 ];
     dstSymMatrix[ 4 ] = srcSymMatrix[ 5 ] * srcSymMatrix[ 3 ] - srcSymMatrix[ 4 ] * srcSymMatrix[ 1 ];
   
     double det = srcSymMatrix[ 0 ] * dstSymMatrix[ 0 ] + srcSymMatrix[ 5 ] * dstSymMatrix[ 5 ] + srcSymMatrix[ 4 ] * dstSymMatrix[ 4 ];
  
     FloatingPoint const invDet = FloatingPoint( 1 ) / det;
   
     dstSymMatrix[ 0 ] *= invDet;
     dstSymMatrix[ 5 ] *= invDet;
     dstSymMatrix[ 4 ] *= invDet;
     dstSymMatrix[ 1 ] = ( srcSymMatrix[ 0 ] * srcSymMatrix[ 2 ] - srcSymMatrix[ 4 ] * srcSymMatrix[ 4 ] ) * invDet;
     dstSymMatrix[ 3 ] = ( srcSymMatrix[ 5 ] * srcSymMatrix[ 4 ] - srcSymMatrix[ 0 ] * srcSymMatrix[ 3 ] ) * invDet;
     dstSymMatrix[ 2 ] = ( srcSymMatrix[ 0 ] * srcSymMatrix[ 1 ] - srcSymMatrix[ 5 ] * srcSymMatrix[ 5 ] ) * invDet;
   
  }
  
  /**
   * @brief Invert the symmetric matrix @p symMatrix overwritting it.
   * @param symMatrix The 3x3 symmetric matrix to take the inverse of and overwrite.
   * @return The determinant.
   * @note @p symMatrix can contain integers but @p dstMatrix must contain floating point values.
  */
  PROXY_HOST_DEVICE  
  void symInvert0( double  symMatrix[6] ) const
  {
      std::remove_reference_t< decltype( symMatrix[ 0 ] ) > temp[ 6 ];
      symInvert( temp, symMatrix );
      
      symMatrix[0]=temp[0];
      symMatrix[1]=temp[1];
      symMatrix[2]=temp[2];
      symMatrix[3]=temp[3];
      symMatrix[4]=temp[4];
      symMatrix[5]=temp[5];
  }
  
  /**
   * @brief Calculates the isoparametric "geometrical" transformation
   *  matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   * @param B Array to store the  the geometrical symetic matrix=detJ*J^{-1}J^{-T}.
  */
  PROXY_HOST_DEVICE 
  void computeBMatrix( int e, int const r,  int const qa, int const qb, int const qc,
                       double const (&X)[8][3],
                       double (& J)[3][3],
                       double (& B)[6] ) const
  {
      jacobianTransformation( e, r, qa, qb, qc, X, J );
      double const detJ = determinant( J );

  
      // compute J^T.J/det(J), using Voigt notation for B
      B[0] = (J[0][0]*J[0][0]+J[1][0]*J[1][0]+J[2][0]*J[2][0])/detJ;
      B[1] = (J[0][1]*J[0][1]+J[1][1]*J[1][1]+J[2][1]*J[2][1])/detJ;
      B[2] = (J[0][2]*J[0][2]+J[1][2]*J[1][2]+J[2][2]*J[2][2])/detJ;
      B[3] = (J[0][1]*J[0][2]+J[1][1]*J[1][2]+J[2][1]*J[2][2])/detJ;
      B[4] = (J[0][0]*J[0][2]+J[1][0]*J[1][2]+J[2][0]*J[2][2])/detJ;
      B[5] = (J[0][0]*J[0][1]+J[1][0]*J[1][1]+J[2][0]*J[2][1])/detJ;
  
      // compute detJ*J^{-1}J^{-T}
      symInvert0( B );
  }
  
  template<typename FUNC>
  PROXY_HOST_DEVICE
  void computeGradPhiBGradPhi( int const e,
                               int const r,
                               int const qa,
                               int const qb,
                               int const qc,
                               double const (&B)[6],
                               FUNC && func ) const
  {
     //const double w = GLBasis.weight<SEMinfo>(qa )*GLBasis.weight<SEMinfo>(qb )*GLBasis.weight<SEMinfo>(qc );
     const double w = SEMQkGLBasisFunctions::weight<SEMinfo>(qa )*
	              SEMQkGLBasisFunctions::weight<SEMinfo>(qb )*
		      SEMQkGLBasisFunctions::weight<SEMinfo>(qc );
     for( int i=0; i<num1dNodes; i++ )
     {
       const int ibc = linearIndex( r,i, qb, qc );
       const int aic = linearIndex( r,qa, i, qc );
       const int abi = linearIndex( r,qa, qb, i );
       const double gia = SEMQkGLBasisFunctions::basisGradientAt(r, i, qa );
       const double gib = SEMQkGLBasisFunctions::basisGradientAt(r, i, qb );
       const double gic = SEMQkGLBasisFunctions::basisGradientAt(r, i, qc );
       for( int j=0; j<num1dNodes; j++ )
       {
         const int jbc = linearIndex( r,j, qb, qc );
         const int ajc = linearIndex( r,qa, j, qc );
         const int abj = linearIndex( r,qa, qb, j );
         const double gja = SEMQkGLBasisFunctions::basisGradientAt(r, j, qa );
         const double gjb = SEMQkGLBasisFunctions::basisGradientAt(r, j, qb );
         const double gjc = SEMQkGLBasisFunctions::basisGradientAt(r, j, qc );
         // diagonal terms
         const double w0 = w * gia * gja;
         func( ibc, jbc, w0 * B[0] );
         const double w1 = w * gib * gjb;
         func( aic, ajc, w1 * B[1] );
         const double w2 = w * gic * gjc;
         func( abi, abj, w2 * B[2] );
         // off-diagonal terms
         const double w3 = w * gib * gjc;
         func( aic, abj, w3 * B[3] );
         func( abj, aic, w3 * B[3] );
         const double w4 = w * gia * gjc;
         func( ibc, abj, w4 * B[4] );
         func( abj, ibc, w4 * B[4] );
         const double w5 = w * gia * gjb;
         func( ibc, ajc, w5 * B[5] );
         func( ajc, ibc, w5 * B[5] );
       }
     }
  }
  
  
  template<typename FUNC>
  PROXY_HOST_DEVICE
  void computeStiffnessTerm( int e,
                             int r,
                             int const q,
                             double const (&X)[8][3],
                             FUNC && func ) const
  {
        int qa, qb, qc;
        multiIndex( r,q, qa, qb, qc );
        double B[6] = {0};
        double J[3][3] = {{0}};
        computeBMatrix( e, r, qa, qb, qc, X, J, B );
        computeGradPhiBGradPhi( e, r, qa, qb, qc, B, func );
  }
  
  
  /**
   * @brief compute  mass Matrix stiffnessVector.
   */
  PROXY_HOST_DEVICE 
  void computeMassMatrixAndStiffnessVector(const int & elementNumber,
                                           const int & order,
                                           const int & nPointsPerElement,
                                           ARRAY_REAL_VIEW const & nodesCoordsX,
                                           ARRAY_REAL_VIEW const & nodesCoordsY,
                                           ARRAY_REAL_VIEW const & nodesCoordsZ,
                                           float massMatrixLocal[],
                                           float pnLocal[],
                                           float Y[]) const
  {
      double X[8][3];
      int I=0;
      for( int k=0;k<order+1;k+=order )
      {
          for ( int j=0; j<order+1;j+=order )
          {
              for( int i=0;i<order+1;i+=order )
              {
                  int l=i+j*(order+1)+k*(order+1)*(order+1);
                  X[I][0]=nodesCoordsX(elementNumber,l);
                  X[I][1]=nodesCoordsZ(elementNumber,l);
                  X[I][2]=nodesCoordsY(elementNumber,l);
                  I++;
              }
          }
      }
      for (int q=0;q<nPointsPerElement;q++)
      {
         Y[q]=0;
      }
      for (int q=0;q<nPointsPerElement;q++)
      {
          massMatrixLocal[q]=computeMassTerm( elementNumber, order,q, X); 
          computeStiffnessTerm(elementNumber, order, q, X, [&] (const int i, const int j, const double val)
                  {
                   Y[i] = Y[i] + val*pnLocal[j];
                  });
      }
  }
  /////////////////////////////////////////////////////////////////////////////////////
  //  end from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
  
};
  
#endif //SEMQKGLINTEGRALSOPTIM_HPP_
