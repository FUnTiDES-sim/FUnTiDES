#ifndef SEMQKGLINTEGRALSCLASSIC_HPP_
#define SEMQKGLINTEGRALSCLASSIC_HPP_

#include "dataType.hpp"
#include "SEMmacros.hpp"
#include "SEMdata.hpp"
#include "SEMQkGLBasisFunctions.hpp"
using namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
class SEMQkGLIntegralsClassic
{
private:
  int order;
  struct SEMinfo infos;
  SEMQkGLBasisFunctions GLBasis;

public:
  PROXY_HOST_DEVICE SEMQkGLIntegralsClassic(){};
  PROXY_HOST_DEVICE ~SEMQkGLIntegralsClassic(){};

  // compute B and M
  PROXY_HOST_DEVICE void computeB( const int & elementNumber,
                                     const int & order,
                                     VECTOR_DOUBLE_VIEW const & weights,
                                     ARRAY_REAL_VIEW const & nodesCoordsX,
                                     ARRAY_REAL_VIEW const & nodesCoordsY,
                                     ARRAY_REAL_VIEW const & nodesCoordsZ,
                                     ARRAY_DOUBLE_VIEW const & dPhi,
                                     float massMatrixLocal[],
                                     float B[][COL] ) const
  {
      for( int i3=0; i3<order+1; i3++ )
      {
        for( int i2=0; i2<order+1; i2++ )
        {
          for( int i1=0; i1<order+1; i1++ )
          {
            int i=i1+i2*(order+1)+i3*(order+1)*(order+1);
            // compute jacobian matrix
            double jac00=0;
            double jac01=0;
            double jac02=0;
            double jac10=0;
            double jac11=0;
            double jac12=0;
            double jac20=0;
            double jac21=0;
            double jac22=0;
   
            for( int j1=0; j1<order+1; j1++ )
            {
              int j=j1+i2*(order+1)+i3*(order+1)*(order+1);
              double X=nodesCoordsX( elementNumber, j );
              double Y=nodesCoordsY( elementNumber, j );
              double Z=nodesCoordsZ( elementNumber, j );
              jac00+=X*dPhi( j1, i1 );
              jac20+=Y*dPhi( j1, i1 );
              jac10+=Z*dPhi( j1, i1 );
              }
            for( int j2=0; j2<order+1; j2++ )
            {
              int j=i1+j2*(order+1)+i3*(order+1)*(order+1);
              double X=nodesCoordsX( elementNumber, j );
              double Y=nodesCoordsY( elementNumber, j );
              double Z=nodesCoordsZ( elementNumber, j );
              jac01+=X*dPhi( j2, i2 );
              jac21+=Y*dPhi( j2, i2 );
              jac11+=Z*dPhi( j2, i2 );
            }
            for( int j3=0; j3<order+1; j3++ )
            {
              int j=i1+i2*(order+1)+j3*(order+1)*(order+1);
              double X=nodesCoordsX( elementNumber, j );
              double Y=nodesCoordsY( elementNumber, j );
              double Z=nodesCoordsZ( elementNumber, j );
              jac02+=X*dPhi( j3, i3 );
              jac22+=Y*dPhi( j3, i3 );
              jac12+=Z*dPhi( j3, i3 );
            }
            // detJ
            double detJ=abs( jac00*(jac11*jac22-jac21*jac12)
                             -jac01*(jac10*jac22-jac20*jac12)
                             +jac02*(jac10*jac21-jac20*jac11));
   
            // inv of jac is equal of the minors of the transposed of jac
            double invJac00=jac11*jac22-jac12*jac21;
            double invJac01=jac02*jac21-jac01*jac22;
            double invJac02=jac01*jac12-jac02*jac11;
            double invJac10=jac12*jac20-jac10*jac22;
            double invJac11=jac00*jac22-jac02*jac20;
            double invJac12=jac02*jac10-jac00*jac12;
            double invJac20=jac10*jac21-jac11*jac20;
            double invJac21=jac01*jac20-jac00*jac21;
            double invJac22=jac00*jac11-jac01*jac10;
    
            double transpInvJac00=invJac00;
            double transpInvJac01=invJac10;
            double transpInvJac02=invJac20;
            double transpInvJac10=invJac01;
            double transpInvJac11=invJac11;
            double transpInvJac12=invJac21;
            double transpInvJac20=invJac02;
            double transpInvJac21=invJac12;
            double transpInvJac22=invJac22;
    
            double detJM1=1./detJ;
    
            // B
            B[i][0]=(invJac00*transpInvJac00+invJac01*transpInvJac10+invJac02*transpInvJac20)*detJM1;    //B11
            B[i][1]=(invJac10*transpInvJac01+invJac11*transpInvJac11+invJac12*transpInvJac21)*detJM1;    //B22
            B[i][2]=(invJac20*transpInvJac02+invJac21*transpInvJac12+invJac22*transpInvJac22)*detJM1;    //B33
            B[i][3]=(invJac00*transpInvJac01+invJac01*transpInvJac11+invJac02*transpInvJac21)*detJM1;    //B12,B21
            B[i][4]=(invJac00*transpInvJac02+invJac01*transpInvJac12+invJac02*transpInvJac22)*detJM1;    //B13,B31
            B[i][5]=(invJac10*transpInvJac02+invJac11*transpInvJac12+invJac12*transpInvJac22)*detJM1;    //B23,B32
    
            //M
            massMatrixLocal[i]=weights[i1]*weights[i2]*weights[i3]*detJ;
  
           }
         }
       }
  }
    
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  PROXY_HOST_DEVICE void gradPhiGradPhi( const int & nPointsPerElement,
                                           const int & order,
                                           VECTOR_DOUBLE_VIEW const & weights,
                                           ARRAY_DOUBLE_VIEW const & dPhi,
                                           float const B[][COL],
                                           float const pnLocal[],
                                           float R[],
                                           float Y[] ) const
  {
      int orderPow2=(order+1)*(order+1);
      for( int i3=0; i3<order+1; i3++ )
      {
        for( int i2=0; i2<order+1; i2++ )
        {
          for( int i1=0; i1<order+1; i1++ )
          {
            for( int j=0; j<nPointsPerElement; j++ )
            {
              R[j]=0;
            }
   
            //B11
            for( int j1=0; j1<order+1; j1++ )
            {
              int j=j1+i2*(order+1)+i3*orderPow2;
              for( int l=0; l<order+1; l++ )
              {
                int ll=l+i2*(order+1)+i3*orderPow2;
                R[j]+=weights[l]*weights[i2]*weights[i3]*(B[ll][0]*dPhi( i1, l )*dPhi( j1, l ));
              }
            }
            //B22
            for( int j2=0; j2<order+1; j2++ )
            {
              int j=i1+j2*(order+1)+i3*orderPow2;
              for( int m=0; m<order+1; m++ )
              {
                int mm=i1+m*(order+1)+i3*orderPow2;
                R[j]+=weights[i1]*weights[m]*weights[i3]*(B[mm][1]*dPhi( i2, m )*dPhi( j2, m ));
              }
            }
            //B33
            for( int j3=0; j3<order+1; j3++ )
            {
              int j=i1+i2*(order+1)+j3*orderPow2;
              for( int n=0; n<order+1; n++ )
              {
                int nn=i1+i2*(order+1)+n*orderPow2;
                R[j]+=weights[i1]*weights[i2]*weights[n]*(B[nn][2]*dPhi( i3, n )*dPhi( j3, n ));
              }
            }
            // B12,B21 (B[][3])
            for( int j2=0; j2<order+1; j2++ )
            {
              for( int j1=0; j1<order+1; j1++ )
              {
                int j=j1+j2*(order+1)+i3*orderPow2;
                int k=j1+i2*(order+1)+i3*orderPow2;
                int l=i1+j2*(order+1)+i3*orderPow2;
                R[j]+=weights[j1]*weights[i2]*weights[i3]*(B[k][3]*dPhi( i1, j1 )*dPhi( j2, i2 ))+
                       weights[i1]*weights[j2]*weights[i3]*(B[l][3]*dPhi( j1, i1 )*dPhi( i2, j2 ));
              }
            }
            // B13,B31 (B[][4])
            for( int j3=0; j3<order+1; j3++ )
            {
              for( int j1=0; j1<order+1; j1++ )
              {
                int j=j1+i2*(order+1)+i3*orderPow2;
                int k=j1+i2*(order+1)+i3*orderPow2;
                int l=j1+i2*(order+1)+j3*orderPow2;
                R[j]+=weights[j1]*weights[i2]*weights[i3]*(B[k][4]*dPhi( j1, i1 )*dPhi( j3, i3 ))+
                       weights[j1]*weights[i2]*weights[j3]*(B[l][4]*dPhi( j1, i1 )*dPhi( i3, j3 ));
              }
            }
            // B23,B32 (B[][5])
            for( int j3=0; j3<order+1; j3++ )
            {
              for( int j2=0; j2<order+1; j2++ )
              {
                int j=i1+j2*(order+1)+j3*orderPow2;
                int k=i1+j2*(order+1)+i3*orderPow2;
                int l=i1+i2*(order+1)+j3*orderPow2;
                R[j]+=weights[i1]*weights[j2]*weights[i3]*(B[k][5]*dPhi( i2, i2 )*dPhi( j3, i3 ))+
                       weights[i1]*weights[i2]*weights[j3]*(B[l][5]*dPhi( j2, i2 )*dPhi( i3, j3 ));
              }
            }
    
            int i=i1+i2*(order+1)+i3*orderPow2;
            Y[i]=0;
            for( int j=0; j<nPointsPerElement; j++ )
            {
              Y[i]+=R[j]*pnLocal[j];
            }
  
          }
        }
      }
  }
  
  // compute stiffnessVector.
  // returns mass matrix and stiffness vector local to an element
  PROXY_HOST_DEVICE void computeMassMatrixAndStiffnessVector(const int & elementNumber,
                                                                      const int & order,
                                                                      const int & nPointsPerElement,
                                                                      ARRAY_REAL_VIEW const & nodesCoordsX,
                                                                      ARRAY_REAL_VIEW const & nodesCoordsY,
                                                                      ARRAY_REAL_VIEW const & nodesCoordsZ,
                                                                      VECTOR_DOUBLE_VIEW const & weights,
                                                                      ARRAY_DOUBLE_VIEW const & dPhi,
                                                                      float massMatrixLocal[],
                                                                      float const pnLocal[],
                                                                      float Y[]) const
  {
      float B[ROW][COL];
      float R[ROW];
      // compute Jacobian, massMatrix and B
      computeB( elementNumber, order, weights, nodesCoordsX,nodesCoordsY, nodesCoordsZ, dPhi, massMatrixLocal, B );
      // compute stifness  matrix ( durufle's optimization)
      gradPhiGradPhi( nPointsPerElement, order, weights, dPhi, B, pnLocal, R, Y );
  }
  
  /////////////////////////////////////////////////////////////////////////////////////
  //  end from first implementation
  /////////////////////////////////////////////////////////////////////////////////////
  
};
  
#endif //SEMQKGLINTEGRALSCLASSIC_HPP_
