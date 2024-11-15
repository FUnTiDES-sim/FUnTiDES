#ifndef SEMQKGLINTEGRALS_HPP_
#define SEMQKGLINTEGRALS_HPP_

#include "dataType.hpp"
#include "SEMmacros.hpp"
#include "SEMdata.hpp"
#include "SEMQkGLBasisFunctions.hpp"
using namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
class SEMQkGLIntegrals
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
  PROXY_HOST_DEVICE SEMQkGLIntegrals(){};
  PROXY_HOST_DEVICE ~SEMQkGLIntegrals(){};



  ////////////////////////////////////////////////////////////////////////////////////
  //  from first implementation
  /////////////////////////////////////////////////////////////////////////////////////
  // V1
  // compute B and M
  PROXY_HOST_DEVICE void computeB( const int & elementNumber,
                                     const int & order,
                                     VECTOR_DOUBLE_VIEW const & weights,
                                     ARRAY_INT_VIEW const & nodesList,
                                     ARRAY_REAL_VIEW const & nodesCoords,
                                     ARRAY_DOUBLE_VIEW const & dPhi,
                                     float massMatrixLocal[],
                                     float B[][COL] ) const
  {
      #ifdef SEM2D
      {
        for( int i2=0; i2<order+1; i2++ )
        {
          for( int i1=0; i1<order+1; i1++ )
          {
            // compute jacobian matrix
            double jac0=0;
            double jac1=0;
            double jac2=0;
            double jac3=0;
            int i=i1+i2*(order+1);
            for( int j1=0; j1<order+1; j1++ )
            {
              int j=j1+i2*(order+1);
              int localToGlobal=nodesList( elementNumber, j );
              double X=nodesCoords( localToGlobal, 0 );
              double Y=nodesCoords( localToGlobal, 1 );
              jac0+=X*dPhi( j1, i1 );
              jac2+=Y*dPhi( j1, i1 );
            }
            for( int j2=0; j2<order+1; j2++ )
            {
              int j=i1+j2*(order+1);
              int localToGlobal=nodesList( elementNumber, j );
              double X=nodesCoords( localToGlobal, 0 );
              double Y=nodesCoords( localToGlobal, 1 );
              jac1+=X*dPhi( j2, i2 );
              jac3+=Y*dPhi( j2, i2 );
            }
            // detJ
            double detJ=abs( jac0*jac3-jac2*jac1 );
            double invJac0=jac3;
            double invJac1=-jac1;
            double invJac2=-jac2;
            double invJac3=jac0;
            double transpInvJac0=jac3;
            double transpInvJac1=-jac2;
            double transpInvJac2=-jac1;
            double transpInvJac3=jac0;
            double detJM1=1./detJ;
            // B
            B[i][0]=(invJac0*transpInvJac0+invJac1*transpInvJac2)*detJM1;
            B[i][1]=(invJac0*transpInvJac1+invJac1*transpInvJac3)*detJM1;
            B[i][2]=(invJac2*transpInvJac0+invJac3*transpInvJac2)*detJM1;
            B[i][3]=(invJac2*transpInvJac1+invJac3*transpInvJac3)*detJM1;
            //M
            massMatrixLocal[i]=weights[i1]*weights[i2]*detJ;
          }
        }
      }
      #else   //3D case
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
                int localToGlobal=nodesList( elementNumber, j );
                double X=nodesCoords( localToGlobal, 0 );
                double Y=nodesCoords( localToGlobal, 2 );
                double Z=nodesCoords( localToGlobal, 1 );
                jac00+=X*dPhi( j1, i1 );
                jac10+=Z*dPhi( j1, i1 );
                jac20+=Y*dPhi( j1, i1 );
              }
              for( int j2=0; j2<order+1; j2++ )
              {
                int j=i1+j2*(order+1)+i3*(order+1)*(order+1);
                int localToGlobal=nodesList( elementNumber, j );
                double X=nodesCoords( localToGlobal, 0 );
                double Y=nodesCoords( localToGlobal, 2 );
                double Z=nodesCoords( localToGlobal, 1 );
                jac01+=X*dPhi( j2, i2 );
                jac11+=Z*dPhi( j2, i2 );
                jac21+=Y*dPhi( j2, i2 );
              }
              for( int j3=0; j3<order+1; j3++ )
              {
                int j=i1+i2*(order+1)+j3*(order+1)*(order+1);
                int localToGlobal=nodesList( elementNumber, j );
                double X=nodesCoords( localToGlobal, 0 );
                double Y=nodesCoords( localToGlobal, 2 );
                double Z=nodesCoords( localToGlobal, 1 );
                jac02+=X*dPhi( j3, i3 );
                jac12+=Z*dPhi( j3, i3 );
                jac22+=Y*dPhi( j3, i3 );
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
      #endif
  }
    
  // V2
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
      #ifdef SEM2D
      {
        for( int i2=0; i2<order+1; i2++ )
        {
          for( int i1=0; i1<order+1; i1++ )
          {
            // compute jacobian matrix
            double jac0=0;
            double jac1=0;
            double jac2=0;
            double jac3=0;
            int i=i1+i2*(order+1);
            for( int j1=0; j1<order+1; j1++ )
            {
              int j=j1+i2*(order+1);
              double X=nodesCoordsX( elementNumber, j );
              double Z=nodesCoordsZ( elementNumber, j );
              jac0+=X*dPhi( j1, i1 );
              jac2+=Z*dPhi( j1, i1 );
            }
            for( int j2=0; j2<order+1; j2++ )
            {
              int j=i1+j2*(order+1);
              double X=nodesCoordsX( elementNumber, j );
              double Z=nodesCoordsZ( elementNumber, j );
              jac1+=X*dPhi( j2, i2 );
              jac3+=Z*dPhi( j2, i2 );
            }
            // detJ
            double detJ=abs( jac0*jac3-jac2*jac1 );
            double invJac0=jac3;
            double invJac1=-jac1;
            double invJac2=-jac2;
            double invJac3=jac0;
            double transpInvJac0=jac3;
            double transpInvJac1=-jac2;
            double transpInvJac2=-jac1;
            double transpInvJac3=jac0;
            double detJM1=1./detJ;
            // B
            B[i][0]=(invJac0*transpInvJac0+invJac1*transpInvJac2)*detJM1;
            B[i][1]=(invJac0*transpInvJac1+invJac1*transpInvJac3)*detJM1;
            B[i][2]=(invJac2*transpInvJac0+invJac3*transpInvJac2)*detJM1;
            B[i][3]=(invJac2*transpInvJac1+invJac3*transpInvJac3)*detJM1;
            //M
            massMatrixLocal[i]=weights[i1]*weights[i2]*detJ;
          }
        }
      }
      #else   //3D case
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
      #endif
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
      #ifdef SEM2D
      {
        // B11
        for( int i2=0; i2<order+1; i2++ )
        {
          for( int i1=0; i1<order+1; i1++ )
          {
            for( int j=0; j<nPointsPerElement; j++ )
            {
              R[j]=0;
            }
            for( int j1=0; j1<order+1; j1++ )
            {
              int j=j1+i2*(order+1);
              for( int m=0; m<order+1; m++ )
              {
                R[j]+=weights[m]*weights[i2]*(B[m+i2*(order+1)][0]*dPhi( i1, m )*dPhi( j1, m ));
              }
            }
            // B21
            for( int j1=0; j1<order+1; j1++ )
            {
              for( int j2=0; j2<order+1; j2++ )
              {
                int j=j1+j2*(order+1);
                R[j]+=weights[i1]*weights[j2]*(B[i1+j2*(order+1)][1]*dPhi( i2, j2 )*dPhi( j1, i1 ));
              }
            }
            // B12
            for( int j1=0; j1<order+1; j1++ )
            {
              for( int j2=0; j2<order+1; j2++ )
              {
                int j=j1+j2*(order+1);
                R[j]+=weights[i2]*weights[j1]*(B[i2+j1*(order+1)][2]*dPhi( i1, j1 )*dPhi( j2, i2 ));
              }
            }
            // B22
            for( int j2=0; j2<order+1; j2++ )
            {
              int j=i1+j2*(order+1);
              for( int n=0; n<order+1; n++ )
              {
                R[j]+=weights[i1]*weights[n]*(B[i1+n*(order+1)][3]*dPhi( i2, n )*dPhi( j2, n ));
              }
            }
            int i=i1+i2*(order+1);
            Y[i]=0;
            for( int j=0; j<nPointsPerElement; j++ )
            {
              Y[i]+=R[j]*pnLocal[j];
            }
          }
        }
      }
      #else
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
      #endif
  }
  
  // V1
  // compute stiffnessVector.
  // returns mass matrix and stiffness vector local to an element
  PROXY_HOST_DEVICE void computeMassMatrixAndStiffnessVector(const int & elementNumber,
                                                                      const int & order,
                                                                      const int & nPointsPerElement,
                                                                      ARRAY_INT_VIEW const & nodesList,
                                                                      ARRAY_REAL_VIEW const & nodesCoords,
                                                                      VECTOR_DOUBLE_VIEW const & weights,
                                                                      ARRAY_DOUBLE_VIEW const & dPhi,
                                                                      float massMatrixLocal[],
                                                                      float const pnLocal[],
                                                                      float Y[]) const
  {
      float B[ROW][COL];
      float R[ROW];
      // compute Jacobian, massMatrix and B
       computeB( elementNumber, order, weights, nodesList, nodesCoords, dPhi, massMatrixLocal, B );
       // compute stifness  matrix ( durufle's optimization)
       gradPhiGradPhi( nPointsPerElement, order, weights, dPhi, B, pnLocal, R, Y );
  }
  
  // V2
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
  
  //computeDs
  PROXY_HOST_DEVICE void computeDs( const int & iFace,
                                     const int & order,
                                     ARRAY_INT_VIEW const & faceInfos,
                                     int numOfBasisFunctionOnFace,
                                     float Js[][6],
                                     ARRAY_REAL_VIEW const & globalNodesCoords,
                                     ARRAY_DOUBLE_VIEW const & dPhi,
                                     float ds[] ) const
  {
      int face=faceInfos( iFace, 1 );
      for( int j=0; j<order+1; j++ )
      {
        Js[0][j]=0;    // x
        Js[1][j]=0;    // y
        for( int i=0; i<order+1; i++ )
        {
          float xi=globalNodesCoords( faceInfos( iFace, 2+i ), 0 );
          float yi=globalNodesCoords( faceInfos( iFace, 2+i ), 1 );
          if( face==0 || face==2 )
          {
            Js[0][j]+=dPhi( i, j )*xi;
            Js[1][j]+=dPhi( i, j )*yi;
          }
          if( face==1 || face==3 )
          {
            Js[0][j]+=dPhi( i, j )*xi;
            Js[1][j]+=dPhi( i, j )*yi;
          }
        }
        ds[j]=sqrt( Js[0][j]*Js[0][j]+Js[1][j]*Js[1][j] );
      }
  }
  /////////////////////////////////////////////////////////////////////////////////////
  //  end from first implementation
  /////////////////////////////////////////////////////////////////////////////////////
  
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
      //const double w3D = GLBasis.weight<SEMinfo>( qa )*GLBasis.weight<SEMinfo>( qb )*GLBasis.weight<SEMinfo>( qc );
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
     const double w = SEMQkGLBasisFunctions::weight<SEMinfo>(qa )*SEMQkGLBasisFunctions::weight<SEMinfo>(qb )*SEMQkGLBasisFunctions::weight<SEMinfo>(qc );
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
          computeStiffnessTerm(elementNumber, order, q, X,[&] (const int i, const int j, const double val)
                  {
                   float localIncrement=val*pnLocal[j];
                   Y[i]+=localIncrement;
                  });
      }
  }
  /////////////////////////////////////////////////////////////////////////////////////
  //  end from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
  
};
  
#endif //SEMQKGLINTEGRALS_HPP_
