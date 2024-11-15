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
  
};
  
#endif //SEMQKGLINTEGRALSCLASSIC_HPP_
