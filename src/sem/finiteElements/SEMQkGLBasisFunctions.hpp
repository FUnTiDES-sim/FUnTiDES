#ifndef SEMQKGLBASISFUNCTIONS_HPP_
#define SEMQKGLBASISFUNCTIONS_HPP_

#include "dataType.hpp"
#include "SEMmacros.hpp"
#include "SEMdata.hpp"
using namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
class SEMQkGLBasisFunctions
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
  PROXY_HOST_DEVICE SEMQkGLBasisFunctions(){};
  PROXY_HOST_DEVICE ~SEMQkGLBasisFunctions(){};

  static  void gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) //const
  {
   switch( order )
   {
     case 1:
       quadraturePoints[0]=-1.;
       quadraturePoints[1]=1.;
       break;
     case 2:
       quadraturePoints[0]=-1.;
       quadraturePoints[1]=0.;
       quadraturePoints[2]=1.;
       break;
     case 3:
       static constexpr double sqrt5 = 2.2360679774997897;
       quadraturePoints[0] = -1.0;
       quadraturePoints[1] = -1./sqrt5;
       quadraturePoints[2] = 1./sqrt5;
       quadraturePoints[3] = 1.;
       break;
     case 4:
       static constexpr double sqrt3_7 = 0.6546536707079771;
       quadraturePoints[0] = -1.0;
       quadraturePoints[1] = -sqrt3_7;
       quadraturePoints[2] = 0.0;
       quadraturePoints[3] = sqrt3_7;
       quadraturePoints[4] = 1.0;
       break;
     case 5:
       static constexpr double sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
       static constexpr double sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
       static constexpr double sqrt_inv21 = 0.218217890235992381;
       quadraturePoints[0] = -1.0;
       quadraturePoints[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
       quadraturePoints[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
       quadraturePoints[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
       quadraturePoints[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
       quadraturePoints[5] = 1.0;
       break;
     default:
       break;
   }
  }

  void gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const
  {
    if( order == 1 )
    {
      weights[0]=1.0;
      weights[1]=1.0;
    }
    if( order == 2 )
    {
      weights[0]=0.33333333;
      weights[1]=1.33333333;
      weights[2]= 0.33333333;
    }
    if( order == 3 )
    {
      weights[0]=0.16666667;
      weights[1]=0.83333333;
      weights[2]=0.83333333;
      weights[3]=0.16666667;
    }
    if( order == 4 )
    {
      weights[0]=0.1;
      weights[1]=0.54444444;
      weights[2]=0.71111111;
      weights[3]=0.54444444;
      weights[4]=0.1;
    }
    if( order == 5 )
    {
      weights[0]=0.06666667;
      weights[1]=0.37847496;
      weights[2]=0.55485838;
      weights[3]=0.55485838;
      weights[4]=0.37847496;
      weights[5]=0.06666667;
    }
  }
  vector< double > shapeFunction1D( int order, double xi ) const
  {
    std::vector< double > shapeFunction( order+1 );
    if( order==1 )
    {
      shapeFunction[0]=0.5*(1.0-xi);
      shapeFunction[1]=0.5*(1.0+xi);
    }
    if( order==2 )
    {
      shapeFunction[0]= -1.0*xi*(0.5 - 0.5*xi);
      shapeFunction[1]=(1.0 - 1.0*xi)*(1.0*xi + 1.0);
      shapeFunction[2]= 1.0*xi*(0.5*xi + 0.5);
    }
    if( order==3 )
    {
      shapeFunction[0]=(0.309016994374947 - 0.690983005625053*xi)*(0.5 - 0.5*xi)
                        *(-1.80901699437495*xi - 0.809016994374947);
      shapeFunction[1]=(0.5 - 1.11803398874989*xi)*(0.690983005625053 - 0.690983005625053*xi)
                        *(1.80901699437495*xi + 1.80901699437495);
      shapeFunction[2]=(1.80901699437495 - 1.80901699437495*xi)
                        *(0.690983005625053*xi + 0.690983005625053)*(1.11803398874989*xi + 0.5);
      shapeFunction[3]=(0.5*xi + 0.5)*(0.690983005625053*xi + 0.309016994374947)
                        *(1.80901699437495*xi - 0.809016994374947);
    }
    if( order==4 )
    {
      shapeFunction[0]=1.0*xi*(0.39564392373896 - 0.60435607626104*xi)*(0.5 - 0.5*xi)
                        *(-2.89564392373896*xi - 1.89564392373896);
      shapeFunction[1]=-1.52752523165195*xi*(0.5 - 0.763762615825973*xi)*(0.60435607626104 - 0.60435607626104*xi)
                        *(2.89564392373896*xi + 2.89564392373896);
      shapeFunction[2]=(1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0);
      shapeFunction[3]= 1.52752523165195*xi*(2.89564392373896 - 2.89564392373896*xi)
                        *(0.60435607626104*xi + 0.60435607626104)*(0.763762615825973*xi + 0.5);
      shapeFunction[4]= 1.0*xi*(0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896)
                      *(2.89564392373896*xi - 1.89564392373896);
    }
    if( order==5 )
    {
      shapeFunction[0]=(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                        *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)
                        *(-1.39905441140358*xi - 0.399054411403579);
      shapeFunction[1]=(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                        *(0.566554479308753 - 0.566554479308753*xi)*(-2.0840983387567*xi - 0.594450529658367)
                        *(4.25632117622354*xi + 4.25632117622354);
      shapeFunction[2]=(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                        *(0.778069933064125 - 0.778069933064125*xi)*(1.39905441140358*xi + 1.39905441140358)
                        *(2.0840983387567*xi + 1.59445052965837);
      shapeFunction[3]=(1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                        *(0.778069933064125*xi + 0.778069933064125)*(0.952120850728289*xi + 0.728425125873928)
                        *(1.75296196636787*xi + 0.5);
      shapeFunction[4]=(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi + 0.566554479308753)
                        *(0.6535475074298*xi + 0.5)*(0.952120850728289*xi + 0.271574874126072)
                        *(2.0840983387567*xi - 0.594450529658367);
      shapeFunction[5]=(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                        *(0.778069933064125*xi + 0.221930066935875)*(1.39905441140358*xi - 0.399054411403579)
                        *(4.25632117622354*xi - 3.25632117622354);
    }
    return shapeFunction;
  }
  vector< double > derivativeShapeFunction1D( int order, double xi ) const
  {
    std::vector< double > derivativeShapeFunction( order+1 );

    if( order == 1 )
    {
      derivativeShapeFunction[0]=-0.5;
      derivativeShapeFunction[1]=0.5;
    }
    if( order == 2 )
    {
      derivativeShapeFunction[0]=1.0*xi - 0.5;
      derivativeShapeFunction[1]=-2.0*xi;
      derivativeShapeFunction[2]=1.0*xi + 0.5;
    }
    if( order == 3 )
    {
      derivativeShapeFunction[0]=-1.80901699437495*(0.309016994374947 - 0.690983005625053*xi)*(0.5 - 0.5*xi)
                                  + (-1.80901699437495*xi - 0.809016994374947)*(0.345491502812526*xi - 0.345491502812526)
                                  + (-1.80901699437495*xi - 0.809016994374947)*(0.345491502812526*xi - 0.154508497187474);

      derivativeShapeFunction[1]=1.80901699437495*(0.5 - 1.11803398874989*xi)*(0.690983005625053 - 0.690983005625053*xi)
                                  + (0.772542485937369*xi - 0.772542485937369)*(1.80901699437495*xi + 1.80901699437495)
                                  + (0.772542485937369*xi - 0.345491502812526)*(1.80901699437495*xi + 1.80901699437495);

      derivativeShapeFunction[2]=(1.80901699437495 - 1.80901699437495*xi)*(0.772542485937369*xi + 0.345491502812526) +
                                  (1.80901699437495 - 1.80901699437495*xi)*(0.772542485937369*xi + 0.772542485937369) -
                                  1.80901699437495*(0.690983005625053*xi + 0.690983005625053)*(1.11803398874989*xi + 0.5);

      derivativeShapeFunction[3]=(0.345491502812526*xi + 0.154508497187474)*(1.80901699437495*xi - 0.809016994374947) +
                                  (0.345491502812526*xi + 0.345491502812526)*(1.80901699437495*xi - 0.809016994374947) +
                                  1.80901699437495*(0.5*xi + 0.5)*(0.690983005625053*xi + 0.309016994374947);
    }
    if( order == 4 )
    {
      derivativeShapeFunction[0]=2.89564392373896*xi*(0.39564392373896 - 0.60435607626104*xi)*(0.5 - 0.5*xi) +
                                  0.5*xi*(0.39564392373896 - 0.60435607626104*xi)*(-2.89564392373896*xi - 1.89564392373896)
                                  + 0.60435607626104*xi*(0.5 - 0.5*xi)*(-2.89564392373896*xi - 1.89564392373896) +
                                  (0.39564392373896 - 0.60435607626104*xi)*(-2.89564392373896*xi - 1.89564392373896)*(0.5*xi - 0.5);

      derivativeShapeFunction[1]=-4.42316915539091*xi*(0.5 - 0.763762615825973*xi)*(0.60435607626104 - 0.60435607626104*xi)
                                  + 0.923169155390906*xi*(0.5 - 0.763762615825973*xi)*(2.89564392373896*xi + 2.89564392373896)
                                  + 1.16666666666667*xi*(0.60435607626104 - 0.60435607626104*xi)*(2.89564392373896*xi + 2.89564392373896)
                                  + (0.60435607626104 - 0.60435607626104*xi)*(1.16666666666667*xi - 0.763762615825973)
                                  *(2.89564392373896*xi + 2.89564392373896);

      derivativeShapeFunction[2]=(1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.52752523165195*xi + 1.0) +
                                  (1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.52752523165195*xi + 1.52752523165195)
                                  - 1.0*(1.0 - 1.52752523165195*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0)
                                  - 1.52752523165195*(1.0 - 1.0*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0);

      derivativeShapeFunction[3]=1.16666666666667*xi*(2.89564392373896 - 2.89564392373896*xi)*(0.60435607626104*xi + 0.60435607626104)
                                  + 0.923169155390906*xi*(2.89564392373896 - 2.89564392373896*xi)*(0.763762615825973*xi + 0.5)
                                  - 4.42316915539091*xi*(0.60435607626104*xi + 0.60435607626104)*(0.763762615825973*xi + 0.5)
                                  + (2.89564392373896 - 2.89564392373896*xi)*(0.60435607626104*xi + 0.60435607626104)
                                  *(1.16666666666667*xi + 0.763762615825973);

      derivativeShapeFunction[4]=2.89564392373896*xi*(0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896)
                                  + 0.60435607626104*xi*(0.5*xi + 0.5)*(2.89564392373896*xi - 1.89564392373896)
                                  + 0.5*xi*(0.60435607626104*xi + 0.39564392373896)*(2.89564392373896*xi - 1.89564392373896)
                                  + (0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896)*(2.89564392373896*xi - 1.89564392373896);
    }
    if( order == 5 )
    {
      derivativeShapeFunction[0]=-1.39905441140358*(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                                  *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)
                                  - 4.25632117622354*(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                                  *(0.5 - 0.5*xi)*(-1.39905441140358*xi - 0.399054411403579) + (0.221930066935875 - 0.778069933064125*xi)
                                  *(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579)
                                  *(0.283277239654376*xi - 0.283277239654376) + (0.221930066935875 - 0.778069933064125*xi)
                                  *(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579)
                                  *(0.283277239654376*xi - 0.216722760345624) - 0.778069933064125*(0.433445520691247 - 0.566554479308753*xi)
                                  *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579);

      derivativeShapeFunction[1]= -2.0840983387567*(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                                  *(0.566554479308753 - 0.566554479308753*xi)*(4.25632117622354*xi + 4.25632117622354)
                                  - 0.566554479308753*(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                                  *(-2.0840983387567*xi - 0.594450529658367)*(4.25632117622354*xi + 4.25632117622354)
                                  +(0.271574874126072 - 0.952120850728289*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                  *(2.12816058811177 - 2.78170809554157*xi)*(-2.0840983387567*xi - 0.594450529658367)
                                  +(0.271574874126072 - 0.952120850728289*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                  *(-2.78170809554157*xi - 2.78170809554157)*(-2.0840983387567*xi - 0.594450529658367)
                                  - 0.952120850728289*(0.5 - 0.6535475074298*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                  *(-2.0840983387567*xi - 0.594450529658367)*(4.25632117622354*xi + 4.25632117622354);

      derivativeShapeFunction[2]= 2.0840983387567*(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                                  *(0.778069933064125 - 0.778069933064125*xi)*(1.39905441140358*xi + 1.39905441140358)
                                  + 1.39905441140358*(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                                  *(0.778069933064125 - 0.778069933064125*xi)*(2.0840983387567*xi + 1.59445052965837)
                                  - 0.952120850728289*(0.5 - 1.75296196636787*xi)*(0.778069933064125 - 0.778069933064125*xi)
                                  *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837)
                                  + (0.728425125873928 - 0.952120850728289*xi)*(1.3639269998358*xi - 1.3639269998358)
                                  *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837)
                                  + (0.728425125873928 - 0.952120850728289*xi)*(1.3639269998358*xi - 0.389034966532063)
                                  *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837);

      derivativeShapeFunction[3]=0.952120850728289*(1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                  *(0.778069933064125*xi + 0.778069933064125)*(1.75296196636787*xi + 0.5)
                                  + (1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                  *(0.952120850728289*xi + 0.728425125873928)*(1.3639269998358*xi + 0.389034966532063)
                                  + (1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                  *(0.952120850728289*xi + 0.728425125873928)*(1.3639269998358*xi + 1.3639269998358)
                                  - 2.0840983387567*(1.39905441140358 - 1.39905441140358*xi)*(0.778069933064125*xi + 0.778069933064125)
                                  *(0.952120850728289*xi + 0.728425125873928)*(1.75296196636787*xi + 0.5)
                                  - 1.39905441140358*(1.59445052965837 - 2.0840983387567*xi)*(0.778069933064125*xi + 0.778069933064125)
                                  *(0.952120850728289*xi + 0.728425125873928)*(1.75296196636787*xi + 0.5);

      derivativeShapeFunction[4]=(2.78170809554157 - 2.78170809554157*xi)*(0.566554479308753*xi + 0.566554479308753)
                                  *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367)
                                  + 2.0840983387567*(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi + 0.566554479308753)
                                  *(0.6535475074298*xi + 0.5)*(0.952120850728289*xi + 0.271574874126072)
                                  + 0.952120850728289*(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi
                                  + 0.566554479308753)*(0.6535475074298*xi + 0.5)*(2.0840983387567*xi - 0.594450529658367)
                                  + 0.566554479308753*(4.25632117622354 - 4.25632117622354*xi)*(0.6535475074298*xi + 0.5)
                                  *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367)
                                  + (-2.78170809554157*xi - 2.12816058811177)*(0.566554479308753*xi + 0.566554479308753)
                                  *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367);

      derivativeShapeFunction[5]=(0.283277239654376*xi + 0.216722760345624)*(0.778069933064125*xi + 0.221930066935875)
                                  *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354)
                                  + (0.283277239654376*xi + 0.283277239654376)*(0.778069933064125*xi + 0.221930066935875)
                                  *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354)
                                  + 4.25632117622354*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                  *(0.778069933064125*xi + 0.221930066935875)*(1.39905441140358*xi - 0.399054411403579)
                                  + 1.39905441140358*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                  *(0.778069933064125*xi + 0.221930066935875)*(4.25632117622354*xi - 3.25632117622354)
                                  + 0.778069933064125*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                  *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354);
    }
    return derivativeShapeFunction  ;
  }

  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints,
                                     arrayDouble const & derivativeBasisFunction1D ) const
  {
    // loop over quadrature points
    for( int i = 0; i < order+1; i++ )
    {
      std::vector< double > tmp( order+1 );
      //extract all basis functions  for current quadrature point
      tmp=derivativeShapeFunction1D( order, quadraturePoints[i] );
      for( int j=0; j<order+1; j++ )
      {
        derivativeBasisFunction1D( j, i )=tmp[j];
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////
  //  end from first implementation
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  //  from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
  
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

  /////////////////////////////////////////////////////////////////////////////////////
  //  end from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
};
  
#endif //SEMQKGLBASISFUNCTIONS_HPP_
