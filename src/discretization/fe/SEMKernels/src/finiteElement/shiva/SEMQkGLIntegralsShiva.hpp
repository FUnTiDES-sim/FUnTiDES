#pragma once

#include "common/macros.hpp"
#include "common/mathUtilites.hpp"

#include "shiva/functions/bases/LagrangeBasis.hpp"
#include "shiva/functions/quadrature/Quadrature.hpp"
#include "shiva/functions/spacing/Spacing.hpp"
#include "shiva/geometry/shapes/NCube.hpp"
#include "shiva/geometry/shapes/InterpolatedShape.hpp"
#include "shiva/geometry/mapping/LinearTransform.hpp"
#include "shiva/common/ShivaMacros.hpp"
#include "shiva/common/pmpl.hpp"
#include "shiva/common/types.hpp"
#include "shiva/discretizations/finiteElementMethod/parentElements/ParentElement.hpp"

#include <stdio.h>

using namespace shiva;
using namespace shiva::functions;
using namespace shiva::geometry;
using namespace shiva::geometry::utilities;
using namespace shiva::discretizations::finiteElementMethod;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
template< int ORDER,
          typename TRANSFORM,
          typename PARENT_ELEMENT >
class SEMQkGLIntegralsShiva
{
public:
  constexpr static bool isShiva = true;

  static constexpr int order = ORDER;
  static constexpr int numSupportPoints1d = ORDER + 1;

  using TransformType = TRANSFORM;
  using ParentElementType = PARENT_ELEMENT;

  using tfloat = typename TransformType::RealType;
  using gfloat = typename ParentElementType::RealType;

  using JacobianType = typename std::remove_reference_t< TransformType >::JacobianType;
  using quadrature = QuadratureGaussLobatto< gfloat, numSupportPoints1d >;
  using basisFunction = LagrangeBasis< gfloat, ORDER, GaussLobattoSpacing >;


  template< int qa, int qb, int qc, typename FUNC1, typename FUNC2 >
  static constexpr inline
  SEMKERNELS_HOST_DEVICE
  void computeGradPhiBGradPhi( tfloat const (&B)[6],
                               FUNC1 && func1,
                               FUNC2 && func2 )
  {
    constexpr gfloat qcoords0 = quadrature::template coordinate< qa >();
    constexpr gfloat qcoords1 = quadrature::template coordinate< qb >();
    constexpr gfloat qcoords2 = quadrature::template coordinate< qc >();
    constexpr gfloat w = quadrature::template weight< qa >() * quadrature::template weight< qb >() * quadrature::template weight< qc >();
    func1( qa, qb, qc);
    forSequence< numSupportPoints1d >( [&] ( auto const ici )
    {
      constexpr int i = decltype(ici)::value;      
      const int ibc = linearIndex<ORDER>( i, qb, qc );
      const int aic = linearIndex<ORDER>( qa, i, qc );
      const int abi = linearIndex<ORDER>( qa, qb, i );
      const gfloat gia = basisFunction::template gradient< i >( qcoords0 );
      const gfloat gib = basisFunction::template gradient< i >( qcoords1 );
      const gfloat gic = basisFunction::template gradient< i >( qcoords2 );
//      printf("i: %d, ibc: %d, aic: %d, abi: %d, gia: %f, gib: %f, gic: %f\n", i, ibc, aic, abi, gia, gib, gic);

      forSequence< numSupportPoints1d >( [&] ( auto const icj )
      {
        constexpr int j = decltype(icj)::value;
        const int jbc = linearIndex<ORDER>( j, qb, qc );
        const int ajc = linearIndex<ORDER>( qa, j, qc );
        const int abj = linearIndex<ORDER>( qa, qb, j );
        const gfloat gja = basisFunction::template gradient< j >( qcoords0 );
        const gfloat gjb = basisFunction::template gradient< j >( qcoords1 );
        const gfloat gjc = basisFunction::template gradient< j >( qcoords2 );

//        printf("j: %d, jbc: %d, ajc: %d, abj: %d, gja: %f, gjb: %f, gjc: %f\n", j, jbc, ajc, abj, gja, gjb, gjc);
        // diagonal terms
        const gfloat w0 = w * gia * gja;
        func2( ibc, jbc, w0 * B[0] );
        const gfloat w1 = w * gib * gjb;
        func2( aic, ajc, w1 * B[1] );
        const gfloat w2 = w * gic * gjc;
        func2( abi, abj, w2 * B[2] );
        // off-diagonal terms
        const gfloat w3 = w * gib * gjc;
        func2( aic, abj, w3 * B[3] );
        func2( abj, aic, w3 * B[3] );
        const gfloat w4 = w * gia * gjc;
        func2( ibc, abj, w4 * B[4] );
        func2( abj, ibc, w4 * B[4] );
        const gfloat w5 = w * gia * gjb;
        func2( ibc, ajc, w5 * B[5] );
        func2( ajc, ibc, w5 * B[5] );
      } );
    } );
  }


  template< int qa, int qb, int qc, typename FUNC1, typename FUNC2 >
  static constexpr inline
  SEMKERNELS_HOST_DEVICE
  void
  computeGradPhiGradPhi( JacobianType & J,
                        FUNC1 && func1,
                        FUNC2 && func2 )
  {
    real_t detJ;
    shiva::mathUtilities::inverse( J, detJ );
    func1( qa, qb, qc, *reinterpret_cast< real_t const (*)[3][3] >(J.data()) );
    constexpr gfloat qcoords0 = quadrature::template coordinate< qa >();
    constexpr gfloat qcoords1 = quadrature::template coordinate< qb >();
    constexpr gfloat qcoords2 = quadrature::template coordinate< qc >();
    constexpr gfloat w = quadrature::template weight< qa >() * quadrature::template weight< qb >() * quadrature::template weight< qc >();
    forSequence< numSupportPoints1d >( [&] ( auto const ici )
    {
      constexpr int i = decltype(ici)::value;      
      const int ibc = linearIndex<ORDER>( i, qb, qc );
      const int aic = linearIndex<ORDER>( qa, i, qc );
      const int abi = linearIndex<ORDER>( qa, qb, i );
      const gfloat gia = basisFunction::template gradient< i >( qcoords0 );
      const gfloat gib = basisFunction::template gradient< i >( qcoords1 );
      const gfloat gic = basisFunction::template gradient< i >( qcoords2 );

      forSequence< numSupportPoints1d >( [&] ( auto const icj )
      {
        constexpr int j = decltype(icj)::value;
        const int jbc = linearIndex<ORDER>( j, qb, qc );
        const int ajc = linearIndex<ORDER>( qa, j, qc );
        const int abj = linearIndex<ORDER>( qa, qb, j );
        const gfloat gja = basisFunction::template gradient< j >( qcoords0 );
        const gfloat gjb = basisFunction::template gradient< j >( qcoords1 );
        const gfloat gjc = basisFunction::template gradient< j >( qcoords2 );

        // diagonal terms
        const real_t w00 = w * gia * gja;
        func2(ibc, jbc, w00 * detJ, 0, 0 );
        const real_t w11 = w * gib * gjb;
        func2(aic, ajc, w11 * detJ, 1, 1 );
        const real_t w22 = w * gic * gjc;
        func2(abi, abj, w22 * detJ, 2, 2 );
        // off-diagonal terms
        const real_t w12 = w * gib * gjc;
        func2(aic, abj, w12 * detJ, 1, 2 );
        func2(abj, aic, w12 * detJ, 2, 1 );
        const real_t w02 = w * gia * gjc;
        func2(ibc, abj, w02 * detJ, 0, 2 );
        func2(abj, ibc, w02 * detJ, 2, 0 );
        const real_t w01 = w * gia * gjb;
        func2(ibc, ajc, w01 * detJ, 0, 1 );
        func2(ajc, ibc, w01 * detJ, 1, 0 );
      });
    });
  }


  template< typename FUNC >
  static constexpr inline
  SEMKERNELS_HOST_DEVICE
  void computeMassTerm( TransformType const & trilinearCell,
                        FUNC && func )
  {

    // this is a compile time quadrature loop over each tensor direction
    forNestedSequence< ORDER + 1,
                       ORDER + 1,
                       ORDER + 1 >( [&] ( auto const icqa,
                                          auto const icqb,
                                          auto const icqc )
    {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;

      JacobianType J{ tfloat(0.0) };

      shiva::geometry::utilities::jacobian< quadrature, qa, qb, qc >( trilinearCell, J );

      tfloat const detJ = determinant( J );
      
      // mass matrix
      constexpr int q = linearIndex<ORDER>( qc, qb, qa );
      constexpr tfloat w3D = quadrature::template weight< qa >() *
                             quadrature::template weight< qb >() *
                             quadrature::template weight< qc >();
      func( q, w3D * detJ );
    } );
  }



  template< typename FUNC1, typename FUNC2 >
  static constexpr inline
  SEMKERNELS_HOST_DEVICE
  void computeStiffnessTerm( TransformType const & trilinearCell,
                             FUNC1 && func1, FUNC2 && func2)
  {

    // this is a compile time quadrature loop over each tensor direction
    forNestedSequence< ORDER + 1,
                       ORDER + 1,
                       ORDER + 1 >( [&] ( auto const icqa,
                                          auto const icqb,
                                          auto const icqc )
    {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;

      JacobianType J{ tfloat(0.0) };

      shiva::geometry::utilities::jacobian< quadrature, qa, qb, qc >( trilinearCell, J );

      tfloat const detJ = determinant( J );
      
      tfloat B[6] = {0};
      computeB( J, B );

      // compute detJ*J^{-1}J^{-T}
      for( int i = 0; i < 6; ++i )
      {
       B[i] *= detJ;
      }

      // compute gradPhiI*B*gradPhiJ and stiffness vector
      computeGradPhiBGradPhi< qa, qb, qc >( B, func1, func2 );
    } );
  }


  template< typename FUNC1, typename FUNC2 >
  static constexpr inline
  SEMKERNELS_HOST_DEVICE
  void computeStiffNessTermwithJac( TransformType const & trilinearCell,
                                    FUNC1 && func1, 
                                    FUNC2 && func2 )
  {

    // this is a compile time quadrature loop over each tensor direction
    forNestedSequence< ORDER + 1,
                       ORDER + 1,
                       ORDER + 1 >( [&] ( auto const icqa,
                                          auto const icqb,
                                          auto const icqc )
    {
      constexpr int qa = decltype(icqa)::value;
      constexpr int qb = decltype(icqb)::value;
      constexpr int qc = decltype(icqc)::value;

      JacobianType J{ tfloat(0.0) };
      shiva::geometry::utilities::jacobian< quadrature, qa, qb, qc >( trilinearCell, J );

      // compute gradPhiI*B*gradPhiJ and stiffness vector
      computeGradPhiGradPhi<qa,qb,qc>( J, func1, func2 );
    } );
  }


};
