#ifndef TENSOROPS_H_
#define TENSOROPS_H_

// #include <commonMacro.hpp>
/**
 * @brief Inverts a 3x3 matrix and returns its determinant.
 *        The matrix is modified in place.
 * @param J The 3x3 matrix to invert.
 * @return The determinant of the original matrix.
 */
PROXY_HOST_DEVICE
double invert3x3(double (&J)[3][3])
{
  // Compute the determinant
  double det =
      J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
      J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
      J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

  double invDet = 1.0 / det;

  double inv[3][3];

  inv[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * invDet;
  inv[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * invDet;
  inv[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * invDet;

  inv[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * invDet;
  inv[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * invDet;
  inv[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * invDet;

  inv[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * invDet;
  inv[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * invDet;
  inv[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * invDet;

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      J[i][j] = inv[i][j];

  return det;
}

/**
* @brief Invert the source matrix @p J and store the result in @p Jinv.
* @tparam Jinv The type of @p Jinv.
* @tparam J The type of @p J.
* @param Jinv The 3v3 matrix to write the inverse to.
* @param srcMatrix The 3x3 matrix to take the inverse of.
* @return The determinant.
* @note @p srcMatrix can contain integers but @p dstMatrix must contain floating point values.
*/
PROXY_HOST_DEVICE
auto invert3x3(double (&Jinv)[3][3], double (&J)[3][3])
{
    Jinv[ 0 ][ 0 ] = J[ 1 ][ 1 ] * J[ 2 ][ 2 ] - J[ 1 ][ 2 ] * J[ 2 ][ 1 ];
    Jinv[ 0 ][ 1 ] = J[ 0 ][ 2 ] * J[ 2 ][ 1 ] - J[ 0 ][ 1 ] * J[ 2 ][ 2 ];
    Jinv[ 0 ][ 2 ] = J[ 0 ][ 1 ] * J[ 1 ][ 2 ] - J[ 0 ][ 2 ] * J[ 1 ][ 1 ];

    auto const det = J[ 0 ][ 0 ] * Jinv[ 0 ][ 0 ] +
                     J[ 1 ][ 0 ] * Jinv[ 0 ][ 1 ] +
                     J[ 2 ][ 0 ] * Jinv[ 0 ][ 2 ];

    auto const invDet = 1.0 / det;

    Jinv[ 0 ][ 0 ] *= invDet;
    Jinv[ 0 ][ 1 ] *= invDet;
    Jinv[ 0 ][ 2 ] *= invDet;
    Jinv[ 1 ][ 0 ] = ( J[ 1 ][ 2 ] * J[ 2 ][ 0 ] - J[ 1 ][ 0 ] * J[ 2 ][ 2 ] ) * invDet;
    Jinv[ 1 ][ 1 ] = ( J[ 0 ][ 0 ] * J[ 2 ][ 2 ] - J[ 0 ][ 2 ] * J[ 2 ][ 0 ] ) * invDet;
    Jinv[ 1 ][ 2 ] = ( J[ 0 ][ 2 ] * J[ 1 ][ 0 ] - J[ 0 ][ 0 ] * J[ 1 ][ 2 ] ) * invDet;
    Jinv[ 2 ][ 0 ] = ( J[ 1 ][ 0 ] * J[ 2 ][ 1 ] - J[ 1 ][ 1 ] * J[ 2 ][ 0 ] ) * invDet;
    Jinv[ 2 ][ 1 ] = ( J[ 0 ][ 1 ] * J[ 2 ][ 0 ] - J[ 0 ][ 0 ] * J[ 2 ][ 1 ] ) * invDet;
    Jinv[ 2 ][ 2 ] = ( J[ 0 ][ 0 ] * J[ 1 ][ 1 ] - J[ 0 ][ 1 ] * J[ 1 ][ 0 ] ) * invDet;


    return det;
}

template< int N >
PROXY_HOST_DEVICE
double determinant(const double (&A)[N][N]);

template<>
PROXY_HOST_DEVICE
double determinant<3>(const double (&A)[3][3])
{
  return
      A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
      A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
      A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

template< int N >
PROXY_HOST_DEVICE
double symDeterminant( double (&B)[3]);

template<>
PROXY_HOST_DEVICE
double symDeterminant<2>( double (&B)[3] )
{
  return B[0] * B[1] - B[2] * B[2];
}

template< int N >
PROXY_HOST_DEVICE
double symDeterminant( double (&B)[6]);

template<>
PROXY_HOST_DEVICE
double symDeterminant<3>( double (&B)[6] )
{
  return B[ 0 ] * B[ 1 ] * B[ 2 ] +
         B[ 5 ] * B[ 4 ] * B[ 3 ] * 2 -
         B[ 0 ] * B[ 3 ] * B[ 3 ] -
         B[ 1 ] * B[ 4 ] * B[ 4 ] -
         B[ 2 ] * B[ 5 ] * B[ 5 ];
}

PROXY_HOST_DEVICE
static auto symInvert( double (&J)[6])
{
  auto const det = symDeterminant<3>(J);
  auto const invDet = 1 / det;

  auto const temp = J[ 0 ];
  J[ 0 ] = J[ 1 ] * invDet;
  J[ 1 ] = temp * invDet;
  J[ 2 ] *= -invDet;

  return det;
}

#endif // TENSOROPS_H_
