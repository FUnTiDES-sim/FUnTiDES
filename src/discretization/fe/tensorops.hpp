#ifndef TENSOROPS_H_
#define TENSOROPS_H_

/**
 * @brief Inverts a 3x3 matrix and returns its determinant.
 *        The matrix is modified in place.
 * @param J The 3x3 matrix to invert.
 * @return The determinant of the original matrix.
 */
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

template< int N >
double determinant(const double (&A)[N][N]);

template<>
inline double determinant<3>(const double (&A)[3][3])
{
  return
      A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
      A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
      A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

template< int N >
double symDeterminant(const double (&B)[3]);

template<>
inline double symDeterminant<2>( const double (&B)[3] )
{
  return B[0] * B[1] - B[2] * B[2];
}

#endif // TENSOROPS_H_
