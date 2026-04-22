#ifndef SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS6GL_H_
#define SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS6GL_H_

/**
 * @file LagrangeBasis6GL.hpp
 *
 * 1D Lagrange basis on Gauss-Lobatto-Legendre nodes for degree 6 (7 nodes).
 * Follows the same implementation pattern as LagrangeBasis5GL.hpp.
 * Polynomial formulas and gradientAt table generated via GL.py.
 *
 * GLL nodes: { -1, -lambda1, -lambda2, 0, lambda2, lambda1, 1 }
 *   lambda1 = 0.8302238962785670
 *   lambda2 = 0.4688487934707142
 */
class LagrangeBasis6GL
{
 public:
  constexpr static int numSupportPoints = 7;

  /// |nodes[5]| = |nodes[1]|
  static constexpr double lambda1 = 0.8302238962785670;

  /// |nodes[4]| = |nodes[2]|
  static constexpr double lambda2 = 0.4688487934707142;

  /**
   * @brief GLL quadrature weight for node q.
   */
  constexpr static double weight(const int q)
  {
    switch (q)
    {
      case 1:
      case 5:
        return 0.2768260473615657;
      case 2:
      case 4:
        return 0.4317453812098626;
      case 3:
        return 0.4876190476190476;
      default:
        return 1.0 / 21.0;
    }
  }

  /**
   * @brief Parent coordinate of support point @p supportPointIndex.
   */
  constexpr static double parentSupportCoord(const int supportPointIndex)
  {
    switch (supportPointIndex)
    {
      case 0:
        return -1.0;
      case 1:
        return -lambda1;
      case 2:
        return -lambda2;
      case 3:
        return 0.0;
      case 4:
        return lambda2;
      case 5:
        return lambda1;
      case 6:
        return 1.0;
      default:
        return 0.0;
    }
  }

  /**
   * @brief Value of basis function @p index at @p xi.
   */
  constexpr static double value(const int index, const double xi)
  {
    switch (index)
    {
      case 0:
        return value0(xi);
      case 1:
        return value1(xi);
      case 2:
        return value2(xi);
      case 3:
        return value3(xi);
      case 4:
        return value4(xi);
      case 5:
        return value5(xi);
      case 6:
        return value6(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Basis function for support point 0.
   */
  constexpr static double value0(const double xi)
  {
    return 2.0625 * xi * xi * xi * xi * xi * xi -
           2.0625 * xi * xi * xi * xi * xi - 1.875 * xi * xi * xi * xi +
           1.875 * xi * xi * xi + 0.3125 * xi * xi - 0.3125 * xi;
  }

  /**
   * @brief Basis function for support point 1.
   */
  constexpr static double value1(const double xi)
  {
    return -4.97286970608696 * xi * xi * xi * xi * xi * xi +
           4.12859526307317 * xi * xi * xi * xi * xi +
           6.06600190251835 * xi * xi * xi * xi -
           5.03613973434199 * xi * xi * xi - 1.09313219643140 * xi * xi +
           0.90754447126882 * xi;
  }

  /**
   * @brief Basis function for support point 2.
   */
  constexpr static double value2(const double xi)
  {
    return 6.21036970608696 * xi * xi * xi * xi * xi * xi -
           2.91172434370594 * xi * xi * xi * xi * xi -
           10.4910019025184 * xi * xi * xi * xi +
           4.91869358429470 * xi * xi * xi + 4.28063219643140 * xi * xi -
           2.00696924058875 * xi;
  }

  /**
   * @brief Basis function for support point 3 (center node).
   */
  constexpr static double value3(const double xi)
  {
    return -6.6 * xi * xi * xi * xi * xi * xi + 12.6 * xi * xi * xi * xi -
           7.0 * xi * xi + 1.0;
  }

  /**
   * @brief Basis function for support point 4.
   */
  constexpr static double value4(const double xi)
  {
    return 6.21036970608696 * xi * xi * xi * xi * xi * xi +
           2.91172434370594 * xi * xi * xi * xi * xi -
           10.4910019025184 * xi * xi * xi * xi -
           4.91869358429470 * xi * xi * xi + 4.28063219643140 * xi * xi +
           2.00696924058875 * xi;
  }

  /**
   * @brief Basis function for support point 5.
   */
  constexpr static double value5(const double xi)
  {
    return -4.97286970608696 * xi * xi * xi * xi * xi * xi -
           4.12859526307317 * xi * xi * xi * xi * xi +
           6.06600190251835 * xi * xi * xi * xi +
           5.03613973434199 * xi * xi * xi - 1.09313219643140 * xi * xi -
           0.90754447126882 * xi;
  }

  /**
   * @brief Basis function for support point 6.
   */
  constexpr static double value6(const double xi)
  {
    return 2.0625 * xi * xi * xi * xi * xi * xi +
           2.0625 * xi * xi * xi * xi * xi - 1.875 * xi * xi * xi * xi -
           1.875 * xi * xi * xi + 0.3125 * xi * xi + 0.3125 * xi;
  }

  /**
   * @brief Gradient of basis function @p index at @p xi.
   */
  constexpr static double gradient(const int index, const double xi)
  {
    switch (index)
    {
      case 0:
        return gradient0(xi);
      case 1:
        return gradient1(xi);
      case 2:
        return gradient2(xi);
      case 3:
        return gradient3(xi);
      case 4:
        return gradient4(xi);
      case 5:
        return gradient5(xi);
      case 6:
        return gradient6(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Gradient of basis function for support point 0.
   */
  constexpr static double gradient0(const double xi)
  {
    return 12.375 * xi * xi * xi * xi * xi - 10.3125 * xi * xi * xi * xi -
           7.5 * xi * xi * xi + 5.625 * xi * xi + 0.625 * xi - 0.3125;
  }

  /**
   * @brief Gradient of basis function for support point 1.
   */
  constexpr static double gradient1(const double xi)
  {
    return -29.8372182365218 * xi * xi * xi * xi * xi +
           20.6429763153658 * xi * xi * xi * xi +
           24.2640076100734 * xi * xi * xi - 15.1084192030260 * xi * xi -
           2.18626439286279 * xi + 0.90754447126882;
  }

  /**
   * @brief Gradient of basis function for support point 2.
   */
  constexpr static double gradient2(const double xi)
  {
    return 37.2622182365217 * xi * xi * xi * xi * xi -
           14.5586217185297 * xi * xi * xi * xi -
           41.9640076100734 * xi * xi * xi + 14.7560807528841 * xi * xi +
           8.56126439286279 * xi - 2.00696924058875;
  }

  /**
   * @brief Gradient of basis function for support point 3 (center node).
   */
  constexpr static double gradient3(const double xi)
  {
    return -39.6 * xi * xi * xi * xi * xi + 50.4 * xi * xi * xi - 14.0 * xi;
  }

  /**
   * @brief Gradient of basis function for support point 4.
   */
  constexpr static double gradient4(const double xi)
  {
    return 37.2622182365217 * xi * xi * xi * xi * xi +
           14.5586217185297 * xi * xi * xi * xi -
           41.9640076100734 * xi * xi * xi - 14.7560807528841 * xi * xi +
           8.56126439286279 * xi + 2.00696924058875;
  }

  /**
   * @brief Gradient of basis function for support point 5.
   */
  constexpr static double gradient5(const double xi)
  {
    return -29.8372182365218 * xi * xi * xi * xi * xi -
           20.6429763153658 * xi * xi * xi * xi +
           24.2640076100734 * xi * xi * xi + 15.1084192030260 * xi * xi -
           2.18626439286279 * xi - 0.90754447126882;
  }

  /**
   * @brief Gradient of basis function for support point 6.
   */
  constexpr static double gradient6(const double xi)
  {
    return 12.375 * xi * xi * xi * xi * xi + 10.3125 * xi * xi * xi * xi -
           7.5 * xi * xi * xi - 5.625 * xi * xi + 0.625 * xi + 0.3125;
  }

  /**
   * @brief Gradient of basis function q evaluated at support point p.
   * Full 7x7 table computed via GL.py.
   * Anti-symmetry: gradientAt(q,p) == -gradientAt(6-q, 6-p).
   */
  constexpr static double gradientAt(const int q, const int p)
  {
    switch (q)
    {
      case 0:
        switch (p)
        {
          case 0:
            return -10.5000000000000000;
          case 1:
            return -2.4429260142442892;
          case 2:
            return 0.6252566655153418;
          case 3:
            return -0.3125000000000000;
          case 4:
            return 0.2260994009425744;
          case 5:
            return -0.2266118703954456;
          case 6:
            return 0.5000000000000000;
        }
        break;
      case 1:
        switch (p)
        {
          case 0:
            return 14.2015766029198272;
          case 1:
            return 0.0;
          case 2:
            return -2.2158042831699714;
          case 3:
            return 0.9075444712688211;
          case 4:
            return -0.6163908355175793;
          case 5:
            return 0.6022471796357809;
          case 6:
            return -1.3173734357024478;
        }
        break;
      case 2:
        switch (p)
        {
          case 0:
            return -5.6689852255454998;
          case 1:
            return 3.4558282142942875;
          case 2:
            return 0.0;
          case 3:
            return -2.0069692405887531;
          case 4:
            return 1.0664419040063753;
          case 5:
            return -0.9613397972887139;
          case 6:
            return 2.0499648130767394;
        }
        break;
      case 3:
        switch (p)
        {
          case 0:
            return 3.2000000000000000;
          case 1:
            return -1.5986066880983891;
          case 2:
            return 2.2666980870859961;
          case 3:
            return 0.0;
          case 4:
            return -2.2666980870859943;
          case 5:
            return 1.5986066880983962;
          case 6:
            return -3.2000000000000000;
        }
        break;
      case 4:
        switch (p)
        {
          case 0:
            return -2.0499648130767412;
          case 1:
            return 0.9613397972887121;
          case 2:
            return -1.0664419040063757;
          case 3:
            return 2.0069692405887531;
          case 4:
            return 0.0;
          case 5:
            return -3.4558282142942858;
          case 6:
            return 5.6689852255455051;
        }
        break;
      case 5:
        switch (p)
        {
          case 0:
            return 1.3173734357024465;
          case 1:
            return -0.6022471796357802;
          case 2:
            return 0.6163908355175800;
          case 3:
            return -0.9075444712688211;
          case 4:
            return 2.2158042831699714;
          case 5:
            return 0.0;
          case 6:
            return -14.2015766029198272;
        }
        break;
      case 6:
        switch (p)
        {
          case 0:
            return -0.5000000000000000;
          case 1:
            return 0.2266118703954456;
          case 2:
            return -0.2260994009425744;
          case 3:
            return 0.3125000000000000;
          case 4:
            return -0.6252566655153419;
          case 5:
            return 2.4429260142442892;
          case 6:
            return 10.5000000000000000;
        }
        break;
    }
    return 0.0;
  }

  /* Tensor product helpers (2D / 3D) */
  struct TensorProduct2D
  {
    constexpr static int numSupportPoints1D =
        LagrangeBasis6GL::numSupportPoints;
    constexpr static int numSupportPoints =
        numSupportPoints1D * numSupportPoints1D;  // 7*7 = 49

    constexpr static int linearIndex(const int i, const int j)
    {
      return i + numSupportPoints1D * j;
    }

    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1)
    {
      i1 = linearIndex / numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    static void value(const double (&coords)[2], double (&N)[numSupportPoints])
    {
      for (int a = 0; a < LagrangeBasis6GL::numSupportPoints; ++a)
      {
        for (int b = 0; b < LagrangeBasis6GL::numSupportPoints; ++b)
        {
          const int lindex = a + LagrangeBasis6GL::numSupportPoints * b;
          N[lindex] = LagrangeBasis6GL::value(a, coords[0]) *
                      LagrangeBasis6GL::value(b, coords[1]);
        }
      }
    }
  };

  struct TensorProduct3D
  {
    constexpr static int numSupportPoints1D =
        LagrangeBasis6GL::numSupportPoints;
    constexpr static int numSupportPoints = numSupportPoints1D *
                                            numSupportPoints1D *
                                            numSupportPoints1D;  // 7^3 = 343

    constexpr static int linearIndex(const int i, const int j, const int k)
    {
      return i + numSupportPoints1D * j +
             numSupportPoints1D * numSupportPoints1D * k;
    }

    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1,
                                     int& i2)
    {
      i2 = linearIndex / (numSupportPoints1D * numSupportPoints1D);
      i1 = (linearIndex % (numSupportPoints1D * numSupportPoints1D)) /
           numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    static void value(const double (&coords)[3], double (&N)[numSupportPoints])
    {
      for (int a = 0; a < LagrangeBasis6GL::numSupportPoints; ++a)
      {
        for (int b = 0; b < LagrangeBasis6GL::numSupportPoints; ++b)
        {
          for (int c = 0; c < LagrangeBasis6GL::numSupportPoints; ++c)
          {
            const int lindex = a + LagrangeBasis6GL::numSupportPoints * b +
                               LagrangeBasis6GL::numSupportPoints *
                                   LagrangeBasis6GL::numSupportPoints * c;
            N[lindex] = LagrangeBasis6GL::value(a, coords[0]) *
                        LagrangeBasis6GL::value(b, coords[1]) *
                        LagrangeBasis6GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif  // SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS6GL_H_
