#ifndef SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS6GL_H_
#define SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS6GL_H_

/**
 * @brief Order 6 Lagrange basis on the 7 Gauss-Lobatto-Legendre nodes -1,
 * -lambda1, -lambda2, 0, lambda2, lambda1 and 1.
 *
 * See docs/design.md, "1D Lagrange bases".
 */
class LagrangeBasis6GL {
 public:
  constexpr static int numSupportPoints = 7;  ///< Number of nodes, order + 1.

  static constexpr double lambda1 = 0.8302238962785670;  ///< Parent coordinate of node 5, opposite of node 1.

  static constexpr double lambda2 = 0.4688487934707142;  ///< Parent coordinate of node 4, opposite of node 2.

  /**
   * @brief Gauss-Lobatto quadrature weight of node @p q on [-1, 1].
   */
  constexpr static double weight(const int q) {
    switch (q) {
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
   * @brief Parent coordinate of node @p supportPointIndex, in [-1, 1].
   */
  constexpr static double parentSupportCoord(const int supportPointIndex) {
    switch (supportPointIndex) {
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
   * @brief Value at @p xi of the basis function of node @p index.
   */
  constexpr static double value(const int index, const double xi) {
    switch (index) {
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
   * @brief Value at @p xi of the basis function of node 0.
   */
  constexpr static double value0(const double xi) {
    return 2.0625 * xi * xi * xi * xi * xi * xi - 2.0625 * xi * xi * xi * xi * xi - 1.875 * xi * xi * xi * xi +
           1.875 * xi * xi * xi + 0.3125 * xi * xi - 0.3125 * xi;
  }

  /**
   * @brief Value at @p xi of the basis function of node 1.
   */
  constexpr static double value1(const double xi) {
    return -4.97286970608696 * xi * xi * xi * xi * xi * xi + 4.12859526307317 * xi * xi * xi * xi * xi +
           6.06600190251835 * xi * xi * xi * xi - 5.03613973434199 * xi * xi * xi - 1.09313219643140 * xi * xi +
           0.90754447126882 * xi;
  }

  /**
   * @brief Value at @p xi of the basis function of node 2.
   */
  constexpr static double value2(const double xi) {
    return 6.21036970608696 * xi * xi * xi * xi * xi * xi - 2.91172434370594 * xi * xi * xi * xi * xi -
           10.4910019025184 * xi * xi * xi * xi + 4.91869358429470 * xi * xi * xi + 4.28063219643140 * xi * xi -
           2.00696924058875 * xi;
  }

  /**
   * @brief Value at @p xi of the basis function of node 3.
   */
  constexpr static double value3(const double xi) {
    return -6.6 * xi * xi * xi * xi * xi * xi + 12.6 * xi * xi * xi * xi - 7.0 * xi * xi + 1.0;
  }

  /**
   * @brief Value at @p xi of the basis function of node 4.
   */
  constexpr static double value4(const double xi) {
    return 6.21036970608696 * xi * xi * xi * xi * xi * xi + 2.91172434370594 * xi * xi * xi * xi * xi -
           10.4910019025184 * xi * xi * xi * xi - 4.91869358429470 * xi * xi * xi + 4.28063219643140 * xi * xi +
           2.00696924058875 * xi;
  }

  /**
   * @brief Value at @p xi of the basis function of node 5.
   */
  constexpr static double value5(const double xi) {
    return -4.97286970608696 * xi * xi * xi * xi * xi * xi - 4.12859526307317 * xi * xi * xi * xi * xi +
           6.06600190251835 * xi * xi * xi * xi + 5.03613973434199 * xi * xi * xi - 1.09313219643140 * xi * xi -
           0.90754447126882 * xi;
  }

  /**
   * @brief Value at @p xi of the basis function of node 6.
   */
  constexpr static double value6(const double xi) {
    return 2.0625 * xi * xi * xi * xi * xi * xi + 2.0625 * xi * xi * xi * xi * xi - 1.875 * xi * xi * xi * xi -
           1.875 * xi * xi * xi + 0.3125 * xi * xi + 0.3125 * xi;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node @p index.
   */
  constexpr static double gradient(const int index, const double xi) {
    switch (index) {
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
   * @brief Derivative at @p xi of the basis function of node 0.
   */
  constexpr static double gradient0(const double xi) {
    return 12.375 * xi * xi * xi * xi * xi - 10.3125 * xi * xi * xi * xi - 7.5 * xi * xi * xi + 5.625 * xi * xi +
           0.625 * xi - 0.3125;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 1.
   */
  constexpr static double gradient1(const double xi) {
    return -29.8372182365218 * xi * xi * xi * xi * xi + 20.6429763153658 * xi * xi * xi * xi +
           24.2640076100734 * xi * xi * xi - 15.1084192030260 * xi * xi - 2.18626439286279 * xi + 0.90754447126882;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 2.
   */
  constexpr static double gradient2(const double xi) {
    return 37.2622182365217 * xi * xi * xi * xi * xi - 14.5586217185297 * xi * xi * xi * xi -
           41.9640076100734 * xi * xi * xi + 14.7560807528841 * xi * xi + 8.56126439286279 * xi - 2.00696924058875;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 3.
   */
  constexpr static double gradient3(const double xi) {
    return -39.6 * xi * xi * xi * xi * xi + 50.4 * xi * xi * xi - 14.0 * xi;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 4.
   */
  constexpr static double gradient4(const double xi) {
    return 37.2622182365217 * xi * xi * xi * xi * xi + 14.5586217185297 * xi * xi * xi * xi -
           41.9640076100734 * xi * xi * xi - 14.7560807528841 * xi * xi + 8.56126439286279 * xi + 2.00696924058875;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 5.
   */
  constexpr static double gradient5(const double xi) {
    return -29.8372182365218 * xi * xi * xi * xi * xi - 20.6429763153658 * xi * xi * xi * xi +
           24.2640076100734 * xi * xi * xi + 15.1084192030260 * xi * xi - 2.18626439286279 * xi - 0.90754447126882;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 6.
   */
  constexpr static double gradient6(const double xi) {
    return 12.375 * xi * xi * xi * xi * xi + 10.3125 * xi * xi * xi * xi - 7.5 * xi * xi * xi - 5.625 * xi * xi +
           0.625 * xi + 0.3125;
  }

  /**
   * @brief Derivative of the basis function of node @p q at node @p p.
   *
   * Tabulated for every node @p p. See docs/design.md, "1D Lagrange bases".
   */
  constexpr static double gradientAt(const int q, const int p) {
    switch (q) {
      case 0:
        switch (p) {
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
        switch (p) {
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
        switch (p) {
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
        switch (p) {
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
        switch (p) {
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
        switch (p) {
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
        switch (p) {
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

  /**
   * @brief Tensor product of the 1D basis on the parent square [-1, 1]^2.
   *
   * See docs/design.md, "1D Lagrange bases".
   */
  struct TensorProduct2D {
    constexpr static int numSupportPoints1D = LagrangeBasis6GL::numSupportPoints;  ///< Number of 1D nodes.
    constexpr static int numSupportPoints = numSupportPoints1D * numSupportPoints1D;  ///< Number of nodes, 7^2.

    /**
     * @brief Index i + 7*j of the node (i, j), with i and j in [0, 6].
     */
    constexpr static int linearIndex(const int i, const int j) { return i + numSupportPoints1D * j; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     */
    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1) {
      i1 = linearIndex / numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    /**
     * @brief Values at one point of all the 2D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1).
     * @param[out] N N[linearIndex(a, b)] = value(a, xi0) * value(b, xi1).
     */
    static void value(const double (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < LagrangeBasis6GL::numSupportPoints; ++a) {
        for (int b = 0; b < LagrangeBasis6GL::numSupportPoints; ++b) {
          const int lindex = a + LagrangeBasis6GL::numSupportPoints * b;
          N[lindex] = LagrangeBasis6GL::value(a, coords[0]) * LagrangeBasis6GL::value(b, coords[1]);
        }
      }
    }
  };

  /**
   * @brief Tensor product of the 1D basis on the parent cube [-1, 1]^3.
   *
   * See docs/design.md, "1D Lagrange bases".
   */
  struct TensorProduct3D {
    constexpr static int numSupportPoints1D = LagrangeBasis6GL::numSupportPoints;  ///< Number of 1D nodes.
    /// Number of nodes, 7^3.
    constexpr static int numSupportPoints = numSupportPoints1D * numSupportPoints1D * numSupportPoints1D;

    /**
     * @brief Index i + 7*j + 49*k of the node (i, j, k), with i, j and k in
     * [0, 6].
     *
     * See docs/design.md, "Hexahedron local numbering".
     */
    constexpr static int linearIndex(const int i, const int j, const int k) {
      return i + numSupportPoints1D * j + numSupportPoints1D * numSupportPoints1D * k;
    }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     * @param[out] i2 Index along xi2.
     */
    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1, int& i2) {
      i2 = linearIndex / (numSupportPoints1D * numSupportPoints1D);
      i1 = (linearIndex % (numSupportPoints1D * numSupportPoints1D)) / numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    /**
     * @brief Values at one point of all the 3D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1, xi2).
     * @param[out] N N[linearIndex(a, b, c)] = value(a, xi0) * value(b, xi1) *
     * value(c, xi2).
     */
    static void value(const double (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < LagrangeBasis6GL::numSupportPoints; ++a) {
        for (int b = 0; b < LagrangeBasis6GL::numSupportPoints; ++b) {
          for (int c = 0; c < LagrangeBasis6GL::numSupportPoints; ++c) {
            const int lindex = a + LagrangeBasis6GL::numSupportPoints * b +
                               LagrangeBasis6GL::numSupportPoints * LagrangeBasis6GL::numSupportPoints * c;
            N[lindex] = LagrangeBasis6GL::value(a, coords[0]) * LagrangeBasis6GL::value(b, coords[1]) *
                        LagrangeBasis6GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif  // SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS6GL_H_
