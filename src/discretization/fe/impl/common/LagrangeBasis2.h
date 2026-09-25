#ifndef LAGRANGEBASIS2_HPP_
#define LAGRANGEBASIS2_HPP_

/**
 * @brief Quadratic (order 2) Lagrange basis on the nodes -1, 0 and 1 of
 * [-1, 1].
 *
 * See docs/design.md, "1D Lagrange bases".
 */
class LagrangeBasis2 {
 public:
  constexpr static int numSupportPoints = 3;  ///< Number of nodes, order + 1.

  /**
   * @brief Gauss-Lobatto quadrature weight of node @p q on [-1, 1].
   */
  PROXY_HOST_DEVICE
  constexpr static real_t weight(const int q) {
    switch (q) {
      case 0:
      case 2:
        return 1.0 / 3.0;
      default:
        return 4.0 / 3.0;
    }
  }

  /**
   * @brief Parent coordinate of node @p supportPointIndex, in [-1, 1].
   */
  PROXY_HOST_DEVICE
  constexpr static double parentSupportCoord(const int supportPointIndex) {
    switch (supportPointIndex) {
      case 0:
        return -1.0;
        break;
      case 2:
        return 1.0;
      case 1:
      default:
        return 0.0;
    }
  }

  /**
   * @brief Value at @p xi of the basis function of node @p index.
   */
  PROXY_HOST_DEVICE
  constexpr static double value(const int index, const double xi) {
    switch (index) {
      case 0:
        return value0(xi);
      case 2:
        return value2(xi);
      case 1:
      default:
        return value1(xi);
    }
  }

  /**
   * @brief Value at @p xi of the basis function of node 0.
   */
  PROXY_HOST_DEVICE
  constexpr static double value0(const double xi) {
    const double xi_div2 = 0.5 * xi;
    return -xi_div2 + xi_div2 * xi;
  }

  /**
   * @brief Value at @p xi of the basis function of node 1.
   */
  PROXY_HOST_DEVICE
  constexpr static double value1(const double xi) { return 1.0 - xi * xi; }

  /**
   * @brief Value at @p xi of the basis function of node 2.
   */
  PROXY_HOST_DEVICE
  constexpr static double value2(const double xi) {
    const double xi_div2 = 0.5 * xi;
    return xi_div2 + xi_div2 * xi;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 0.
   */
  PROXY_HOST_DEVICE
  constexpr static double gradient0(const double xi) { return -0.5 + xi; }

  /**
   * @brief Derivative at @p xi of the basis function of node 1.
   */
  PROXY_HOST_DEVICE
  constexpr static double gradient1(const double xi) { return -2 * xi; }

  /**
   * @brief Derivative at @p xi of the basis function of node 2.
   */
  PROXY_HOST_DEVICE
  constexpr static double gradient2(const double xi) { return 0.5 + xi; }

  /**
   * @brief Derivative at @p xi of the basis function of node @p index.
   */
  PROXY_HOST_DEVICE
  constexpr static double gradient(const int index, const double xi) {
    switch (index) {
      case 0:
        return gradient0(xi);
      case 2:
        return gradient2(xi);
      case 1:
      default:
        return gradient1(xi);
    }
  }

  /**
   * @brief Derivative of the basis function of node @p q at node @p p.
   *
   * @pre p <= (numSupportPoints - 1) / 2: other values of @p p return a
   * meaningless result. See docs/design.md, "1D Lagrange bases".
   */
  PROXY_HOST_DEVICE
  constexpr static double gradientAt(const int q, const int p) {
    switch (q) {
      case 0:
        return p == 0 ? -1.5 : -0.5;
      case 1:
        return p == 0 ? 2.0 : 0.0;
      case 2:
        return p == 0 ? -0.5 : 0.5;
      default:
        return 0;
    }
  }

  /**
   * @brief Tensor product of the 1D basis on the parent square [-1, 1]^2.
   *
   * See docs/design.md, "1D Lagrange bases".
   */
  struct TensorProduct2D {
    constexpr static int numSupportPoints = 9;  ///< Number of nodes, 3^2.

    /**
     * @brief Index i + 3*j of the node (i, j), with i and j in [0, 2].
     */
    PROXY_HOST_DEVICE
    constexpr static int linearIndex(const int i, const int j) { return i + 3 * j; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     */
    PROXY_HOST_DEVICE
    constexpr static void multiIndex(const int linearIndex, int &i0, int &i1) {
      // (x * 22) >> 6 == x / 3 for x in [0, 26].
      i1 = ((linearIndex * 22) >> 6);

      i0 = linearIndex - i1 * 3;
    }

    /**
     * @brief Values at one point of all the 2D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1).
     * @param[out] N N[linearIndex(a, b)] = value(a, xi0) * value(b, xi1).
     */
    PROXY_HOST_DEVICE
    static void value(double const (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          const int lindex = LagrangeBasis2::TensorProduct2D::linearIndex(a, b);
          N[lindex] = LagrangeBasis2::value(a, coords[0]) * LagrangeBasis2::value(b, coords[1]);
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
    constexpr static int numSupportPoints = 27;  ///< Number of nodes, 3^3.

    /**
     * @brief Index i + 3*j + 9*k of the node (i, j, k), with i, j and k in
     * [0, 2].
     *
     * See docs/design.md, "Hexahedron local numbering".
     */
    constexpr static int linearIndex(const int i, const int j, const int k) { return i + 3 * j + 9 * k; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     * @param[out] i2 Index along xi2.
     */
    constexpr static void multiIndex(const int linearIndex, int &i0, int &i1, int &i2) {
      // Divisions by shifts: (x * 29) >> 8 == x / 9 and (x * 22) >> 6 == x / 3
      // for x in [0, 26].
      i2 = (linearIndex * 29) >> 8;

      i1 = ((linearIndex * 22) >> 6) - i2 * 3;

      i0 = linearIndex - i1 * 3 - i2 * 9;
    }

    /**
     * @brief Values at one point of all the 3D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1, xi2).
     * @param[out] N N[linearIndex(a, b, c)] = value(a, xi0) * value(b, xi1) *
     * value(c, xi2).
     */
    PROXY_HOST_DEVICE
    static void value(const double (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          for (int c = 0; c < 3; ++c) {
            const int lindex = LagrangeBasis2::TensorProduct3D::linearIndex(a, b, c);
            N[lindex] = LagrangeBasis2::value(a, coords[0]) * LagrangeBasis2::value(b, coords[1]) *
                        LagrangeBasis2::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif /* LAGRANGEBASIS2_HPP_ */
