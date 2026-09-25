#ifndef _LAGRANGEBASIS3GL_HPP_
#define _LAGRANGEBASIS3GL_HPP_

/**
 * @brief Cubic (order 3) Lagrange basis on the 4 Gauss-Lobatto-Legendre nodes
 * -1, -1/sqrt(5), 1/sqrt(5) and 1.
 *
 * See docs/design.md, "1D Lagrange bases".
 */
class LagrangeBasis3GL {
 public:
  constexpr static int numSupportPoints = 4;  ///< Number of nodes, order + 1.

  constexpr static double sqrt5 = 2.2360679774997897;  ///< sqrt(5).

  /**
   * @brief Gauss-Lobatto quadrature weight of node @p q on [-1, 1].
   */
  inline constexpr static double weight(const int q) {
    switch (q) {
      case 1:
      case 2:
        return 5.0 / 6.0;
      default:
        return 1.0 / 6.0;
    }
  }

  /**
   * @brief Parent coordinate of node @p supportPointIndex, in [-1, 1].
   */
  inline constexpr static double parentSupportCoord(const int supportPointIndex) {
    double result = 0.0;

    switch (supportPointIndex) {
      case 0:
        result = -1.0;
        break;
      case 1:
        result = -1.0 / sqrt5;
        break;
      case 2:
        result = 1.0 / sqrt5;
        break;
      case 3:
        result = 1.0;
        break;
      default:
        break;
    }

    return result;
  }

  /**
   * @brief Value at @p xi of the basis function of node @p index.
   */
  inline constexpr static double value(const int index, const double xi) {
    double result = 0.0;

    switch (index) {
      case 0:
        result = LagrangeBasis3GL::value0(xi);
        break;
      case 1:
        result = LagrangeBasis3GL::value1(xi);
        break;
      case 2:
        result = LagrangeBasis3GL::value2(xi);
        break;
      case 3:
        result = LagrangeBasis3GL::value3(xi);
        break;
      default:
        break;
    }

    return result;
  }

  /**
   * @brief Value at @p xi of the basis function of node 0.
   */
  inline constexpr static double value0(const double xi) {
    return -(5.0 / 8.0) * (xi * xi * xi - xi * xi - (1.0 / 5.0) * xi + 1.0 / 5.0);
  }

  /**
   * @brief Value at @p xi of the basis function of node 1.
   */
  inline constexpr static double value1(const double xi) {
    return (5.0 * sqrt5 / 8.0) * (xi * xi * xi - (1.0 / sqrt5) * xi * xi - xi + 1.0 / sqrt5);
  }

  /**
   * @brief Value at @p xi of the basis function of node 2.
   */
  inline constexpr static double value2(const double xi) {
    return -(5.0 * sqrt5 / 8.0) * (xi * xi * xi + (1.0 / sqrt5) * xi * xi - xi - 1.0 / sqrt5);
  }

  /**
   * @brief Value at @p xi of the basis function of node 3.
   */
  inline constexpr static double value3(const double xi) {
    return (5.0 / 8.0) * (xi * xi * xi + xi * xi - (1.0 / 5.0) * xi - 1.0 / 5.0);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node @p index.
   */
  inline constexpr static double gradient(const int index, const double xi) {
    double result = 0.0;

    switch (index) {
      case 0:
        result = LagrangeBasis3GL::gradient0(xi);
        break;
      case 1:
        result = LagrangeBasis3GL::gradient1(xi);
        break;
      case 2:
        result = LagrangeBasis3GL::gradient2(xi);
        break;
      case 3:
        result = LagrangeBasis3GL::gradient3(xi);
        break;
      default:
        break;
    }

    return result;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 0.
   */
  inline constexpr static double gradient0(const double xi) {
    return -(5.0 / 8.0) * (3.0 * xi * xi - 2.0 * xi - (1.0 / 5.0));
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 1.
   */
  inline constexpr static double gradient1(const double xi) {
    return (5.0 * sqrt5 / 8.0) * (3.0 * xi * xi - (2.0 / sqrt5) * xi - 1.0);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 2.
   */
  inline constexpr static double gradient2(const double xi) {
    return -(5.0 * sqrt5 / 8.0) * (3.0 * xi * xi + (2.0 / sqrt5) * xi - 1.0);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 3.
   */
  inline constexpr static double gradient3(const double xi) {
    return (5.0 / 8.0) * (3.0 * xi * xi + 2.0 * xi - (1.0 / 5.0));
    ;
  }

  /**
   * @brief Derivative of the basis function of node @p q at node @p p.
   *
   * @pre p <= (numSupportPoints - 1) / 2: other values of @p p return a
   * meaningless result. See docs/design.md, "1D Lagrange bases".
   */
  constexpr static double gradientAt(const int q, const int p) {
    switch (q) {
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
  }

  /**
   * @brief Tensor product of the 1D basis on the parent square [-1, 1]^2.
   *
   * See docs/design.md, "1D Lagrange bases".
   */
  struct TensorProduct2D {
    constexpr static int numSupportPoints = 16;  ///< Number of nodes, 4^2.

    /**
     * @brief Index i + 4*j of the node (i, j), with i and j in [0, 3].
     */
    inline constexpr static int linearIndex(const int i, const int j) { return i + 4 * j; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     */
    inline constexpr static void multiIndex(int const linearIndex, int &i0, int &i1) {
      i1 = linearIndex / 4;

      i0 = linearIndex % 4;
    }

    /**
     * @brief Values at one point of all the 2D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1).
     * @param[out] N N[linearIndex(a, b)] = value(a, xi0) * value(b, xi1).
     */
    inline static void value(const double (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
          const int lindex = LagrangeBasis3GL::TensorProduct2D::linearIndex(a, b);
          N[lindex] = LagrangeBasis3GL::value(a, coords[0]) * LagrangeBasis3GL::value(b, coords[1]);
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
    constexpr static int numSupportPoints = 64;  ///< Number of nodes, 4^3.

    /**
     * @brief Index i + 4*j + 16*k of the node (i, j, k), with i, j and k in
     * [0, 3].
     *
     * See docs/design.md, "Hexahedron local numbering".
     */
    inline constexpr static int linearIndex(const int i, const int j, const int k) { return i + 4 * j + 16 * k; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     * @param[out] i2 Index along xi2.
     */
    inline constexpr static void multiIndex(int const linearIndex, int &i0, int &i1, int &i2) {
      i2 = linearIndex / 16;

      i1 = (linearIndex % 16) / 4;

      i0 = (linearIndex % 16) % 4;
    }

    /**
     * @brief Values at one point of all the 3D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1, xi2).
     * @param[out] N N[linearIndex(a, b, c)] = value(a, xi0) * value(b, xi1) *
     * value(c, xi2).
     */
    PROXY_HOST_DEVICE
    static void value(const double (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
          for (int c = 0; c < 4; ++c) {
            const int lindex = LagrangeBasis3GL::TensorProduct3D::linearIndex(a, b, c);
            N[lindex] = LagrangeBasis3GL::value(a, coords[0]) * LagrangeBasis3GL::value(b, coords[1]) *
                        LagrangeBasis3GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif /* _LAGRANGEBASIS3GL_HPP_  */
