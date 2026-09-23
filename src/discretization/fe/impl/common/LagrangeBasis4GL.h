#ifndef _LAGRANGEBASIS4GL_HPP_
#define _LAGRANGEBASIS4GL_HPP_

/**
 * @brief Order 4 Lagrange basis on the 5 Gauss-Lobatto-Legendre nodes -1,
 * -sqrt(3/7), 0, sqrt(3/7) and 1.
 *
 * See docs/design.md, "1D Lagrange bases".
 */
class LagrangeBasis4GL {
 public:
  constexpr static int numSupportPoints = 5;  ///< Number of nodes, order + 1.

  constexpr static double sqrt3_7 = 0.6546536707079771;  ///< sqrt(3/7).

  /**
   * @brief Gauss-Lobatto quadrature weight of node @p q on [-1, 1].
   */
  constexpr static double weight(const int q) {
    switch (q) {
      case 0:
      case 4:
        return 1.0 / 10.0;
      case 1:
      case 3:
        return 49.0 / 90.0;
      default:
        return 32.0 / 45.0;
    }
  }

  /**
   * @brief Parent coordinate of node @p supportPointIndex, in [-1, 1].
   */
  constexpr static double parentSupportCoord(const int supportPointIndex) {
    double result = 0.0;

    switch (supportPointIndex) {
      case 0:
        result = -1.0;
        break;
      case 1:
        result = -sqrt3_7;
        break;
      case 2:
        result = 0.0;
        break;
      case 3:
        result = sqrt3_7;
        break;
      case 4:
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
  constexpr static double value(const int index, const double xi) {
    double result = 0.0;

    switch (index) {
      case 0:
        result = LagrangeBasis4GL::value0(xi);
        break;
      case 1:
        result = LagrangeBasis4GL::value1(xi);
        break;
      case 2:
        result = LagrangeBasis4GL::value2(xi);
        break;
      case 3:
        result = LagrangeBasis4GL::value3(xi);
        break;
      case 4:
        result = LagrangeBasis4GL::value4(xi);
        break;
      default:
        break;
    }

    return result;
  }

  /**
   * @brief Value at @p xi of the basis function of node 0.
   */
  constexpr static double value0(const double xi) { return (1.0 / 8.0) * (-1.0 + xi) * xi * (-3.0 + 7.0 * xi * xi); }

  /**
   * @brief Value at @p xi of the basis function of node 1.
   */
  constexpr static double value1(const double xi) { return (49.0 / 24.0) * (sqrt3_7 - xi) * xi * (-1.0 + xi * xi); }

  /**
   * @brief Value at @p xi of the basis function of node 2.
   */
  constexpr static double value2(const double xi) {
    return (1.0 / 3.0) * (3.0 - 10.0 * xi * xi + 7.0 * xi * xi * xi * xi);
  }

  /**
   * @brief Value at @p xi of the basis function of node 3.
   */
  constexpr static double value3(const double xi) { return -(49.0 / 24.0) * (sqrt3_7 + xi) * xi * (-1.0 + xi * xi); }

  /**
   * @brief Value at @p xi of the basis function of node 4.
   */
  constexpr static double value4(const double xi) { return (1.0 / 8.0) * (1.0 + xi) * xi * (-3.0 + 7.0 * xi * xi); }

  /**
   * @brief Derivative at @p xi of the basis function of node @p index.
   */
  constexpr static double gradient(const int index, const double xi) {
    double result = 0.0;

    switch (index) {
      case 0:
        result = LagrangeBasis4GL::gradient0(xi);
        break;
      case 1:
        result = LagrangeBasis4GL::gradient1(xi);
        break;
      case 2:
        result = LagrangeBasis4GL::gradient2(xi);
        break;
      case 3:
        result = LagrangeBasis4GL::gradient3(xi);
        break;
      case 4:
        result = LagrangeBasis4GL::gradient4(xi);
        break;
      default:
        break;
    }

    return result;
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 0.
   */
  constexpr static double gradient0(const double xi) {
    return (1.0 / 8.0) * (3.0 + xi * (-6.0 + 7.0 * xi * (-3.0 + 4.0 * xi)));
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 1.
   */
  constexpr static double gradient1(const double xi) {
    return (49.0 / 24.0) * (-sqrt3_7 + xi * (2.0 + 3.0 * sqrt3_7 * xi - 4.0 * xi * xi));
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 2.
   */
  constexpr static double gradient2(const double xi) { return (4.0 / 3.0) * xi * (-5.0 + 7.0 * xi * xi); }

  /**
   * @brief Derivative at @p xi of the basis function of node 3.
   */
  constexpr static double gradient3(const double xi) {
    return (49.0 / 24.0) * (sqrt3_7 + xi * (2.0 - 3.0 * sqrt3_7 * xi - 4.0 * xi * xi));
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 4.
   */
  constexpr static double gradient4(const double xi) {
    return (1.0 / 8.0) * (-3.0 + xi * (-6.0 + 7.0 * xi * (3.0 + 4.0 * xi)));
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
        switch (p) {
          case 0:
            return -5.0000000000000000000;
          case 1:
            return -1.2409902530309828578;
          case 2:
            return 0.37500000000000000000;
        }
        break;
      case 1:
        switch (p) {
          case 0:
            return 6.7565024887242400038;
          case 1:
            return 0.0;
          case 2:
            return -1.3365845776954533353;
        }
        break;
      case 2:
        switch (p) {
          case 0:
            return -2.6666666666666666667;
          case 1:
            return 1.7457431218879390501;
          case 2:
            return 0.0;
        }
        break;
      case 3:
        switch (p) {
          case 0:
            return 1.4101641779424266628;
          case 1:
            return -0.7637626158259733344;
          case 2:
            return 1.3365845776954533353;
        }
        break;
      case 4:
        switch (p) {
          case 0:
            return -0.50000000000000000000;
          case 1:
            return 0.25900974696901714215;
          case 2:
            return -0.37500000000000000000;
        }
        break;
    }
    return 0;
  }

  /**
   * @brief Tensor product of the 1D basis on the parent square [-1, 1]^2.
   *
   * See docs/design.md, "1D Lagrange bases".
   */
  struct TensorProduct2D {
    constexpr static int numSupportPoints = 25;  ///< Number of nodes, 5^2.

    /**
     * @brief Index i + 5*j of the node (i, j), with i and j in [0, 4].
     */
    constexpr static int linearIndex(const int i, const int j) { return i + 5 * j; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     */
    constexpr static void multiIndex(int const linearIndex, int& i0, int& i1) {
      i1 = linearIndex / 5;

      i0 = linearIndex % 5;
    }

    /**
     * @brief Values at one point of all the 2D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1).
     * @param[out] N N[linearIndex(a, b)] = value(a, xi0) * value(b, xi1).
     */
    static void value(const double (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 5; ++a) {
        for (int b = 0; b < 5; ++b) {
          const int lindex = LagrangeBasis4GL::TensorProduct2D::linearIndex(a, b);
          N[lindex] = LagrangeBasis4GL::value(a, coords[0]) * LagrangeBasis4GL::value(b, coords[1]);
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
    constexpr static int numSupportPoints = 125;  ///< Number of nodes, 5^3.

    /**
     * @brief Index i + 5*j + 25*k of the node (i, j, k), with i, j and k in
     * [0, 4].
     *
     * See docs/design.md, "Hexahedron local numbering".
     */
    constexpr static int linearIndex(const int i, const int j, const int k) { return i + 5 * j + 25 * k; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     * @param[out] i2 Index along xi2.
     */
    constexpr static void multiIndex(int const linearIndex, int& i0, int& i1, int& i2) {
      i2 = linearIndex / 25;

      i1 = (linearIndex % 25) / 5;

      i0 = (linearIndex % 25) % 5;
    }

    /**
     * @brief Values at one point of all the 3D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1, xi2).
     * @param[out] N N[linearIndex(a, b, c)] = value(a, xi0) * value(b, xi1) *
     * value(c, xi2).
     */
    PROXY_HOST_DEVICE
    static void value(const double (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 5; ++a) {
        for (int b = 0; b < 5; ++b) {
          for (int c = 0; c < 5; ++c) {
            const int lindex = LagrangeBasis4GL::TensorProduct3D::linearIndex(a, b, c);
            N[lindex] = LagrangeBasis4GL::value(a, coords[0]) * LagrangeBasis4GL::value(b, coords[1]) *
                        LagrangeBasis4GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif /* _LAGRANGEBASIS4GL_HPP_ */
