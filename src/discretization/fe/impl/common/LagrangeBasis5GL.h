#ifndef SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS5GL_H_
#define SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS5GL_H_

/**
 * @brief Order 5 Lagrange basis on the 6 Gauss-Lobatto-Legendre nodes -1,
 * -sqrt((7 + 2 sqrt(7))/21), -sqrt((7 - 2 sqrt(7))/21), their opposites, and 1.
 *
 * See docs/design.md, "1D Lagrange bases".
 */
class LagrangeBasis5GL {
 public:
  constexpr static int numSupportPoints = 6;  ///< Number of nodes, order + 1.

  static constexpr double sqrt_7_ = 2.64575131106459059;  ///< sqrt(7).

  static constexpr double sqrt__7_plus_2sqrt7__ = 3.50592393273573196;  ///< sqrt(7 + 2 sqrt(7)).

  static constexpr double sqrt__7_mins_2sqrt7__ = 1.30709501485960033;  ///< sqrt(7 - 2 sqrt(7)).

  static constexpr double sqrt__7_plus_sqrt7_div2__ = 2.884939454396278;  ///< sqrt(7 + sqrt(7)/2).

  static constexpr double sqrt__7_mins_sqrt7_div2__ = 2.382671682055189;  ///< sqrt(7 - sqrt(7)/2).

  static constexpr double sqrt_inv21 = 0.218217890235992381;  ///< sqrt(1/21).

  /**
   * @brief Gauss-Lobatto quadrature weight of node @p q on [-1, 1].
   */
  constexpr static double weight(const int q) {
    switch (q) {
      case 1:
      case 4:
        return (1.0 / 30.0) * (14.0 - sqrt_7_);
      case 2:
      case 3:
        return (1.0 / 30.0) * (14.0 + sqrt_7_);
      default:
        return 1.0 / 15.0;
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
        result = -sqrt_inv21 * sqrt__7_plus_2sqrt7__;
        break;

      case 2:
        result = -sqrt_inv21 * sqrt__7_mins_2sqrt7__;
        break;

      case 3:
        result = sqrt_inv21 * sqrt__7_mins_2sqrt7__;
        break;

      case 4:
        result = sqrt_inv21 * sqrt__7_plus_2sqrt7__;
        break;

      case 5:
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
        return result = LagrangeBasis5GL::value0(xi);
        break;

      case 1:
        return result = LagrangeBasis5GL::value1(xi);
        break;

      case 2:
        return result = LagrangeBasis5GL::value2(xi);
        break;

      case 3:
        return result = LagrangeBasis5GL::value3(xi);
        break;

      case 4:
        return result = LagrangeBasis5GL::value4(xi);
        break;

      case 5:
        return result = LagrangeBasis5GL::value5(xi);
        break;

      default:
        break;
    }

    return result;
  }

  /**
   * @brief Value at @p xi of the basis function of node 0.
   */
  constexpr static double value0(const double xi) {

    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return (-21.0 / 16.0) *
           (xi * xi * xi * xi * xi - xi * xi * xi * xi - (lambda3 * lambda3 + lambda4 * lambda4) * xi * xi * xi +
            (lambda3 * lambda3 + lambda4 * lambda4) * xi * xi + lambda3 * lambda3 * lambda4 * lambda4 * xi -
            lambda3 * lambda3 * lambda4 * lambda4);
  }

  /**
   * @brief Value at @p xi of the basis function of node 1.
   */
  constexpr static double value1(const double xi) {

    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.

    return ((21.0 / 16.0) * sqrt__7_mins_sqrt7_div2__) *
           (xi * xi * xi * xi * xi - lambda4 * xi * xi * xi * xi - (lambda3 * lambda3 + 1) * xi * xi * xi +
            lambda4 * (lambda3 * lambda3 + 1) * xi * xi + lambda3 * lambda3 * xi - lambda4 * lambda3 * lambda3);
  }

  /**
   * @brief Value at @p xi of the basis function of node 2.
   */
  constexpr static double value2(const double xi) {

    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return ((-21.0 / 16.0) * sqrt__7_plus_sqrt7_div2__) *
           (xi * xi * xi * xi * xi - lambda3 * xi * xi * xi * xi - (lambda4 * lambda4 + 1) * xi * xi * xi +
            lambda3 * (lambda4 * lambda4 + 1) * xi * xi + lambda4 * lambda4 * xi - lambda3 * lambda4 * lambda4);
  }

  /**
   * @brief Value at @p xi of the basis function of node 3.
   */
  constexpr static double value3(const double xi) {

    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return ((21.0 / 16.0) * sqrt__7_plus_sqrt7_div2__) *
           (xi * xi * xi * xi * xi + lambda3 * xi * xi * xi * xi - (lambda4 * lambda4 + 1) * xi * xi * xi -
            lambda3 * (lambda4 * lambda4 + 1) * xi * xi + lambda4 * lambda4 * xi + lambda3 * lambda4 * lambda4);
  }

  /**
   * @brief Value at @p xi of the basis function of node 4.
   */
  constexpr static double value4(const double xi) {

    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return ((-21.0 / 16.0) * sqrt__7_mins_sqrt7_div2__) *
           (xi * xi * xi * xi * xi + lambda4 * xi * xi * xi * xi - (lambda3 * lambda3 + 1) * xi * xi * xi -
            lambda4 * (lambda3 * lambda3 + 1) * xi * xi + lambda3 * lambda3 * xi + lambda4 * lambda3 * lambda3);
  }

  /**
   * @brief Value at @p xi of the basis function of node 5.
   */
  constexpr static double value5(const double xi) {

    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.

    return (21.0 / 16.0) *
           (xi * xi * xi * xi * xi + xi * xi * xi * xi - (lambda4 * lambda4 + lambda3 * lambda3) * xi * xi * xi -
            (lambda3 * lambda3 + lambda4 * lambda4) * xi * xi + lambda3 * lambda3 * lambda4 * lambda4 * xi +
            lambda4 * lambda4 * lambda3 * lambda3);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node @p index.
   */
  constexpr static double gradient(const int index, const double xi) {
    double result = 0.0;

    switch (index) {
      case 0:
        result = LagrangeBasis5GL::gradient0(xi);
        break;

      case 1:
        result = LagrangeBasis5GL::gradient1(xi);
        break;

      case 2:
        result = LagrangeBasis5GL::gradient2(xi);
        break;

      case 3:
        result = LagrangeBasis5GL::gradient3(xi);
        break;

      case 4:
        result = LagrangeBasis5GL::gradient4(xi);
        break;

      case 5:
        result = LagrangeBasis5GL::gradient5(xi);
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
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return (-21.0 / 16.0) *
           (5.0 * xi * xi * xi * xi - 4.0 * xi * xi * xi - 3.0 * (lambda3 * lambda3 + lambda4 * lambda4) * xi * xi +
            2.0 * (lambda3 * lambda3 + lambda4 * lambda4) * xi + lambda3 * lambda3 * lambda4 * lambda4);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 1.
   */
  constexpr static double gradient1(const double xi) {
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.

    return (21.0 / 16.0) * sqrt__7_mins_sqrt7_div2__ *
           (5.0 * xi * xi * xi * xi - 4.0 * lambda4 * xi * xi * xi - 3.0 * (lambda3 * lambda3 + 1.0) * xi * xi +
            2.0 * lambda4 * (lambda3 * lambda3 + 1.0) * xi + lambda3 * lambda3);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 2.
   */
  constexpr static double gradient2(const double xi) {
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return (-21.0 / 16.0) * sqrt__7_plus_sqrt7_div2__ *
           (5.0 * xi * xi * xi * xi - 4.0 * lambda3 * xi * xi * xi - 3.0 * (lambda4 * lambda4 + 1.0) * xi * xi +
            2.0 * lambda3 * (lambda4 * lambda4 + 1.0) * xi + lambda4 * lambda4);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 3.
   */
  constexpr static double gradient3(const double xi) {
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return (21.0 / 16.0) * sqrt__7_plus_sqrt7_div2__ *
           (5.0 * xi * xi * xi * xi + 4.0 * lambda3 * xi * xi * xi - 3.0 * (lambda4 * lambda4 + 1.0) * xi * xi -
            2 * lambda3 * (lambda4 * lambda4 + 1.0) * xi + lambda4 * lambda4);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 4.
   */
  constexpr static double gradient4(const double xi) {
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return (-21.0 / 16.0) * sqrt__7_mins_sqrt7_div2__ *
           (5.0 * xi * xi * xi * xi + 4.0 * lambda4 * xi * xi * xi - 3.0 * (lambda3 * lambda3 + 1.0) * xi * xi -
            2.0 * lambda4 * (lambda3 * lambda3 + 1.0) * xi + lambda3 * lambda3);
  }

  /**
   * @brief Derivative at @p xi of the basis function of node 5.
   */
  constexpr static double gradient5(const double xi) {
    double lambda4 = LagrangeBasis5GL::parentSupportCoord(4);  ///< Parent coordinate of node 1, opposite of node 4.
    double lambda3 = LagrangeBasis5GL::parentSupportCoord(3);  ///< Parent coordinate of node 2, opposite of node 3.

    return (21.0 / 16.0) *
           (5.0 * xi * xi * xi * xi + 4.0 * xi * xi * xi - 3.0 * (lambda3 * lambda3 + lambda4 * lambda4) * xi * xi -
            2.0 * (lambda3 * lambda3 + lambda4 * lambda4) * xi + lambda3 * lambda3 * lambda4 * lambda4);
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
            return -7.5000000000000000000;
          case 1:
            return -1.7863649483390948939;
          case 2:
            return 0.48495104785356916930;
        }
        break;
      case 1:
        switch (p) {
          case 0:
            return 10.14141593631966928023;
          case 1:
            return 0.0;
          case 2:
            return -1.72125695283023338321;
        }
        break;
      case 2:
        switch (p) {
          case 0:
            return -4.03618727030534800527;
          case 1:
            return 2.5234267774294554319088;
          case 2:
            return 0.0;
        }
        break;
      case 3:
        switch (p) {
          case 0:
            return 2.2446846481761668242712;
          case 1:
            return -1.1528281585359293413318;
          case 2:
            return 1.7529619663678659788775;
        }
        break;
      case 4:
        switch (p) {
          case 0:
            return -1.3499133141904880992312;
          case 1:
            return 0.6535475074298001672007;
          case 2:
            return -0.7863566722232407374395;
        }
        break;
      case 5:
        switch (p) {
          case 0:
            return 0.500000000000000000000;
          case 1:
            return -0.2377811779842313638052;
          case 2:
            return 0.2697006108320389724720;
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
    constexpr static int numSupportPoints = 36;  ///< Number of nodes, 6^2.

    /**
     * @brief Index i + 6*j of the node (i, j), with i and j in [0, 5].
     */
    constexpr static int linearIndex(const int i, const int j) { return i + 6 * j; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     */
    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1) {
      i1 = linearIndex / 6;

      i0 = linearIndex % 6;
    }

    /**
     * @brief Values at one point of all the 2D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1).
     * @param[out] N N[linearIndex(a, b)] = value(a, xi0) * value(b, xi1).
     */
    static void value(const double (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 6; ++a) {
        for (int b = 0; b < 6; ++b) {
          const int lindex = LagrangeBasis5GL::TensorProduct2D::linearIndex(a, b);
          N[lindex] = LagrangeBasis5GL::value(a, coords[0]) * LagrangeBasis5GL::value(b, coords[1]);
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
    constexpr static int numSupportPoints = 216;  ///< Number of nodes, 6^3.

    /**
     * @brief Index i + 6*j + 36*k of the node (i, j, k), with i, j and k in
     * [0, 5].
     *
     * See docs/design.md, "Hexahedron local numbering".
     */
    constexpr static int linearIndex(const int i, const int j, const int k) { return i + 6 * j + 36 * k; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     * @param[out] i2 Index along xi2.
     */
    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1, int& i2) {
      i2 = linearIndex / 36;

      i1 = (linearIndex % 36) / 6;

      i0 = (linearIndex % 36) % 6;
    }

    /**
     * @brief Values at one point of all the 3D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1, xi2).
     * @param[out] N N[linearIndex(a, b, c)] = value(a, xi0) * value(b, xi1) *
     * value(c, xi2).
     */
    PROXY_HOST_DEVICE
    static void value(const double (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 6; ++a) {
        for (int b = 0; b < 6; ++b) {
          for (int c = 0; c < 6; ++c) {
            const int lindex = LagrangeBasis5GL::TensorProduct3D::linearIndex(a, b, c);
            N[lindex] = LagrangeBasis5GL::value(a, coords[0]) * LagrangeBasis5GL::value(b, coords[1]) *
                        LagrangeBasis5GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif  // SRC__DISCRETIZATION_FE_MAKUTU_INCLUDE_LAGRANGEBASIS5GL_H_
