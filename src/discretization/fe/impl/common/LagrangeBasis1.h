#ifndef _LAGRANGEBASIS1_HPP_
#define _LAGRANGEBASIS1_HPP_

/**
 * @brief Linear (order 1) Lagrange basis on the nodes -1 and 1 of [-1, 1].
 *
 * See docs/design.md, "1D Lagrange bases".
 */
class LagrangeBasis1 {
 public:
  constexpr static int numSupportPoints = 2;  ///< Number of nodes, order + 1.

  /**
   * @brief Gauss-Lobatto quadrature weight of node @p q on [-1, 1].
   */
  constexpr static double weight(const int q) { return 1.0; }

  /**
   * @brief Parent coordinate of node @p supportPointIndex, in [-1, 1].
   */
  constexpr static double parentSupportCoord(const int supportPointIndex) {
    return -1.0 + 2.0 * (supportPointIndex & 1);
  }

  /**
   * @brief Value at @p xi of the basis function of node @p index.
   */
  constexpr static double value(const int index, const double xi) { return 0.5 + 0.5 * xi * parentSupportCoord(index); }

  /**
   * @brief Value at @p xi of the basis function of node 0.
   */
  constexpr static double value0(const double xi) { return 0.5 - 0.5 * xi; }

  /**
   * @brief Value at @p xi of the basis function of node 1.
   */
  constexpr static double value1(const double xi) { return 0.5 + 0.5 * xi; }

  /**
   * @brief Bubble function 1 - xi^2, which vanishes at both nodes.
   */
  constexpr static double valueBubble(const double xi) { return 1.0 - pow(xi, 2); }

  /**
   * @brief Derivative at @p xi of the basis function of node @p index.
   */
  constexpr static double gradient(const int index, const double xi) { return 0.5 * parentSupportCoord(index); }

  /**
   * @brief Derivative at @p xi of the basis function of node 0.
   */
  constexpr static double gradient0(const double xi) { return -0.5; }

  /**
   * @brief Derivative at @p xi of the basis function of node 1.
   */
  constexpr static double gradient1(const double xi) { return 0.5; }

  /**
   * @brief Derivative of the basis function of node @p q at any node; it is
   * constant.
   */
  constexpr static double gradientAt(const int q, const int) { return q == 0 ? -0.5 : 0.5; }

  /**
   * @brief Tensor product of the 1D basis on the parent square [-1, 1]^2.
   *
   * See docs/design.md, "1D Lagrange bases".
   */
  struct TensorProduct2D {
    constexpr static int numSupportPoints = 4;  ///< Number of nodes, 2^2.

    /**
     * @brief Index i + 2*j of the node (i, j), with i and j in [0, 1].
     */
    constexpr static int linearIndex(const int i, const int j) { return i + 2 * j; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     */
    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1) {
      i0 = (linearIndex & 1);
      i1 = (linearIndex & 2) >> 1;
    }

    /**
     * @brief Values at one point of all the 2D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1).
     * @param[out] N N[linearIndex(a, b)] = value(a, xi0) * value(b, xi1).
     */
    static void value(double const (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
          const int lindex = LagrangeBasis1::TensorProduct2D::linearIndex(a, b);
          N[lindex] = LagrangeBasis1::value(a, coords[0]) * LagrangeBasis1::value(b, coords[1]);
        }
      }
    }

    /**
     * @brief Parent coordinate along xi0 of node @p linearIndex.
     */
    constexpr static double parentCoords0(int const linearIndex) { return -1.0 + 2.0 * (linearIndex & 1); }

    /**
     * @brief Parent coordinate along xi1 of node @p linearIndex.
     */
    constexpr static double parentCoords1(int const linearIndex) { return -1.0 + (linearIndex & 2); }
  };

  /**
   * @brief Tensor product of the 1D basis on the parent cube [-1, 1]^3.
   *
   * See docs/design.md, "1D Lagrange bases".
   */
  struct TensorProduct3D {
    constexpr static int numSupportPoints = 8;  ///< Number of nodes, 2^3.

    constexpr static int numSupportFaces = 6;  ///< Number of faces of the cube.

    /**
     * @brief Index i + 2*j + 4*k of the node (i, j, k), with i, j and k in
     * [0, 1].
     *
     * See docs/design.md, "Hexahedron local numbering".
     */
    constexpr static int linearIndex(const int i, const int j, const int k) { return i + 2 * j + 4 * k; }

    /**
     * @brief Inverse of linearIndex().
     * @param[in] linearIndex Node index.
     * @param[out] i0 Index along xi0.
     * @param[out] i1 Index along xi1.
     * @param[out] i2 Index along xi2.
     */
    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1, int& i2) {
      i0 = (linearIndex & 1);
      i1 = (linearIndex & 2) >> 1;
      i2 = (linearIndex & 4) >> 2;
    }

    /**
     * @brief Values at one point of all the 3D basis functions.
     * @param[in] coords Parent coordinates (xi0, xi1, xi2).
     * @param[out] N N[linearIndex(a, b, c)] = value(a, xi0) * value(b, xi1) *
     * value(c, xi2).
     */
    PROXY_HOST_DEVICE
    static void value(double const (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
          for (int c = 0; c < 2; ++c) {
            const int lindex = LagrangeBasis1::TensorProduct3D::linearIndex(a, b, c);
            N[lindex] = LagrangeBasis1::value(a, coords[0]) * LagrangeBasis1::value(b, coords[1]) *
                        LagrangeBasis1::value(c, coords[2]);
          }
        }
      }
    }

    /**
     * @brief Parent coordinate along xi0 of node @p linearIndex.
     */
    constexpr static double parentCoords0(int const linearIndex) { return -1.0 + 2.0 * (linearIndex & 1); }

    /**
     * @brief Parent coordinate along xi1 of node @p linearIndex.
     */
    constexpr static double parentCoords1(int const linearIndex) { return -1.0 + (linearIndex & 2); }

    /**
     * @brief Parent coordinate along xi2 of node @p linearIndex.
     */
    constexpr static double parentCoords2(int const linearIndex) { return -1.0 + 0.5 * (linearIndex & 4); }
  };
};

#endif /* _LAGRANGEBASIS1_HPP_ */
