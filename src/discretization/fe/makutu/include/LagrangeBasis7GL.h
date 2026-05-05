#ifndef _LAGRANGEBASIS7GL_HPP_
#define _LAGRANGEBASIS7GL_HPP_

/**
 * @file LagrangeBasis7GL.hpp
 *
 * 1D Lagrange basis on Gauss-Lobatto-Legendre nodes for degree 7 (8 nodes).
 * Follows the same implementation pattern as LagrangeBasis6GL.hpp.
 * Polynomial formulas and gradientAt table generated via computGL.py.
 *
 * GLL nodes: { -1, -lambda1, -lambda2, -lambda3, lambda3, lambda2, lambda1, 1 }
 *   lambda1 = 0.8717401485096067
 *   lambda2 = 0.5917001814331423
 *   lambda3 = 0.2092992179024789
 */
class LagrangeBasis7GL {
 public:
  constexpr static int numSupportPoints = 8;

  /// |nodes[6]| = |nodes[1]|
  static constexpr double lambda1 = 0.8717401485096067;

  /// |nodes[5]| = |nodes[2]|
  static constexpr double lambda2 = 0.5917001814331423;

  /// |nodes[4]| = |nodes[3]|
  static constexpr double lambda3 = 0.2092992179024789;

  /**
   * @brief GLL quadrature weight for node q.
   */
  constexpr static double weight(const int q) {
    switch (q) {
      case 1:
      case 6:
        return 0.2107042271435060;
      case 2:
      case 5:
        return 0.3411226924835044;
      case 3:
      case 4:
        return 0.4124587946587037;
      default:
        return 1.0 / 28.0;
    }
  }

  /**
   * @brief Parent coordinate of support point @p supportPointIndex.
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
        return -lambda3;
      case 4:
        return lambda3;
      case 5:
        return lambda2;
      case 6:
        return lambda1;
      case 7:
        return 1.0;
      default:
        return 0.0;
    }
  }

  /**
   * @brief Value of basis function @p index at @p xi.
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
      case 7:
        return value7(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Basis function for support point 0.
   */
  constexpr static double value0(const double xi) {
    return -3.3515625 * xi * xi * xi * xi * xi * xi * xi + 3.3515625 * xi * xi * xi * xi * xi * xi +
           3.8671875 * xi * xi * xi * xi * xi - 3.8671875 * xi * xi * xi * xi - 1.0546875 * xi * xi * xi +
           1.0546875 * xi * xi + 0.0390625 * xi - 0.0390625;
  }

  /**
   * @brief Basis function for support point 1.
   */
  constexpr static double value1(const double xi) {
    return 8.14072271825387 * xi * xi * xi * xi * xi * xi * xi - 7.09659483138615 * xi * xi * xi * xi * xi * xi -
           11.347477684014 * xi * xi * xi * xi * xi + 9.89205188147183 * xi * xi * xi * xi +
           3.33160871212585 * xi * xi * xi - 2.90429707348449 * xi * xi - 0.124853746365692 * xi + 0.108840023398809;
  }

  /**
   * @brief Basis function for support point 2.
   */
  constexpr static double value2(const double xi) {
    return -10.3581368289505 * xi * xi * xi * xi * xi * xi * xi + 6.12891144099929 * xi * xi * xi * xi * xi * xi +
           18.6833551584202 * xi * xi * xi * xi * xi - 11.0549446370171 * xi * xi * xi * xi -
           8.67003714121216 * xi * xi * xi + 5.13006254948732 * xi * xi + 0.34481881174243 * xi - 0.204029353469556;
  }

  /**
   * @brief Basis function for support point 3.
   */
  constexpr static double value3(const double xi) {
    return 11.3898137484866 * xi * xi * xi * xi * xi * xi * xi - 2.38387910961314 * xi * xi * xi * xi * xi * xi -
           24.0329625019858 * xi * xi * xi * xi * xi + 5.03008025554523 * xi * xi * xi * xi +
           15.6735080468926 * xi * xi * xi - 3.28045297600283 * xi * xi - 3.0303592933934 * xi + 0.634251830070747;
  }

  /**
   * @brief Basis function for support point 4.
   */
  constexpr static double value4(const double xi) {
    return -11.3898137484866 * xi * xi * xi * xi * xi * xi * xi - 2.38387910961314 * xi * xi * xi * xi * xi * xi +
           24.0329625019858 * xi * xi * xi * xi * xi + 5.03008025554523 * xi * xi * xi * xi -
           15.6735080468926 * xi * xi * xi - 3.28045297600283 * xi * xi + 3.0303592933934 * xi + 0.634251830070747;
  }

  /**
   * @brief Basis function for support point 5.
   */
  constexpr static double value5(const double xi) {
    return 10.3581368289505 * xi * xi * xi * xi * xi * xi * xi + 6.12891144099929 * xi * xi * xi * xi * xi * xi -
           18.6833551584202 * xi * xi * xi * xi * xi - 11.0549446370171 * xi * xi * xi * xi +
           8.67003714121216 * xi * xi * xi + 5.13006254948732 * xi * xi - 0.34481881174243 * xi - 0.204029353469556;
  }

  /**
   * @brief Basis function for support point 6.
   */
  constexpr static double value6(const double xi) {
    return -8.14072271825387 * xi * xi * xi * xi * xi * xi * xi - 7.09659483138615 * xi * xi * xi * xi * xi * xi +
           11.347477684014 * xi * xi * xi * xi * xi + 9.89205188147183 * xi * xi * xi * xi -
           3.33160871212585 * xi * xi * xi - 2.90429707348449 * xi * xi + 0.124853746365692 * xi + 0.108840023398809;
  }

  /**
   * @brief Basis function for support point 7.
   */
  constexpr static double value7(const double xi) {
    return 3.3515625 * xi * xi * xi * xi * xi * xi * xi + 3.3515625 * xi * xi * xi * xi * xi * xi -
           3.8671875 * xi * xi * xi * xi * xi - 3.8671875 * xi * xi * xi * xi + 1.0546875 * xi * xi * xi +
           1.0546875 * xi * xi - 0.0390625 * xi - 0.0390625;
  }

  /**
   * @brief Gradient of basis function @p index at @p xi.
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
      case 7:
        return gradient7(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Gradient of basis function for support point 0.
   */
  constexpr static double gradient0(const double xi) {
    return -23.4609375 * xi * xi * xi * xi * xi * xi + 20.109375 * xi * xi * xi * xi * xi +
           19.3359375 * xi * xi * xi * xi - 15.46875 * xi * xi * xi - 3.1640625 * xi * xi + 2.109375 * xi + 0.0390625;
  }

  /**
   * @brief Gradient of basis function for support point 1.
   */
  constexpr static double gradient1(const double xi) {
    return 56.9850590277771 * xi * xi * xi * xi * xi * xi - 42.5795689883169 * xi * xi * xi * xi * xi -
           56.7373884200701 * xi * xi * xi * xi + 39.5682075258873 * xi * xi * xi + 9.99482613637756 * xi * xi -
           5.80859414696897 * xi - 0.124853746365692;
  }

  /**
   * @brief Gradient of basis function for support point 2.
   */
  constexpr static double gradient2(const double xi) {
    return -72.5069578026532 * xi * xi * xi * xi * xi * xi + 36.7734686459958 * xi * xi * xi * xi * xi +
           93.4167757921009 * xi * xi * xi * xi - 44.2197785480682 * xi * xi * xi - 26.0101114236365 * xi * xi +
           10.2601250989746 * xi + 0.34481881174243;
  }

  /**
   * @brief Gradient of basis function for support point 3.
   */
  constexpr static double gradient3(const double xi) {
    return 79.7286962394062 * xi * xi * xi * xi * xi * xi - 14.3032746576789 * xi * xi * xi * xi * xi -
           120.164812509929 * xi * xi * xi * xi + 20.1203210221809 * xi * xi * xi + 47.0205241406778 * xi * xi -
           6.56090595200567 * xi - 3.0303592933934;
  }

  /**
   * @brief Gradient of basis function for support point 4.
   */
  constexpr static double gradient4(const double xi) {
    return -79.7286962394062 * xi * xi * xi * xi * xi * xi - 14.3032746576789 * xi * xi * xi * xi * xi +
           120.164812509929 * xi * xi * xi * xi + 20.1203210221809 * xi * xi * xi - 47.0205241406778 * xi * xi -
           6.56090595200567 * xi + 3.0303592933934;
  }

  /**
   * @brief Gradient of basis function for support point 5.
   */
  constexpr static double gradient5(const double xi) {
    return 72.5069578026532 * xi * xi * xi * xi * xi * xi + 36.7734686459958 * xi * xi * xi * xi * xi -
           93.4167757921010 * xi * xi * xi * xi - 44.2197785480682 * xi * xi * xi + 26.0101114236365 * xi * xi +
           10.2601250989746 * xi - 0.34481881174243;
  }

  /**
   * @brief Gradient of basis function for support point 6.
   */
  constexpr static double gradient6(const double xi) {
    return -56.9850590277771 * xi * xi * xi * xi * xi * xi - 42.5795689883169 * xi * xi * xi * xi * xi +
           56.7373884200701 * xi * xi * xi * xi + 39.5682075258873 * xi * xi * xi - 9.99482613637755 * xi * xi -
           5.80859414696897 * xi + 0.124853746365692;
  }

  /**
   * @brief Gradient of basis function for support point 7.
   */
  constexpr static double gradient7(const double xi) {
    return 23.4609375 * xi * xi * xi * xi * xi * xi + 20.109375 * xi * xi * xi * xi * xi -
           19.3359375 * xi * xi * xi * xi - 15.46875 * xi * xi * xi + 3.1640625 * xi * xi + 2.109375 * xi - 0.0390625;
  }

  /**
   * @brief Gradient of basis function q evaluated at support point p.
   * Full 8x8 table computed via computGL.py.
   * Anti-symmetry: gradientAt(q,p) == -gradientAt(7-q, 7-p).
   */
  constexpr static double gradientAt(const int q, const int p) {
    switch (q) {
      case 0:
        switch (p) {
          case 0:
            return -14.0000000000000000;
          case 1:
            return -3.2099157030029941;
          case 2:
            return 0.7924766813205160;
          case 3:
            return -0.3721504357285950;
          case 4:
            return 0.2433307127237914;
          case 5:
            return -0.2032845689005931;
          case 6:
            return 0.2199575147712993;
          case 7:
            return -0.5000000000000000;
        }
        break;
      case 1:
        switch (p) {
          case 0:
            return 18.9375986071174012;
          case 1:
            return 0.0;
          case 2:
            return -2.8064757947364267;
          case 3:
            return 1.0789446887904528;
          case 4:
            return -0.6611573509003108;
          case 5:
            return 0.5370395861576664;
          case 6:
            return -0.5735654149402336;
          case 7:
            return 1.2976873883202700;
        }
        break;
      case 2:
        switch (p) {
          case 0:
            return -7.5692898193484837;
          case 1:
            return 4.5435850645665532;
          case 2:
            return 0.0;
          case 3:
            return -2.3781872335155048;
          case 4:
            return 1.1353580168811113;
          case 5:
            return -0.8450225565065095;
          case 6:
            return 0.8694480983314570;
          case 7:
            return -1.9416594255441311;
        }
        break;
      case 3:
        switch (p) {
          case 0:
            return 4.2979081642651522;
          case 1:
            return -2.1120612143145365;
          case 2:
            return 2.8755174059725039;
          case 3:
            return 0.0;
          case 4:
            return -2.3889243591582399;
          case 5:
            return 1.3727858318060253;
          case 6:
            return -1.2942320509134930;
          case 7:
            return 2.8101889892579379;
        }
        break;
      case 4:
        switch (p) {
          case 0:
            return -2.8101889892579237;
          case 1:
            return 1.2942320509134930;
          case 2:
            return -1.3727858318060253;
          case 3:
            return 2.3889243591582399;
          case 4:
            return 0.0;
          case 5:
            return -2.8755174059725039;
          case 6:
            return 2.1120612143145365;
          case 7:
            return -4.2979081642651664;
        }
        break;
      case 5:
        switch (p) {
          case 0:
            return 1.9416594255441169;
          case 1:
            return -0.8694480983314818;
          case 2:
            return 0.8450225565065104;
          case 3:
            return -1.1353580168811113;
          case 4:
            return 2.3781872335155052;
          case 5:
            return 0.0;
          case 6:
            return -4.5435850645665781;
          case 7:
            return 7.5692898193484552;
        }
        break;
      case 6:
        switch (p) {
          case 0:
            return -1.2976873883202629;
          case 1:
            return 0.5735654149402336;
          case 2:
            return -0.5370395861576669;
          case 3:
            return 0.6611573509003108;
          case 4:
            return -1.0789446887904528;
          case 5:
            return 2.8064757947364272;
          case 6:
            return 0.0;
          case 7:
            return -18.9375986071174012;
        }
        break;
      case 7:
        switch (p) {
          case 0:
            return 0.5000000000000000;
          case 1:
            return -0.2199575147713002;
          case 2:
            return 0.2032845689005938;
          case 3:
            return -0.2433307127237914;
          case 4:
            return 0.3721504357285950;
          case 5:
            return -0.7924766813205157;
          case 6:
            return 3.2099157030029950;
          case 7:
            return 14.0000000000000000;
        }
        break;
    }
    return 0.0;
  }

  /* Tensor product helpers (2D / 3D) */
  struct TensorProduct2D {
    constexpr static int numSupportPoints1D = LagrangeBasis7GL::numSupportPoints;
    constexpr static int numSupportPoints = numSupportPoints1D * numSupportPoints1D;  // 8*8 = 64

    constexpr static int linearIndex(const int i, const int j) { return i + numSupportPoints1D * j; }

    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1) {
      i1 = linearIndex / numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    static void value(const double (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < LagrangeBasis7GL::numSupportPoints; ++a) {
        for (int b = 0; b < LagrangeBasis7GL::numSupportPoints; ++b) {
          const int lindex = a + LagrangeBasis7GL::numSupportPoints * b;
          N[lindex] = LagrangeBasis7GL::value(a, coords[0]) * LagrangeBasis7GL::value(b, coords[1]);
        }
      }
    }
  };

  struct TensorProduct3D {
    constexpr static int numSupportPoints1D = LagrangeBasis7GL::numSupportPoints;
    constexpr static int numSupportPoints = numSupportPoints1D * numSupportPoints1D * numSupportPoints1D;  // 8^3 = 512

    constexpr static int linearIndex(const int i, const int j, const int k) {
      return i + numSupportPoints1D * j + numSupportPoints1D * numSupportPoints1D * k;
    }

    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1, int& i2) {
      i2 = linearIndex / (numSupportPoints1D * numSupportPoints1D);
      i1 = (linearIndex % (numSupportPoints1D * numSupportPoints1D)) / numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    static void value(const double (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < LagrangeBasis7GL::numSupportPoints; ++a) {
        for (int b = 0; b < LagrangeBasis7GL::numSupportPoints; ++b) {
          for (int c = 0; c < LagrangeBasis7GL::numSupportPoints; ++c) {
            const int lindex = a + LagrangeBasis7GL::numSupportPoints * b +
                               LagrangeBasis7GL::numSupportPoints * LagrangeBasis7GL::numSupportPoints * c;
            N[lindex] = LagrangeBasis7GL::value(a, coords[0]) * LagrangeBasis7GL::value(b, coords[1]) *
                        LagrangeBasis7GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif /* _LAGRANGEBASIS7GL_HPP_ */
