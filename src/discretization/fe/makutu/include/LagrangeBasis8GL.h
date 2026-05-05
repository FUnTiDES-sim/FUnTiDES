#ifndef _LAGRANGEBASIS8GL_HPP_
#define _LAGRANGEBASIS8GL_HPP_

/**
 * @file LagrangeBasis8GL.hpp
 *
 * 1D Lagrange basis on Gauss-Lobatto-Legendre nodes for degree 8 (9 nodes).
 * Follows the same implementation pattern as LagrangeBasis7GL.hpp.
 * Polynomial formulas and gradientAt table generated via computGL.py.
 *
 * GLL nodes: { -1, -lambda1, -lambda2, -lambda3, 0, lambda3, lambda2, lambda1,
 * 1 } lambda1 = 0.8997579954114602 lambda2 = 0.6771862795107377 lambda3 =
 * 0.3631174638261782
 */
class LagrangeBasis8GL {
 public:
  constexpr static int numSupportPoints = 9;

  /// |nodes[7]| = |nodes[1]|
  static constexpr double lambda1 = 0.8997579954114602;

  /// |nodes[6]| = |nodes[2]|
  static constexpr double lambda2 = 0.6771862795107377;

  /// |nodes[5]| = |nodes[3]|
  static constexpr double lambda3 = 0.3631174638261782;

  /**
   * @brief GLL quadrature weight for node q.
   */
  constexpr static double weight(const int q) {
    switch (q) {
      case 1:
      case 7:
        return 0.1654953615608056;
      case 2:
      case 6:
        return 0.2745387125001617;
      case 3:
      case 5:
        return 0.3464285109730462;
      case 4:
        return 0.3715192743764172;
      default:
        return 1.0 / 36.0;
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
        return 0.0;
      case 5:
        return lambda3;
      case 6:
        return lambda2;
      case 7:
        return lambda1;
      case 8:
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
      case 8:
        return value8(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Basis function for support point 0.
   */
  constexpr static double value0(const double xi) {
    return 5.5859375 * xi * xi * xi * xi * xi * xi * xi * xi - 5.5859375 * xi * xi * xi * xi * xi * xi * xi -
           7.8203125 * xi * xi * xi * xi * xi * xi + 7.8203125 * xi * xi * xi * xi * xi +
           3.0078125 * xi * xi * xi * xi - 3.0078125 * xi * xi * xi - 0.2734375 * xi * xi + 0.2734375 * xi;
  }

  /**
   * @brief Basis function for support point 1.
   */
  constexpr static double value1(const double xi) {
    return -13.6345320004901 * xi * xi * xi * xi * xi * xi * xi * xi +
           12.2677791811344 * xi * xi * xi * xi * xi * xi * xi + 21.6848443970084 * xi * xi * xi * xi * xi * xi -
           19.5111121254617 * xi * xi * xi * xi * xi - 8.87473674361952 * xi * xi * xi * xi +
           7.98511534224353 * xi * xi * xi + 0.824424347101285 * xi * xi - 0.741782397916255 * xi;
  }

  /**
   * @brief Basis function for support point 2.
   */
  constexpr static double value2(const double xi) {
    return 17.5609949844537 * xi * xi * xi * xi * xi * xi * xi * xi -
           11.8920648580289 * xi * xi * xi * xi * xi * xi * xi - 34.0932448057799 * xi * xi * xi * xi * xi * xi +
           23.0874776064749 * xi * xi * xi * xi * xi + 18.4067902908633 * xi * xi * xi * xi -
           12.4648258348041 * xi * xi * xi - 1.87454046953711 * xi * xi + 1.26941308635815 * xi;
  }

  /**
   * @brief Basis function for support point 3.
   */
  constexpr static double value3(const double xi) {
    return -19.7266861982493 * xi * xi * xi * xi * xi * xi * xi * xi +
           7.16310426200316 * xi * xi * xi * xi * xi * xi * xi + 44.7429986230572 * xi * xi * xi * xi * xi * xi -
           16.2469641839827 * xi * xi * xi * xi * xi - 32.3398660472438 * xi * xi * xi * xi +
           11.7431701395535 * xi * xi * xi + 7.32355362243583 * xi * xi - 2.65931021757392 * xi;
  }

  /**
   * @brief Basis function for support point 4 (centre node, even polynomial).
   */
  constexpr static double value4(const double xi) {
    return 20.4285714285714 * xi * xi * xi * xi * xi * xi * xi * xi - 49.0285714285714 * xi * xi * xi * xi * xi * xi +
           39.6 * xi * xi * xi * xi - 12.0 * xi * xi + 1.0;
  }

  /**
   * @brief Basis function for support point 5.
   */
  constexpr static double value5(const double xi) {
    return -19.7266861982493 * xi * xi * xi * xi * xi * xi * xi * xi -
           7.16310426200316 * xi * xi * xi * xi * xi * xi * xi + 44.7429986230572 * xi * xi * xi * xi * xi * xi +
           16.2469641839827 * xi * xi * xi * xi * xi - 32.3398660472438 * xi * xi * xi * xi -
           11.7431701395535 * xi * xi * xi + 7.32355362243583 * xi * xi + 2.65931021757392 * xi;
  }

  /**
   * @brief Basis function for support point 6.
   */
  constexpr static double value6(const double xi) {
    return 17.5609949844537 * xi * xi * xi * xi * xi * xi * xi * xi +
           11.8920648580289 * xi * xi * xi * xi * xi * xi * xi - 34.0932448057799 * xi * xi * xi * xi * xi * xi -
           23.0874776064749 * xi * xi * xi * xi * xi + 18.4067902908633 * xi * xi * xi * xi +
           12.4648258348041 * xi * xi * xi - 1.87454046953711 * xi * xi - 1.26941308635815 * xi;
  }

  /**
   * @brief Basis function for support point 7.
   */
  constexpr static double value7(const double xi) {
    return -13.6345320004901 * xi * xi * xi * xi * xi * xi * xi * xi -
           12.2677791811344 * xi * xi * xi * xi * xi * xi * xi + 21.6848443970084 * xi * xi * xi * xi * xi * xi +
           19.5111121254617 * xi * xi * xi * xi * xi - 8.87473674361952 * xi * xi * xi * xi -
           7.98511534224353 * xi * xi * xi + 0.824424347101285 * xi * xi + 0.741782397916255 * xi;
  }

  /**
   * @brief Basis function for support point 8.
   */
  constexpr static double value8(const double xi) {
    return 5.5859375 * xi * xi * xi * xi * xi * xi * xi * xi + 5.5859375 * xi * xi * xi * xi * xi * xi * xi -
           7.8203125 * xi * xi * xi * xi * xi * xi - 7.8203125 * xi * xi * xi * xi * xi +
           3.0078125 * xi * xi * xi * xi + 3.0078125 * xi * xi * xi - 0.2734375 * xi * xi - 0.2734375 * xi;
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
      case 8:
        return gradient8(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Gradient of basis function for support point 0.
   */
  constexpr static double gradient0(const double xi) {
    return 44.6875 * xi * xi * xi * xi * xi * xi * xi - 39.1015625 * xi * xi * xi * xi * xi * xi -
           46.921875 * xi * xi * xi * xi * xi + 39.1015625 * xi * xi * xi * xi + 12.03125 * xi * xi * xi -
           9.0234375 * xi * xi - 0.546875 * xi + 0.2734375;
  }

  /**
   * @brief Gradient of basis function for support point 1.
   */
  constexpr static double gradient1(const double xi) {
    return -109.076256003921 * xi * xi * xi * xi * xi * xi * xi + 85.8744542679408 * xi * xi * xi * xi * xi * xi +
           130.10906638205 * xi * xi * xi * xi * xi - 97.5555606273084 * xi * xi * xi * xi -
           35.4989469744781 * xi * xi * xi + 23.9553460267306 * xi * xi + 1.64884869420257 * xi - 0.741782397916255;
  }

  /**
   * @brief Gradient of basis function for support point 2.
   */
  constexpr static double gradient2(const double xi) {
    return 140.48795987563 * xi * xi * xi * xi * xi * xi * xi - 83.2444540062024 * xi * xi * xi * xi * xi * xi -
           204.559468834679 * xi * xi * xi * xi * xi + 115.437388032374 * xi * xi * xi * xi +
           73.6271611634532 * xi * xi * xi - 37.3944775044123 * xi * xi - 3.74908093907423 * xi + 1.26941308635815;
  }

  /**
   * @brief Gradient of basis function for support point 3.
   */
  constexpr static double gradient3(const double xi) {
    return -157.813489585994 * xi * xi * xi * xi * xi * xi * xi + 50.1417298340221 * xi * xi * xi * xi * xi * xi +
           268.457991738343 * xi * xi * xi * xi * xi - 81.2348209199137 * xi * xi * xi * xi -
           129.359464188975 * xi * xi * xi + 35.2295104186605 * xi * xi + 14.6471072448717 * xi - 2.65931021757392;
  }

  /**
   * @brief Gradient of basis function for support point 4 (centre node, odd
   * polynomial).
   */
  constexpr static double gradient4(const double xi) {
    return 163.428571428571 * xi * xi * xi * xi * xi * xi * xi - 294.171428571428 * xi * xi * xi * xi * xi +
           158.4 * xi * xi * xi - 24.0 * xi;
  }

  /**
   * @brief Gradient of basis function for support point 5.
   */
  constexpr static double gradient5(const double xi) {
    return -157.813489585994 * xi * xi * xi * xi * xi * xi * xi - 50.1417298340221 * xi * xi * xi * xi * xi * xi +
           268.457991738343 * xi * xi * xi * xi * xi + 81.2348209199137 * xi * xi * xi * xi -
           129.359464188975 * xi * xi * xi - 35.2295104186605 * xi * xi + 14.6471072448717 * xi + 2.65931021757392;
  }

  /**
   * @brief Gradient of basis function for support point 6.
   */
  constexpr static double gradient6(const double xi) {
    return 140.48795987563 * xi * xi * xi * xi * xi * xi * xi + 83.2444540062024 * xi * xi * xi * xi * xi * xi -
           204.559468834679 * xi * xi * xi * xi * xi - 115.437388032374 * xi * xi * xi * xi +
           73.6271611634532 * xi * xi * xi + 37.3944775044123 * xi * xi - 3.74908093907423 * xi - 1.26941308635815;
  }

  /**
   * @brief Gradient of basis function for support point 7.
   */
  constexpr static double gradient7(const double xi) {
    return -109.076256003921 * xi * xi * xi * xi * xi * xi * xi - 85.8744542679408 * xi * xi * xi * xi * xi * xi +
           130.10906638205 * xi * xi * xi * xi * xi + 97.5555606273084 * xi * xi * xi * xi -
           35.4989469744781 * xi * xi * xi - 23.9553460267306 * xi * xi + 1.64884869420257 * xi + 0.741782397916255;
  }

  /**
   * @brief Gradient of basis function for support point 8.
   */
  constexpr static double gradient8(const double xi) {
    return 44.6875 * xi * xi * xi * xi * xi * xi * xi + 39.1015625 * xi * xi * xi * xi * xi * xi -
           46.921875 * xi * xi * xi * xi * xi - 39.1015625 * xi * xi * xi * xi + 12.03125 * xi * xi * xi +
           9.0234375 * xi * xi - 0.546875 * xi - 0.2734375;
  }

  /**
   * @brief Gradient of basis function q evaluated at support point p.
   * Full 9x9 table computed via computGL.py.
   * Anti-symmetry: gradientAt(q,p) == -gradientAt(8-q, 8-p).
   * Corner formulas: gradientAt(0,0) = -N*(N-1)/4 = -18, gradientAt(8,8) = +18.
   */
  constexpr static double gradientAt(const int q, const int p) {
    switch (q) {
      case 0:
        switch (p) {
          case 0:
            return -18.0000000000000000;
          case 1:
            return -4.0870137020336728;
          case 2:
            return 0.9853600900745079;
          case 3:
            return -0.4446134492810899;
          case 4:
            return 0.2734375;
          case 5:
            return -0.2077345120355975;
          case 6:
            return 0.1896555919783545;
          case 7:
            return -0.2156540187025087;
          case 8:
            return 0.5;
        }
        break;
      case 1:
        switch (p) {
          case 0:
            return 24.3497451715930495;
          case 1:
            return 0.0;
          case 2:
            return -3.4883587534344560;
          case 3:
            return 1.2879607500639068;
          case 4:
            return -0.7417823979162546;
          case 5:
            return 0.5473001605340513;
          case 6:
            return -0.4923509383155125;
          case 7:
            return 0.5557049812837214;
          case 8:
            return -1.2848306326996024;
        }
        break;
      case 2:
        switch (p) {
          case 0:
            return -9.7387016572115428;
          case 1:
            return 5.7868058166373082;
          case 2:
            return 0.0;
          case 3:
            return -2.8344589120794206;
          case 4:
            return 1.2694130863581499;
          case 5:
            return -0.8557261850926734;
          case 6:
            return 0.7383492771903919;
          case 7:
            return -0.8167563817413388;
          case 8:
            return 1.8744408734470928;
        }
        break;
      case 3:
        switch (p) {
          case 0:
            return 5.5449639069493628;
          case 1:
            return -2.6960654403140580;
          case 2:
            return 3.5766809401256214;
          case 3:
            return 0.0;
          case 4:
            return -2.6593102175739181;
          case 5:
            return 1.3769648937605121;
          case 6:
            return -1.0798038112826394;
          case 7:
            return 1.1456537384550955;
          case 8:
            return -2.5907456765594361;
        }
        break;
      case 4:
        switch (p) {
          case 0:
            return -3.6571428571431284;
          case 1:
            return 1.6652216450052322;
          case 2:
            return -1.7178321571951187;
          case 3:
            return 2.8519159684628885;
          case 4:
            return 0.0;
          case 5:
            return -2.8519159684628899;
          case 6:
            return 1.7178321571951045;
          case 7:
            return -1.6652216450052890;
          case 8:
            return 3.6571428571430147;
        }
        break;
      case 5:
        switch (p) {
          case 0:
            return 2.5907456765594290;
          case 1:
            return -1.1456537384550778;
          case 2:
            return 1.0798038112826447;
          case 3:
            return -1.3769648937605121;
          case 4:
            return 2.6593102175739181;
          case 5:
            return 0.0;
          case 6:
            return -3.5766809401256232;
          case 7:
            return 2.6960654403141042;
          case 8:
            return -5.5449639069493131;
        }
        break;
      case 6:
        switch (p) {
          case 0:
            return -1.8744408734470359;
          case 1:
            return 0.8167563817413530;
          case 2:
            return -0.7383492771904061;
          case 3:
            return 0.8557261850926716;
          case 4:
            return -1.2694130863581499;
          case 5:
            return 2.8344589120794215;
          case 6:
            return 0.0;
          case 7:
            return -5.7868058166372940;
          case 8:
            return 9.7387016572115712;
        }
        break;
      case 7:
        switch (p) {
          case 0:
            return 1.2848306326995953;
          case 1:
            return -0.5557049812837036;
          case 2:
            return 0.4923509383155089;
          case 3:
            return -0.5473001605340508;
          case 4:
            return 0.7417823979162546;
          case 5:
            return -1.2879607500639065;
          case 6:
            return 3.4883587534344596;
          case 7:
            return 0.0;
          case 8:
            return -24.3497451715930424;
        }
        break;
      case 8:
        switch (p) {
          case 0:
            return -0.5;
          case 1:
            return 0.2156540187025087;
          case 2:
            return -0.1896555919783545;
          case 3:
            return 0.2077345120355977;
          case 4:
            return -0.2734375;
          case 5:
            return 0.4446134492810899;
          case 6:
            return -0.9853600900745070;
          case 7:
            return 4.0870137020336728;
          case 8:
            return 18.0000000000000000;
        }
        break;
    }
    return 0.0;
  }

  /* Tensor product helpers (2D / 3D) */
  struct TensorProduct2D {
    constexpr static int numSupportPoints1D = LagrangeBasis8GL::numSupportPoints;
    constexpr static int numSupportPoints = numSupportPoints1D * numSupportPoints1D;  // 9*9 = 81

    constexpr static int linearIndex(const int i, const int j) { return i + numSupportPoints1D * j; }

    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1) {
      i1 = linearIndex / numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    static void value(const double (&coords)[2], double (&N)[numSupportPoints]) {
      for (int a = 0; a < LagrangeBasis8GL::numSupportPoints; ++a) {
        for (int b = 0; b < LagrangeBasis8GL::numSupportPoints; ++b) {
          const int lindex = a + LagrangeBasis8GL::numSupportPoints * b;
          N[lindex] = LagrangeBasis8GL::value(a, coords[0]) * LagrangeBasis8GL::value(b, coords[1]);
        }
      }
    }
  };

  struct TensorProduct3D {
    constexpr static int numSupportPoints1D = LagrangeBasis8GL::numSupportPoints;
    constexpr static int numSupportPoints = numSupportPoints1D * numSupportPoints1D * numSupportPoints1D;  // 9^3 = 729

    constexpr static int linearIndex(const int i, const int j, const int k) {
      return i + numSupportPoints1D * j + numSupportPoints1D * numSupportPoints1D * k;
    }

    constexpr static void multiIndex(const int linearIndex, int& i0, int& i1, int& i2) {
      i2 = linearIndex / (numSupportPoints1D * numSupportPoints1D);
      i1 = (linearIndex % (numSupportPoints1D * numSupportPoints1D)) / numSupportPoints1D;
      i0 = linearIndex % numSupportPoints1D;
    }

    static void value(const double (&coords)[3], double (&N)[numSupportPoints]) {
      for (int a = 0; a < LagrangeBasis8GL::numSupportPoints; ++a) {
        for (int b = 0; b < LagrangeBasis8GL::numSupportPoints; ++b) {
          for (int c = 0; c < LagrangeBasis8GL::numSupportPoints; ++c) {
            const int lindex = a + LagrangeBasis8GL::numSupportPoints * b +
                               LagrangeBasis8GL::numSupportPoints * LagrangeBasis8GL::numSupportPoints * c;
            N[lindex] = LagrangeBasis8GL::value(a, coords[0]) * LagrangeBasis8GL::value(b, coords[1]) *
                        LagrangeBasis8GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif /* _LAGRANGEBASIS8GL_HPP_ */
