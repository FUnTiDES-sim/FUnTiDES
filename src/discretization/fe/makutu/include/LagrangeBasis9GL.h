#ifndef _LAGRANGEBASIS9GL_HPP_
#define _LAGRANGEBASIS9GL_HPP_

/**
 * @file LagrangeBasis9GL.hpp
 *
 * 1D Lagrange basis on Gauss-Lobatto-Legendre nodes for degree 9 (10 nodes).
 * Follows the same implementation pattern as LagrangeBasis8GL.hpp.
 * Polynomial formulas and gradientAt table generated via computGL.py.
 *
 * GLL nodes: { -1, -lambda1, -lambda2, -lambda3, -lambda4,
 *                  lambda4,  lambda3,  lambda2,  lambda1, 1 }
 *   lambda1 = 0.9195339081664589
 *   lambda2 = 0.7387738651055051
 *   lambda3 = 0.4779249498104445
 *   lambda4 = 0.1652789576663870
 */
class LagrangeBasis9GL
{
 public:
  constexpr static int numSupportPoints = 10;

  /// |nodes[8]| = |nodes[1]|
  static constexpr double lambda1 = 0.9195339081664589;

  /// |nodes[7]| = |nodes[2]|
  static constexpr double lambda2 = 0.7387738651055051;

  /// |nodes[6]| = |nodes[3]|
  static constexpr double lambda3 = 0.4779249498104445;

  /// |nodes[5]| = |nodes[4]|
  static constexpr double lambda4 = 0.1652789576663870;

  /**
   * @brief GLL quadrature weight for node q.
   */
  constexpr static double weight(const int q)
  {
    switch (q)
    {
      case 1:
      case 8:
        return 0.1333059908510701;
      case 2:
      case 7:
        return 0.2248893420631264;
      case 3:
      case 6:
        return 0.2920426836796838;
      case 4:
      case 5:
        return 0.3275397611838974;
      default:
        return 1.0 / 45.0;
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
        return -lambda3;
      case 4:
        return -lambda4;
      case 5:
        return lambda4;
      case 6:
        return lambda3;
      case 7:
        return lambda2;
      case 8:
        return lambda1;
      case 9:
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
      case 7:
        return value7(xi);
      case 8:
        return value8(xi);
      case 9:
        return value9(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Basis function for support point 0.
   */
  constexpr static double value0(const double xi)
  {
    return -9.49609375 * xi * xi * xi * xi * xi * xi * xi * xi * xi +
           9.49609375 * xi * xi * xi * xi * xi * xi * xi * xi +
           15.640625 * xi * xi * xi * xi * xi * xi * xi -
           15.640625 * xi * xi * xi * xi * xi * xi -
           7.8203125 * xi * xi * xi * xi * xi + 7.8203125 * xi * xi * xi * xi +
           1.203125 * xi * xi * xi - 1.203125 * xi * xi - 0.02734375 * xi +
           0.02734375;
  }

  /**
   * @brief Basis function for support point 1.
   */
  constexpr static double value1(const double xi)
  {
    return 23.2581991069276 * xi * xi * xi * xi * xi * xi * xi * xi * xi -
           21.3867027217068 * xi * xi * xi * xi * xi * xi * xi * xi -
           41.9000228289113 * xi * xi * xi * xi * xi * xi * xi +
           38.5284917441326 * xi * xi * xi * xi * xi * xi +
           22.033178498462 * xi * xi * xi * xi * xi -
           20.26025473402 * xi * xi * xi * xi -
           3.47055997155687 * xi * xi * xi + 3.19129757417176 * xi * xi +
           0.0792051950785524 * xi - 0.0728318625776681;
  }

  /**
   * @brief Basis function for support point 2.
   */
  constexpr static double value2(const double xi)
  {
    return -30.2089539641742 * xi * xi * xi * xi * xi * xi * xi * xi * xi +
           22.3175856809073 * xi * xi * xi * xi * xi * xi * xi * xi +
           63.4772291071541 * xi * xi * xi * xi * xi * xi * xi -
           46.8953178936799 * xi * xi * xi * xi * xi * xi -
           39.9888510087652 * xi * xi * xi * xi * xi +
           29.5427180208736 * xi * xi * xi * xi +
           6.87995289293189 * xi * xi * xi - 5.0827293904551 * xi * xi -
           0.159377027146516 * xi + 0.117743582354057;
  }

  /**
   * @brief Basis function for support point 3.
   */
  constexpr static double value3(const double xi)
  {
    return 34.4250370034963 * xi * xi * xi * xi * xi * xi * xi * xi * xi -
           16.4525840821186 * xi * xi * xi * xi * xi * xi * xi * xi -
           83.2619975287327 * xi * xi * xi * xi * xi * xi * xi +
           39.7929859900369 * xi * xi * xi * xi * xi * xi +
           66.0320305883065 * xi * xi * xi * xi * xi -
           31.5583549047981 * xi * xi * xi * xi -
           17.629048439256 * xi * xi * xi + 8.42536209053729 * xi * xi +
           0.433978376185841 * xi - 0.207409093657436;
  }

  /**
   * @brief Basis function for support point 4.
   */
  constexpr static double value4(const double xi)
  {
    return -36.457196112531 * xi * xi * xi * xi * xi * xi * xi * xi * xi +
           6.02560737291818 * xi * xi * xi * xi * xi * xi * xi * xi +
           95.5084365449144 * xi * xi * xi * xi * xi * xi * xi -
           15.7855348404897 * xi * xi * xi * xi * xi * xi -
           87.4617030627869 * xi * xi * xi * xi * xi +
           14.4555791179445 * xi * xi * xi * xi +
           32.2533814922412 * xi * xi * xi - 5.33080527425396 * xi * xi -
           3.84291886183778 * xi + 0.635153623881047;
  }

  /**
   * @brief Basis function for support point 5.
   */
  constexpr static double value5(const double xi)
  {
    return 36.457196112531 * xi * xi * xi * xi * xi * xi * xi * xi * xi +
           6.02560737291818 * xi * xi * xi * xi * xi * xi * xi * xi -
           95.5084365449144 * xi * xi * xi * xi * xi * xi * xi -
           15.7855348404897 * xi * xi * xi * xi * xi * xi +
           87.4617030627869 * xi * xi * xi * xi * xi +
           14.4555791179445 * xi * xi * xi * xi -
           32.2533814922412 * xi * xi * xi - 5.33080527425395 * xi * xi +
           3.84291886183778 * xi + 0.635153623881047;
  }

  /**
   * @brief Basis function for support point 6.
   */
  constexpr static double value6(const double xi)
  {
    return -34.4250370034963 * xi * xi * xi * xi * xi * xi * xi * xi * xi -
           16.4525840821186 * xi * xi * xi * xi * xi * xi * xi * xi +
           83.2619975287328 * xi * xi * xi * xi * xi * xi * xi +
           39.7929859900369 * xi * xi * xi * xi * xi * xi -
           66.0320305883065 * xi * xi * xi * xi * xi -
           31.5583549047981 * xi * xi * xi * xi +
           17.629048439256 * xi * xi * xi + 8.4253620905373 * xi * xi -
           0.433978376185841 * xi - 0.207409093657436;
  }

  /**
   * @brief Basis function for support point 7.
   */
  constexpr static double value7(const double xi)
  {
    return 30.2089539641742 * xi * xi * xi * xi * xi * xi * xi * xi * xi +
           22.3175856809073 * xi * xi * xi * xi * xi * xi * xi * xi -
           63.4772291071541 * xi * xi * xi * xi * xi * xi * xi -
           46.8953178936799 * xi * xi * xi * xi * xi * xi +
           39.9888510087652 * xi * xi * xi * xi * xi +
           29.5427180208736 * xi * xi * xi * xi -
           6.87995289293189 * xi * xi * xi - 5.0827293904551 * xi * xi +
           0.159377027146516 * xi + 0.117743582354057;
  }

  /**
   * @brief Basis function for support point 8.
   */
  constexpr static double value8(const double xi)
  {
    return -23.2581991069276 * xi * xi * xi * xi * xi * xi * xi * xi * xi -
           21.3867027217068 * xi * xi * xi * xi * xi * xi * xi * xi +
           41.9000228289113 * xi * xi * xi * xi * xi * xi * xi +
           38.5284917441326 * xi * xi * xi * xi * xi * xi -
           22.033178498462 * xi * xi * xi * xi * xi -
           20.26025473402 * xi * xi * xi * xi +
           3.47055997155687 * xi * xi * xi + 3.19129757417176 * xi * xi -
           0.0792051950785524 * xi - 0.0728318625776681;
  }

  /**
   * @brief Basis function for support point 9.
   */
  constexpr static double value9(const double xi)
  {
    return 9.49609375 * xi * xi * xi * xi * xi * xi * xi * xi * xi +
           9.49609375 * xi * xi * xi * xi * xi * xi * xi * xi -
           15.640625 * xi * xi * xi * xi * xi * xi * xi -
           15.640625 * xi * xi * xi * xi * xi * xi +
           7.8203125 * xi * xi * xi * xi * xi + 7.8203125 * xi * xi * xi * xi -
           1.203125 * xi * xi * xi - 1.203125 * xi * xi + 0.02734375 * xi +
           0.02734375;
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
      case 7:
        return gradient7(xi);
      case 8:
        return gradient8(xi);
      case 9:
        return gradient9(xi);
      default:
        return 0.0;
    }
  }

  /**
   * @brief Gradient of basis function for support point 0.
   */
  constexpr static double gradient0(const double xi)
  {
    return -85.46484375 * xi * xi * xi * xi * xi * xi * xi * xi +
           75.96875 * xi * xi * xi * xi * xi * xi * xi +
           109.484375 * xi * xi * xi * xi * xi * xi -
           93.84375 * xi * xi * xi * xi * xi - 39.1015625 * xi * xi * xi * xi +
           31.28125 * xi * xi * xi + 3.609375 * xi * xi - 2.40625 * xi -
           0.02734375;
  }

  /**
   * @brief Gradient of basis function for support point 1.
   */
  constexpr static double gradient1(const double xi)
  {
    return 209.323791962348 * xi * xi * xi * xi * xi * xi * xi * xi -
           171.093621773654 * xi * xi * xi * xi * xi * xi * xi -
           293.300159802379 * xi * xi * xi * xi * xi * xi +
           231.170950464796 * xi * xi * xi * xi * xi +
           110.16589249231 * xi * xi * xi * xi -
           81.0410189360798 * xi * xi * xi - 10.4116799146706 * xi * xi +
           6.38259514834352 * xi + 0.0792051950785524;
  }

  /**
   * @brief Gradient of basis function for support point 2.
   */
  constexpr static double gradient2(const double xi)
  {
    return -271.880585677568 * xi * xi * xi * xi * xi * xi * xi * xi +
           178.540685447258 * xi * xi * xi * xi * xi * xi * xi +
           444.340603750078 * xi * xi * xi * xi * xi * xi -
           281.371907362079 * xi * xi * xi * xi * xi -
           199.944255043826 * xi * xi * xi * xi +
           118.170872083495 * xi * xi * xi + 20.6398586787957 * xi * xi -
           10.1654587809102 * xi - 0.159377027146516;
  }

  /**
   * @brief Gradient of basis function for support point 3.
   */
  constexpr static double gradient3(const double xi)
  {
    return 309.825333031467 * xi * xi * xi * xi * xi * xi * xi * xi -
           131.620672656949 * xi * xi * xi * xi * xi * xi * xi -
           582.833982701129 * xi * xi * xi * xi * xi * xi +
           238.757915940222 * xi * xi * xi * xi * xi +
           330.160152941533 * xi * xi * xi * xi -
           126.233419619193 * xi * xi * xi - 52.8871453177679 * xi * xi +
           16.8507241810746 * xi + 0.433978376185841;
  }

  /**
   * @brief Gradient of basis function for support point 4.
   */
  constexpr static double gradient4(const double xi)
  {
    return -328.114765012779 * xi * xi * xi * xi * xi * xi * xi * xi +
           48.2048589833454 * xi * xi * xi * xi * xi * xi * xi +
           668.559055814401 * xi * xi * xi * xi * xi * xi -
           94.7132090429383 * xi * xi * xi * xi * xi -
           437.308515313934 * xi * xi * xi * xi +
           57.8223164717778 * xi * xi * xi + 96.7601444767235 * xi * xi -
           10.6616105485079 * xi - 3.84291886183778;
  }

  /**
   * @brief Gradient of basis function for support point 5.
   */
  constexpr static double gradient5(const double xi)
  {
    return 328.114765012779 * xi * xi * xi * xi * xi * xi * xi * xi +
           48.2048589833454 * xi * xi * xi * xi * xi * xi * xi -
           668.559055814401 * xi * xi * xi * xi * xi * xi -
           94.7132090429383 * xi * xi * xi * xi * xi +
           437.308515313934 * xi * xi * xi * xi +
           57.8223164717778 * xi * xi * xi - 96.7601444767235 * xi * xi -
           10.6616105485079 * xi + 3.84291886183778;
  }

  /**
   * @brief Gradient of basis function for support point 6.
   */
  constexpr static double gradient6(const double xi)
  {
    return -309.825333031467 * xi * xi * xi * xi * xi * xi * xi * xi -
           131.620672656949 * xi * xi * xi * xi * xi * xi * xi +
           582.833982701129 * xi * xi * xi * xi * xi * xi +
           238.757915940221 * xi * xi * xi * xi * xi -
           330.160152941533 * xi * xi * xi * xi -
           126.233419619193 * xi * xi * xi + 52.8871453177679 * xi * xi +
           16.8507241810746 * xi - 0.433978376185841;
  }

  /**
   * @brief Gradient of basis function for support point 7.
   */
  constexpr static double gradient7(const double xi)
  {
    return 271.880585677568 * xi * xi * xi * xi * xi * xi * xi * xi +
           178.540685447258 * xi * xi * xi * xi * xi * xi * xi -
           444.340603750078 * xi * xi * xi * xi * xi * xi -
           281.371907362079 * xi * xi * xi * xi * xi +
           199.944255043826 * xi * xi * xi * xi +
           118.170872083495 * xi * xi * xi - 20.6398586787957 * xi * xi -
           10.1654587809102 * xi + 0.159377027146516;
  }

  /**
   * @brief Gradient of basis function for support point 8.
   */
  constexpr static double gradient8(const double xi)
  {
    return -209.323791962348 * xi * xi * xi * xi * xi * xi * xi * xi -
           171.093621773654 * xi * xi * xi * xi * xi * xi * xi +
           293.300159802379 * xi * xi * xi * xi * xi * xi +
           231.170950464796 * xi * xi * xi * xi * xi -
           110.16589249231 * xi * xi * xi * xi -
           81.0410189360798 * xi * xi * xi + 10.4116799146706 * xi * xi +
           6.38259514834352 * xi - 0.0792051950785524;
  }

  /**
   * @brief Gradient of basis function for support point 9.
   */
  constexpr static double gradient9(const double xi)
  {
    return 85.46484375 * xi * xi * xi * xi * xi * xi * xi * xi +
           75.96875 * xi * xi * xi * xi * xi * xi * xi -
           109.484375 * xi * xi * xi * xi * xi * xi -
           93.84375 * xi * xi * xi * xi * xi + 39.1015625 * xi * xi * xi * xi +
           31.28125 * xi * xi * xi - 3.609375 * xi * xi - 2.40625 * xi +
           0.02734375;
  }

  /**
   * @brief Gradient of basis function q evaluated at support point p.
   * Full 10x10 table computed via computGL.py.
   * Anti-symmetry: gradientAt(q,p) == -gradientAt(9-q, 9-p).
   * Corner formulas: gradientAt(0,0) = -N*(N-1)/4 = -22.5, gradientAt(9,9) =
   * +22.5.
   */
  constexpr static double gradientAt(const int q, const int p)
  {
    switch (q)
    {
      case 0:
        switch (p)
        {
          case 0:
            return -22.5;
          case 1:
            return -5.0740647029780490;
          case 2:
            return 1.2033519928522214;
          case 3:
            return -0.5283693768202707;
          case 4:
            return 0.3120472556084114;
          case 5:
            return -0.2235279447424541;
          case 6:
            return 0.1866457893937360;
          case 7:
            return -0.1807865854892450;
          case 8:
            return 0.2127027580091791;
          case 9:
            return -0.5;
        }
        break;
      case 1:
        switch (p)
        {
          case 0:
            return 30.4381450292821043;
          case 1:
            return 0.0;
          case 2:
            return -4.2592973549651916;
          case 3:
            return 1.5299026381816052;
          case 4:
            return -0.8458135734064254;
          case 5:
            return 0.5880821430451697;
          case 6:
            return -0.4834623263339513;
          case 7:
            return 0.4642749589081347;
          case 8:
            return -0.5437537382357931;
          case 9:
            return 1.2759548360925805;
        }
        break;
      case 2:
        switch (p)
        {
          case 0:
            return -12.1779467074297827;
          case 1:
            return 7.1855028697058856;
          case 2:
            return 0.0;
          case 3:
            return -3.3641258682978226;
          case 4:
            return 1.4448503156016588;
          case 5:
            return -0.9165551803364349;
          case 6:
            return 0.7212373127215983;
          case 7:
            return -0.6767970871960785;
          case 8:
            return 0.7832392931379388;
          case 9:
            return -1.8295639319032944;
        }
        break;
      case 3:
        switch (p)
        {
          case 0:
            return 6.9437884851340499;
          case 1:
            return -3.3516638627466904;
          case 2:
            return 4.3686745570102090;
          case 3:
            return 0.0;
          case 4:
            return -3.0202179581993462;
          case 5:
            return 1.4680555093899925;
          case 6:
            return -1.0461893655025047;
          case 7:
            return 0.9366032131393922;
          case 8:
            return -1.0591544636454699;
          case 9:
            return 2.4528841754426622;
        }
        break;
      case 4:
        switch (p)
        {
          case 0:
            return -4.5993547611031005;
          case 1:
            return 2.0782079940363687;
          case 2:
            return -2.1043501794131387;
          case 3:
            return 3.3873181012024447;
          case 4:
            return 0.0;
          case 5:
            return -3.0251884877519748;
          case 6:
            return 1.6464940839870605;
          case 7:
            return -1.3349154838782162;
          case 8:
            return 1.4449484487513828;
          case 9:
            return -3.2946430337490966;
        }
        break;
      case 5:
        switch (p)
        {
          case 0:
            return 3.2946430337492529;
          case 1:
            return -1.4449484487515178;
          case 2:
            return 1.3349154838782340;
          case 3:
            return -1.6464940839870614;
          case 4:
            return 3.0251884877519739;
          case 5:
            return 0.0;
          case 6:
            return -3.3873181012024567;
          case 7:
            return 2.1043501794131316;
          case 8:
            return -2.0782079940365179;
          case 9:
            return 4.5993547611030863;
        }
        break;
      case 6:
        switch (p)
        {
          case 0:
            return -2.4528841754425912;
          case 1:
            return 1.0591544636455126;
          case 2:
            return -0.9366032131393744;
          case 3:
            return 1.0461893655024967;
          case 4:
            return -1.4680555093899941;
          case 5:
            return 3.0202179581993471;
          case 6:
            return 0.0;
          case 7:
            return -4.3686745570101913;
          case 8:
            return 3.3516638627468467;
          case 9:
            return -6.9437884851339504;
        }
        break;
      case 7:
        switch (p)
        {
          case 0:
            return 1.8295639319031682;
          case 1:
            return -0.7832392931379495;
          case 2:
            return 0.6767970871960758;
          case 3:
            return -0.7212373127216001;
          case 4:
            return 0.9165551803364348;
          case 5:
            return -1.4448503156016588;
          case 6:
            return 3.3641258682978226;
          case 7:
            return 0.0;
          case 8:
            return -7.1855028697058749;
          case 9:
            return 12.1779467074297809;
        }
        break;
      case 8:
        switch (p)
        {
          case 0:
            return -1.2759548360924953;
          case 1:
            return 0.5437537382357931;
          case 2:
            return -0.4642749589081276;
          case 3:
            return 0.4834623263339504;
          case 4:
            return -0.5880821430451697;
          case 5:
            return 0.8458135734064252;
          case 6:
            return -1.5299026381816043;
          case 7:
            return 4.2592973549651987;
          case 8:
            return 0.0;
          case 9:
            return -30.4381450292821114;
        }
        break;
      case 9:
        switch (p)
        {
          case 0:
            return 0.5;
          case 1:
            return -0.2127027580091769;
          case 2:
            return 0.1807865854892459;
          case 3:
            return -0.1866457893937387;
          case 4:
            return 0.2235279447424537;
          case 5:
            return -0.3120472556084113;
          case 6:
            return 0.5283693768202747;
          case 7:
            return -1.2033519928522045;
          case 8:
            return 5.0740647029781041;
          case 9:
            return 22.5;
        }
        break;
    }
    return 0.0;
  }

  /* Tensor product helpers (2D / 3D) */
  struct TensorProduct2D
  {
    constexpr static int numSupportPoints1D =
        LagrangeBasis9GL::numSupportPoints;
    constexpr static int numSupportPoints =
        numSupportPoints1D * numSupportPoints1D;  // 10*10 = 100

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
      for (int a = 0; a < LagrangeBasis9GL::numSupportPoints; ++a)
      {
        for (int b = 0; b < LagrangeBasis9GL::numSupportPoints; ++b)
        {
          const int lindex = a + LagrangeBasis9GL::numSupportPoints * b;
          N[lindex] = LagrangeBasis9GL::value(a, coords[0]) *
                      LagrangeBasis9GL::value(b, coords[1]);
        }
      }
    }
  };

  struct TensorProduct3D
  {
    constexpr static int numSupportPoints1D =
        LagrangeBasis9GL::numSupportPoints;
    constexpr static int numSupportPoints = numSupportPoints1D *
                                            numSupportPoints1D *
                                            numSupportPoints1D;  // 10^3 = 1000

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
      for (int a = 0; a < LagrangeBasis9GL::numSupportPoints; ++a)
      {
        for (int b = 0; b < LagrangeBasis9GL::numSupportPoints; ++b)
        {
          for (int c = 0; c < LagrangeBasis9GL::numSupportPoints; ++c)
          {
            const int lindex = a + LagrangeBasis9GL::numSupportPoints * b +
                               LagrangeBasis9GL::numSupportPoints *
                                   LagrangeBasis9GL::numSupportPoints * c;
            N[lindex] = LagrangeBasis9GL::value(a, coords[0]) *
                        LagrangeBasis9GL::value(b, coords[1]) *
                        LagrangeBasis9GL::value(c, coords[2]);
          }
        }
      }
    }
  };
};

#endif /* _LAGRANGEBASIS9GL_HPP_ */
