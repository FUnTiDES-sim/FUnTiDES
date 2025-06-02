//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.hpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#ifndef SEM_SOLVER_HPP_
#define SEM_SOLVER_HPP_

// #include "SEMQkGL.hpp"
#include <BasisFunctions.hpp>
#include <Integrals.hpp>
#include <model.hpp>
#ifdef USE_CALIPER
#include <caliper/cali.h>
#endif
#include <cmath>

class SEMsolver {
public:
  PROXY_HOST_DEVICE SEMsolver() {};
  PROXY_HOST_DEVICE ~SEMsolver(){};

  /**
   * @brief computeFEInit function:
   * init all FE components for computing mass and stiffness matrices
   */
  void computeFEInit(SEMinfo &myInfo, Mesh mesh);

  /**
   * @brief computeOneStep function:
   * init all FE components for computing mass and stiffness matrices
   */

  void computeOneStep(const int &timeSample, const int &order,
                      const int &nPointsPerElement, const int &i1,
                      const int &i2, SEMinfo &myInfo,
                      const arrayReal &myRHSTerm, arrayReal const &myPnGlobal,
                      const vectorInt &myRhsElement);

  void outputPnValues(Mesh mesh, const int &indexTimeStep, int &i1,
                      int &myElementSource, const arrayReal &pnGlobal);

  void initFEarrays(SEMinfo &myInfo, Mesh mesh);

  void allocateFEarrays(SEMinfo &myInfo);

  /**
   * @brief Compute coefficients for the taper layers. In this computation the
   * choice of the taper length and the coefficient of reflection (r)
   * highly depends on the model. Usually R will be between 10^{-3} and 1
   * and you need to find a compromise with sizeT.
   *
   * @param[in] vMin Min wavespeed (P-wavespeed for acoustic, S-wavespeed for
   * elastic)
   * @param[in] r desired reflectivity of the Taper
   */
  void getSpongeValues(Mesh &mesh, SEMinfo &myInfo, const float vMin,
                       const float r);

  void applyTaperCoeff(SEMinfo myinfo);

  // Define dummy structure for computing sponge value
  struct MinMax3D {
    float min_x, max_x;
    float min_y, max_y;
    float min_z, max_z;

#ifdef USE_KOKKOS
    KOKKOS_INLINE_FUNCTION
    MinMax3D()
        : min_x(Kokkos::reduction_identity<float>::max()),
          max_x(Kokkos::reduction_identity<float>::min()),
          min_y(Kokkos::reduction_identity<float>::max()),
          max_y(Kokkos::reduction_identity<float>::min()),
          min_z(Kokkos::reduction_identity<float>::max()),
          max_z(Kokkos::reduction_identity<float>::min()) {}
#else
    MinMax3D()
        : min_x(std::numeric_limits<double>::max()),
          max_x(std::numeric_limits<double>::lowest()),
          min_y(std::numeric_limits<double>::max()),
          max_y(std::numeric_limits<double>::lowest()),
          min_z(std::numeric_limits<double>::max()),
          max_z(std::numeric_limits<double>::lowest()) {}
#endif // USE_KOKKOS

    KOKKOS_INLINE_FUNCTION
    void operator+=(const MinMax3D &rhs) {
      min_x = min_x < rhs.min_x ? min_x : rhs.min_x;
      max_x = max_x > rhs.max_x ? max_x : rhs.max_x;

      min_y = min_y < rhs.min_y ? min_y : rhs.min_y;
      max_y = max_y > rhs.max_y ? max_y : rhs.max_y;

      min_z = min_z < rhs.min_z ? min_z : rhs.min_z;
      max_z = max_z > rhs.max_z ? max_z : rhs.max_z;
    }
  };

private:
  int order;
  // SEMQkGL myQk;
  SEMQkGLBasisFunctions myQkBasis;
  SEMQkGLIntegrals myQkIntegrals;

  // shared arrays
  arrayInt globalNodesList;
  arrayReal globalNodesCoordsX;
  arrayReal globalNodesCoordsY;
  arrayReal globalNodesCoordsZ;
  vectorInt listOfInteriorNodes;
  vectorInt listOfDampingNodes;
  vectorInt listOfSpongeNodes;
  arrayReal spongeTaperCoeff;

  // get model
  vectorReal model;
  double vMin; // min wavespeed in model

  // get quadrature points and weights
  vectorDouble quadraturePoints;
  vectorDouble weights;

  // get basis function and corresponding derivatives
  arrayDouble derivativeBasisFunction1D;

  // shared arrays
  vectorReal massMatrixGlobal;
  vectorReal yGlobal;
  vectorReal dampingValues;
  vectorReal dampingDistanceValues;
};
#endif // SEM_SOLVER_HPP_
