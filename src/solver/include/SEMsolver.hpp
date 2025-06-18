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
#include "dataType.hpp"
#include <BasisFunctions.hpp>
#include <Integrals.hpp>
#include <KokkosExp_InterOp.hpp>
#include <cmath>
#include <model.hpp>

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
   * @brief Compute one step of the spectral element wave equation solver.
   *
   * This function advances the pressure field `pnGlobal` by one time step using
   * a second-order explicit scheme. It resets global accumulators, applies RHS
   * forcing, computes local element contributions to mass and stiffness
   * matrices, updates pressure for interior nodes, and applies sponge damping.
   *
   * @param timeSample   Index of the current time step in `rhsTerm`
   * @param order        Spectral element interpolation order
   * @param nPointsPerElement Number of quadrature points per element
   * @param i1           Index for pressure at previous time step
   * @param i2           Index for pressure at current time step
   * @param myInfo       Structure containing mesh and solver configuration
   * @param rhsTerm      External forcing term, function of space and time
   * @param pnGlobal     2D array storing the global pressure field [node][time]
   * @param rhsElement   List of elements with a non-zero forcing term
   */
  void computeOneStep(const int &timeSample, const int &order,
                      const int &nPointsPerElement, const int &i1,
                      const int &i2, SEMinfo &myInfo, const arrayReal &rhsTerm,
                      const arrayReal &pnGlobal, const vectorInt &rhsElement);

  void computeOneStep_wrapper(int t, int order, int npts, int i1, int i2,
                              SEMinfo info,
                              Kokkos::Experimental::python_view_type_t<
                                  Kokkos::View<float **, Layout, MemSpace>>
                                  rhsTerm,
                              Kokkos::Experimental::python_view_type_t<
                                  Kokkos::View<float **, Layout, MemSpace>>
                                  pnGlobal,
                              Kokkos::Experimental::python_view_type_t<
                                  Kokkos::View<int *, Layout, MemSpace>>
                                  rhsElement);

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
  void initSpongeValues(Mesh &mesh, SEMinfo &myInfo);

  void spongeUpdate(const arrayReal &pnGlobal, const int i1, const int i2);

  /**
   * @brief Reset the global mass matrix and stiffness vector to zero.
   *
   * @param numNodes Total number of global nodes.
   */
  void resetGlobalVectors(int numNodes);

  /**
   * @brief Apply the external forcing term to the pressure field.
   *
   * @param timeSample Current time index into `rhsTerm`
   * @param i2 Index of the current time step in `pnGlobal`
   * @param rhsTerm Right-hand side values (forcing)
   * @param rhsElement Elements affected by the forcing term
   * @param myInfo Solver and mesh configuration
   * @param pnGlobal Pressure field array to update
   */
  void applyRHSTerm(int timeSample, int i2, const arrayReal &rhsTerm,
                    const vectorInt &rhsElement, SEMinfo &myInfo,
                    const arrayReal &pnGlobal);

  /**
   * @brief Compute local element contributions to the global mass and stiffness
   * system.
   *
   * @param order Polynomial interpolation order of the elements
   * @param nPointsPerElement Number of quadrature points per element
   * @param myInfo Solver configuration and mesh info
   * @param i2 Index of the current time step in `pnGlobal`
   * @param pnGlobal Global pressure field (used as input)
   */
  void computeElementContributions(int order, int nPointsPerElement,
                                   SEMinfo &myInfo, int i2,
                                   const arrayReal &pnGlobal);

  /**
   * @brief Update the pressure field for interior nodes using the time
   * integration scheme.
   *
   * @param i1 Index for pressure at the previous time step
   * @param i2 Index for pressure at the current time step
   * @param myInfo Solver and mesh configuration
   * @param pnGlobal Pressure field array (updated in-place)
   */
  void updatePressureField(int i1, int i2, SEMinfo &myInfo,
                           const arrayReal &pnGlobal);

  // Getters for shared arrays
  const arrayInt &getGlobalNodesList() const { return globalNodesList; }
  const arrayReal &getGlobalNodesCoordsX() const { return globalNodesCoordsX; }
  const arrayReal &getGlobalNodesCoordsY() const { return globalNodesCoordsY; }
  const arrayReal &getGlobalNodesCoordsZ() const { return globalNodesCoordsZ; }
  const vectorInt &getListOfInteriorNodes() const {
    return listOfInteriorNodes;
#ifdef USE_CALIPER
    mgr.flush();
#endif // USE_CALIPER
  }
  const vectorInt &getListOfDampingNodes() const { return listOfDampingNodes; }
  const vectorReal &getSpongeTaperCoeff() const { return spongeTaperCoeff; }

  void setGlobalNodesList(const arrayInt &list) { globalNodesList = list; }
  void setGlobalNodesCoordsX(const arrayReal &x) { globalNodesCoordsX = x; }
  void setGlobalNodesCoordsY(const arrayReal &y) { globalNodesCoordsY = y; }
  void setGlobalNodesCoordsZ(const arrayReal &z) { globalNodesCoordsZ = z; }
  void setListOfInteriorNodes(const vectorInt &nodes) {
    listOfInteriorNodes = nodes;
  }
  void setListOfDampingNodes(const vectorInt &nodes) {
    listOfDampingNodes = nodes;
  }
  void setSpongeTaperCoeff(const vectorReal &coeff) {
    spongeTaperCoeff = coeff;
  }

  void setOrder(int o) { order = o; }
  int getOrder() const { return order; }

  void setVMin(double v) { vMin = v; }
  double getVMin() const { return vMin; }

  void setModel(const vectorReal &m) { model = m; }
  void setModel_wrapper(Kokkos::Experimental::python_view_type_t<
                        Kokkos::View<float *, Layout, MemSpace>>
                            m) {
    setModel(m);
  }
  const vectorReal &getModel() const { return model; }

  void setMassMatrixGlobal(const vectorReal &m) { massMatrixGlobal = m; }
  const vectorReal &getMassMatrixGlobal() const { return massMatrixGlobal; }

private:
  int order;
  SEMinfo *myInfo;
  Mesh myMesh;
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
  // sponge boundaries data
  vectorReal spongeTaperCoeff;

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
