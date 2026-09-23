#ifndef FUNTIDES_SOLVER_FE_API_INCLUDE_SOLVER_H_
#define FUNTIDES_SOLVER_FE_API_INCLUDE_SOLVER_H_
#include <array>
#include <cmath>

#include "model.h"
#include "sem_enums.h"

namespace solver {
namespace fe {
/**
 * @brief Finite element wave propagation solver: assembles the operators of
 * one mesh (sub)domain once, then advances a wavefield in time.
 *
 * Typical use: optional set*() configuration, computeFEInit(), then per time
 * step either computeOneStep() or, in a distributed run, computeForces(),
 * synchronization of the force vectors, updateSolutionForward() (or
 * updateSolutionBackward()), followed by Wavefield::swap() by the caller.
 * Global vectors returned by reference have one value per global mesh node.
 *
 * @see docs/design.md, "Time levels and the split time step".
 */
class Solver {
 public:
  Solver() = default;
  virtual ~Solver() = default;

  /**
   * @brief Per-run state passed to the time-stepping methods (wavefield and
   * sources). Each solver requires the concrete data type that matches it.
   */
  struct DataStruct {
    PROXY_HOST_DEVICE
    virtual ~DataStruct() = default;

    /** @brief Print a debug summary to standard output. */
    virtual void print() const = 0;
  };

  /**
   * @brief Advance one time step in a non-distributed run: computeForces()
   * followed by the forward solution update.
   *
   * @param[in] dt Time step.
   * @todo VERIFY: is dt in seconds, and does the model impose a time unit?
   * @param[in] timeSample Time sample index into the source terms (second
   * index of Rhs::getTerm()).
   * @param[in,out] data Solver data; the new time level is written into its
   * wavefield.
   * @throws std::runtime_error if data is flagged as distributed.
   */
  virtual void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) = 0;

  /**
   * @brief Bind the solver to a mesh and build every time-independent
   * operator (mass, damping and sponge coefficients). Call once before time
   * stepping, after the set*() configuration methods.
   *
   * @param[in] mesh Mesh and model of the local subdomain; copied by the solver.
   * @param[in] sponge_size Thickness of the absorbing sponge layers along x, y
   * and z, in mesh coordinate units.
   * @param[in] surface_sponge When true, the free-surface side of the domain is
   * excluded from the sponge.
   * @todo VERIFY: which boundary is the free surface for surface_sponge (the
   * sponge code excludes the x = 0 side)?
   * @param[in] taper_delta_ Decay length of the sponge taper, in mesh coordinate
   * units.
   * @throws std::runtime_error if mesh is not of the mesh type the solver was
   * instantiated for.
   */
  virtual void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                             const bool surface_sponge, const float taper_delta_) = 0;

  /** @brief Fill the arrays allocated by allocateFEarrays(). */
  virtual void initFEarrays() = 0;

  /** @brief Allocate the global arrays sized by the bound mesh. */
  virtual void allocateFEarrays() = 0;

  /**
   * @brief Compute the sponge taper coefficients from the parameters given to
   * computeFEInit().
   */
  virtual void initSpongeValues() = 0;

  /**
   * @brief Zero the force vectors (see getForceVector()) before a new
   * accumulation.
   *
   * @param[in] numNodes Number of entries to zero, normally the number of
   * global mesh nodes.
   */
  virtual void resetGlobalVectors(int numNodes) = 0;

  /** @brief Assemble the lumped (diagonal) global mass matrix. */
  virtual void computeGlobalMassMatrix() = 0;

  /**
   * @brief Assemble the lumped damping matrices of the absorbing boundary
   * faces (see getDampingMatrix()).
   */
  virtual void computeDampingMatrix() = 0;

  /**
   * @brief Print a diagnostic value of a nodal field at one element.
   *
   * @param[in] t Time step index, for the message.
   * @param[in] e Element to probe.
   * @param[in] field Field with one value per global mesh node.
   * @param[in] fieldName Field name, for the message.
   */
  virtual void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) = 0;

  /**
   * @brief Print a diagnostic value of an element-wise field at one element.
   *
   * @param[in] t Time step index, for the message.
   * @param[in] e Element to probe.
   * @param[in] field Field of shape (number of elements, DOFs per element).
   * @param[in] fieldName Field name, for the message.
   */
  virtual void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) = 0;

  /**
   * @brief Number of solution components.
   * @return Bound of the component index of getForceVector() and
   * getDampingMatrix().
   */
  virtual int getNumComponents() const = 0;

  /**
   * @brief Lumped mass matrix of the acoustic (fluid) part, to be summed at
   * partition boundaries by a distributed driver.
   *
   * A solver with a single mass matrix returns it from both this method and
   * getMassMatrixElastic().
   *
   * @return One value per global mesh node.
   */
  virtual vectorReal& getMassMatrixAcoustic() = 0;

  /**
   * @brief Lumped mass matrix of the elastic (solid) part, to be summed at
   * partition boundaries by a distributed driver.
   *
   * A solver with a single mass matrix returns it from both this method and
   * getMassMatrixAcoustic().
   *
   * @return One value per global mesh node.
   */
  virtual vectorReal& getMassMatrixElastic() = 0;

  /**
   * @brief Acoustic/elastic interface coupling coefficient, the integral over
   * the interface of the basis function times the interface normal.
   *
   * Assembled from the locally owned acoustic element faces only, so a
   * distributed driver must sum it at partition boundaries.
   *
   * The argument selects the direction of the normal (0 = x, 1 = y, 2 = z).
   *
   * @return One value per global mesh node, or an empty view when the solver
   * has no such interface.
   */
  virtual vectorReal& getInterfaceCouplingCoeff(int) {
    static vectorReal empty;
    return empty;
  }

  /**
   * @brief Lumped damping matrix of one component, to be summed at partition
   * boundaries by a distributed driver after computeFEInit().
   *
   * @param[in] c Component index, in [0, getNumComponents()).
   * @return One value per global mesh node.
   */
  virtual vectorReal& getDampingMatrix(int c) = 0;

  /**
   * @brief Force vector of one component, filled by computeForces() and to be
   * summed at partition boundaries by a distributed driver at each time step.
   *
   * @param[in] component Component index, in [0, getNumComponents()).
   * @return One value per global mesh node.
   */
  virtual vectorReal& getForceVector(int component) = 0;

  /**
   * @brief First phase of a time step: fill the force vectors with the local
   * stiffness and source contributions.
   *
   * In a distributed run, the caller must sum the force vectors at partition
   * boundaries before calling updateSolutionForward() or
   * updateSolutionBackward().
   *
   * @param[in] dt Time step.
   * @param[in] timeSample Time sample index into the source terms.
   * @param[in] data Solver data providing the current fields and the sources.
   */
  virtual void computeForces(const float& dt, const int& timeSample, DataStruct& data) = 0;

  /**
   * @brief Second phase of a forward time step: compute the next time level
   * from the mass matrix and the assembled force vectors, and write it into
   * the previous buffer of the wavefield.
   *
   * @param[in] dt Time step.
   * @param[in,out] data Solver data whose wavefield has no previous-previous
   * buffer.
   * @throws std::runtime_error if the wavefield has a previous-previous buffer.
   */
  virtual void updateSolutionForward(const float& dt, DataStruct& data) = 0;

  /**
   * @brief Second phase of a backward (adjoint) time step: same as
   * updateSolutionForward(), but writes the next time level into the
   * previous-previous buffer.
   *
   * @param[in] dt Time step.
   * @param[in,out] data Solver data whose wavefield has a previous-previous
   * buffer.
   * @throws std::runtime_error if the wavefield has no previous-previous
   * buffer.
   */
  virtual void updateSolutionBackward(const float& dt, DataStruct& data) = 0;

  /**
   * @brief Select the anisotropy model of the medium.
   * @param[in] type Anisotropy model.
   */
  virtual void setAnisotropyType(model::AnisotropyType type) = 0;

  /**
   * @brief Set the coordinate of the plane z = constant that splits the mesh
   * into two element domains, for solvers that couple two discretizations.
   * No-op by default.
   * @todo VERIFY: which domain lies above the plane, and must this be called
   * before computeFEInit()?
   */
  virtual void setZBoundary(float) {}

  /**
   * @brief Give the split of the mesh into two element domains explicitly,
   * one tag per element of the mesh passed to computeFEInit(), instead of the
   * plane of setZBoundary(). No-op by default.
   *
   * Call before computeFEInit(). Needed on deformed meshes, where a z
   * threshold no longer separates the domains along a plane. The tag values
   * are defined by each solver.
   */
  virtual void setElementTags(const vectorInt&) {}

  /**
   * @brief Declare how the mesh builder filled the acoustic/elastic interface
   * nodes. Ignored by the solvers that have no such interface.
   */
  virtual void setInterfacePropertyConvention(utils::enums::interfacePropertyConvention) {}

  /**
   * @brief Enable viscoelastic attenuation with standard linear solids (SLS),
   * or disable it with an empty reference_frequencies.
   *
   * @param[in] reference_frequencies Relaxation angular frequency of each SLS,
   * in radians per unit of dt.
   * @param[in] anelasticity_coefficients Anelasticity coefficient of each SLS,
   * same size as reference_frequencies. When empty, the solver derives them
   * from the minimum quality factor of the model.
   * @throws std::runtime_error if the two vectors differ in size.
   * @todo VERIFY: must this be called before computeFEInit()?
   */
  virtual void setSLSAttenuation(const vectorReal& reference_frequencies,
                                 const vectorReal& anelasticity_coefficients = vectorReal()) = 0;
};
}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_API_INCLUDE_SOLVER_H_
