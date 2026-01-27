#ifndef SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_
#define SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_

#include "data_type.h"
#include "parallel_topology.h"
#include "sem_macros.h"

/**
 * @namespace model
 * @brief Namespace containing model and mesh representation classes
 */
namespace model
{

/**
 * @struct ModelDataBase
 * @brief Base structure for model data storage
 *
 * This structure serves as a base for storing model-specific data
 * that can be passed to ModelApi implementations.
 *
 * @tparam FloatType Floating-point type for physical quantities (e.g., float,
 * double)
 * @tparam ScalarType Integer type for indices (e.g., int, size_t)
 */
template <typename FloatType, typename ScalarType>
struct ModelDataBase
{
  PROXY_HOST_DEVICE ModelDataBase() = default;
  PROXY_HOST_DEVICE ~ModelDataBase() = default;
  PROXY_HOST_DEVICE ModelDataBase(const ModelDataBase&) = default;
  PROXY_HOST_DEVICE ModelDataBase& operator=(const ModelDataBase&) = default;
};

/**
 * @enum CubicFace
 * @brief Local face identifiers for cubic elements
 *
 * Convention: even indices = minus faces, odd indices = plus faces.
 * Face normals point outward from the element.
 */
enum class CubicFace : int
{
  kXMinus = 0,  ///< Face at x = x_min (left face, normal = [-1, 0, 0])
  kXPlus = 1,   ///< Face at x = x_max (right face, normal = [+1, 0, 0])
  kYMinus = 2,  ///< Face at y = y_min (front face, normal = [0, -1, 0])
  kYPlus = 3,   ///< Face at y = y_max (back face, normal = [0, +1, 0])
  kZMinus = 4,  ///< Face at z = z_min (bottom face, normal = [0, 0, -1])
  kZPlus = 5    ///< Face at z = z_max (top face, normal = [0, 0, +1])
};

/**
 * @class ModelApi
 * @brief Abstract base class representing a 3D mesh and material model.
 *
 * This class provides a unified interface for accessing mesh geometry,
 * material properties, and boundary conditions. Implementations can
 * represent structured or unstructured meshes with various material models
 * (acoustic, elastic, anisotropic, etc.).
 *
 * @tparam FloatType Floating-point type for physical quantities (e.g., float,
 * double)
 * @tparam ScalarType Integer type for indices (e.g., int, size_t)
 */
template <typename FloatType, typename ScalarType>
class ModelApi
{
 public:
  PROXY_HOST_DEVICE ModelApi() = default;

  /**
   * @brief Construct ModelApi from model data
   * @param data Model-specific data structure
   */
  PROXY_HOST_DEVICE ModelApi(const ModelDataBase<ScalarType, FloatType>& data)
  {
  }

  PROXY_HOST_DEVICE ModelApi(const ModelApi&) = default;
  PROXY_HOST_DEVICE ModelApi& operator=(const ModelApi&) = default;
  PROXY_HOST_DEVICE ~ModelApi() = default;

  /**
   * @brief Get the coordinate of a node along a specific dimension
   * @param dofGlobal Global node index
   * @param dim Dimension index (0=x, 1=y, 2=z)
   * @return Coordinate value
   */
  PROXY_HOST_DEVICE
  virtual FloatType nodeCoord(ScalarType dofGlobal, int dim) const = 0;

  /**
   * @brief Get the global node index for a given element and local indices
   * @param e Element index
   * @param i Local index in x-direction (0 to order)
   * @param j Local index in y-direction (0 to order)
   * @param k Local index in z-direction (0 to order)
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  virtual ScalarType globalNodeIndex(ScalarType e, int i, int j,
                                     int k) const = 0;

  /**
   * @brief Get P-wave velocity at a node
   * @param n Node index
   * @return P-wave velocity (Vp) in m/s
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVpOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get P-wave velocity for an element
   * @param e Element index
   * @return P-wave velocity (Vp) in m/s
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVpOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get density at a node
   * @param n Node index
   * @return Density (rho) in kg/m³
   */
  PROXY_HOST_DEVICE virtual FloatType getModelRhoOnNodes(
      ScalarType n) const = 0;

  /**
   * @brief Get density for an element
   * @param e Element index
   * @return Density (rho) in kg/m³
   */
  PROXY_HOST_DEVICE virtual FloatType getModelRhoOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get S-wave velocity at a node
   * @param n Node index
   * @return S-wave velocity (Vs) in m/s
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVsOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get S-wave velocity for an element
   * @param e Element index
   * @return S-wave velocity (Vs) in m/s
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVsOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get Thomsen delta parameter at a node (anisotropy)
   * @param n Node index
   * @return Delta parameter (dimensionless)
   */
  PROXY_HOST_DEVICE virtual FloatType getModelDeltaOnNodes(
      ScalarType n) const = 0;

  /**
   * @brief Get Thomsen delta parameter for an element (anisotropy)
   * @param e Element index
   * @return Delta parameter (dimensionless)
   */
  PROXY_HOST_DEVICE virtual FloatType getModelDeltaOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get Thomsen epsilon parameter at a node (anisotropy)
   * @param n Node index
   * @return Epsilon parameter (dimensionless)
   */
  PROXY_HOST_DEVICE virtual FloatType getModelEpsilonOnNodes(
      ScalarType n) const = 0;

  /**
   * @brief Get Thomsen epsilon parameter for an element (anisotropy)
   * @param e Element index
   * @return Epsilon parameter (dimensionless)
   */
  PROXY_HOST_DEVICE virtual FloatType getModelEpsilonOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get Thomsen gamma parameter at a node (anisotropy)
   * @param n Node index
   * @return Gamma parameter (dimensionless)
   */
  PROXY_HOST_DEVICE virtual FloatType getModelGammaOnNodes(
      ScalarType n) const = 0;

  /**
   * @brief Get Thomsen gamma parameter for an element (anisotropy)
   * @param e Element index
   * @return Gamma parameter (dimensionless)
   */
  PROXY_HOST_DEVICE virtual FloatType getModelGammaOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get theta angle at a node (anisotropy orientation)
   * @param n Node index
   * @return Theta angle in radians
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelThetaOnNodes(
      ScalarType n) const = 0;

  /**
   * @brief Get theta angle for an element (anisotropy orientation)
   * @param e Element index
   * @return Theta angle in radians
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelThetaOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get phi angle at a node (anisotropy orientation)
   * @param n Node index
   * @return Phi angle in radians
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelPhiOnNodes(
      ScalarType n) const = 0;

  /**
   * @brief Get phi angle for an element (anisotropy orientation)
   * @param e Element index
   * @return Phi angle in radians
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelPhiOnElement(
      ScalarType e) const = 0;

  /**
   * @brief Get the elasticity tensor for an element
   * @param e Element index
   * @param[out] CTTI 6x6 elasticity tensor in Voigt notation
   */
  PROXY_HOST_DEVICE
  virtual void getCTensorOnElement(ScalarType e,
                                   FloatType CTTI[6][6]) const = 0;

  /**
   * @brief Initialize and compute elasticity tensors for all elements
   *
   * This method should be called after model properties are set and
   * before running simulations requiring elasticity tensors.
   */
  virtual void initElasticityTensors() = 0;

  /**
   * @brief Get the total number of elements in the mesh
   * @return Number of elements
   */
  PROXY_HOST_DEVICE virtual ScalarType getNumberOfElements() const = 0;

  /**
   * @brief Get the total number of nodes in the mesh
   * @return Number of nodes
   */
  PROXY_HOST_DEVICE virtual ScalarType getNumberOfNodes() const = 0;

  /**
   * @brief Get the number of quadrature/collocation points per element
   * @return Number of points per element ((order+1)^3 for 3D hex elements)
   */
  PROXY_HOST_DEVICE virtual int getNumberOfPointsPerElement() const = 0;

  /**
   * @brief Get the polynomial order of the spectral element basis
   * @return Polynomial order (e.g., 4 for P4 elements)
   */
  PROXY_HOST_DEVICE virtual int getOrder() const = 0;

  /**
   * @brief Compute the outward unit normal vector of an element face.
   * @param e Element index
   * @param face Face identifier (using CubicFace enum)
   * @param[out] v Output array (size 3) holding the normal vector [nx, ny, nz]
   */
  PROXY_HOST_DEVICE
  virtual void faceNormal(ScalarType e, CubicFace face,
                          FloatType v[3]) const = 0;

  /**
   * @brief Get the domain size along a specific dimension
   * @param dim Dimension index (0=x, 1=y, 2=z)
   * @return Domain extent in the specified dimension
   */
  PROXY_HOST_DEVICE virtual FloatType domainSize(int dim) const = 0;

  /**
   * @brief Get the minimum grid spacing in the mesh
   * @return Minimum spacing between nodes in meters
   */
  PROXY_HOST_DEVICE virtual FloatType getMinSpacing() const = 0;

  /**
   * @brief Get the maximum wave speed in the model
   * @return Maximum speed (typically max Vp) in m/s
   */
  virtual FloatType getMaxSpeed() const = 0;

  /**
   * @brief Check if material properties are defined on nodes
   * @return true if properties are at nodes, false if at element centers
   */
  PROXY_HOST_DEVICE virtual bool isModelOnNodes() const = 0;

  /**
   * @brief Check if the model includes elastic (shear) properties
   * @return true for elastic models, false for acoustic models
   */
  PROXY_HOST_DEVICE virtual bool isElastic() const = 0;

  /**
   * @brief Build face connectivity for absorbing boundary conditions
   *
   * This method constructs data structures that identify which element
   * faces lie on absorbing boundaries and their connectivity information.
   */
  virtual void buildFaceConnectivity() = 0;

  /**
   * @brief Check if node is on free surface
   * @param n Node index
   * @return True if on free surface (top boundary)
   */
  PROXY_HOST_DEVICE
  virtual bool isFreeSurface(ScalarType n) const = 0;

  /**
   * @brief Enable/disable free surface on top boundary
   * @param enable True to enable free surface BC, false for absorbing
   * everywhere
   */
  virtual void setFreeSurfaceEnabled(bool enable) = 0;

  /**
   * @brief Initialize free surface detection
   *
   * Marks nodes on top boundary (Z+) as free surface.
   * Called once during mesh setup.
   */
  virtual void initFreeSurface() = 0;
};

}  // namespace model
#endif  // SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_