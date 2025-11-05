#ifndef SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_
#define SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_

#include <data_type.h>
#include <sem_macros.h>

namespace model
{

template <typename FloatType, typename ScalarType>
struct ModelDataBase
{
  // GPU-compatible special member functions
  PROXY_HOST_DEVICE ModelDataBase() = default;
  PROXY_HOST_DEVICE ~ModelDataBase() = default;
  PROXY_HOST_DEVICE ModelDataBase(const ModelDataBase&) = default;
  PROXY_HOST_DEVICE ModelDataBase& operator=(const ModelDataBase&) = default;
};

/**
 * @enum BoundaryFlag
 * @brief Flags representing the boundary condition type of a mesh node.
 */
enum BoundaryFlag : uint8_t
{
  InteriorNode = 0,  ///< Node inside the domain
  Damping = 1 << 0,  ///< Node in damping boundary zone
  Sponge = 1 << 1,   ///< Node in sponge layer
  Surface = 1 << 2,  ///< Node on a free surface
  Ghost = 1 << 3     ///< Ghost node for halo/exchange
};

/**
 * @brief Abstract base class representing a structured 3D mesh.
 */
template <typename FloatType, typename ScalarType>
class ModelApi
{
 public:
  /**
   * @brief Default constructor.
   */
  PROXY_HOST_DEVICE ModelApi() = default;

  /**
   * @brief Constructor from ModelDataBase.
   * @param data ModelDataBase structure containing all the mesh data
   */
  PROXY_HOST_DEVICE ModelApi(const ModelDataBase<ScalarType, FloatType>& data)
  {
  }

  /**
   * @brief Copy constructor.
   */
  PROXY_HOST_DEVICE ModelApi(const ModelApi&) = default;

  /**
   * @brief Assignment operator.
   */
  PROXY_HOST_DEVICE ModelApi& operator=(const ModelApi&) = default;

  /**
   * @brief Destructor.
   */
  PROXY_HOST_DEVICE ~ModelApi() = default;

  /**
   * @brief Get the coordinate of a global node in the given dimension.
   * @param dofGlobal Global node index
   * @param dim Dimension index (0 = x, 1 = y, 2 = z)
   * @return Coordinate value in the specified dimension
   */
  PROXY_HOST_DEVICE
  virtual FloatType nodeCoord(ScalarType dofGlobal, int dim) const = 0;

  /**
   * @brief Get the global node index for a local element-node triplet.
   * @param e Element index
   * @param i Local i-index in the element
   * @param j Local j-index in the element
   * @param k Local k-index in the element
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  virtual ScalarType globalNodeIndex(ScalarType e, int i, int j,
                                     int k) const = 0;

  /**
   * @brief Get the P-wave velocity value at a global node.
   * @param n Global node index
   * @return Model P-wave velocity value at the node
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelVpOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average P-wave velocity value on a given element.
   * @param e Element index
   * @return Model P-wave velocity value for the element
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelVpOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the density value at a global node.
   * @param n Global node index
   * @return Model density value at the node
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelRhoOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average density value on a given element.
   * @param e Element index
   * @return Model density value for the element
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelRhoOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the S-wave velocity value at a global node.
   * @param n Global node index
   * @return Model P-wave velocity value at the node
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelVsOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average S-wave velocity value on a given element.
   * @param e Element index
   * @return Model P-wave velocity value for the element
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelVsOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the average Thomsen parameter delta value on a given nodes.
   * @param n Global node index
   * @return Model Thomsen paramter delta value for the node
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelDeltaOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average Thomsen parameter delta value on a given element.
   * @param e Element index
   * @return Model Thomsen paramter delta value for the element
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelDeltaOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the average Thomsen parameter epsilon value on a given node.
   * @param n Global node index
   * @return Model Thomsen paramter epsilon value for the node
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelEpsilonOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average Thomsen parameter epsilon value on a given element.
   * @param e Element index
   * @return Model Thomsen paramter epsilon value for the element
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelEpsilonOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the average Thomsen parameter gamma value on a given node.
   * @param n Global node index
   * @return Model Thomsen paramter gamma value for the node
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelGammaOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average Thomsen parameter gamma value on a given element.
   * @param e Element index
   * @return Model Thomsen paramter gamma value for the element
   */
  PROXY_HOST_DEVICE
  virtual FloatType getModelGammaOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the average anisotropic parameter theta value on a given node.
   * @param n Global node index
   * @return Model anisotropic paramter theta value for the node
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getModelThetaOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average anisotropic parameter theta value on a given
   * element.
   * @param e Element index
   * @return Model anisotropic paramter theta value for the element
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getModelThetaOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the average anisotropic parameter phi value on a given node.
   * @param n Global node index
   * @return Model anisotropic paramter phi value for the node
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getModelPhiOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the average anisotropic parameter phi value on a given element.
   * @param e Element index
   * @return Model anisotropic paramter phi value for the element
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getModelPhiOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the total number of elements in the mesh.
   * @return Total element count
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getNumberOfElements() const = 0;

  /**
   * @brief Get the total number of global nodes in the mesh.
   * @return Total node count
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getNumberOfNodes() const = 0;

  /**
   * @brief Get the number of interpolation points per element.
   * @return Number of interpolation points in one element
   */
  PROXY_HOST_DEVICE
  virtual int getNumberOfPointsPerElement() const = 0;

  /**
   * @brief Get the polynomial order of the elements.
   * @return ORDER
   */
  PROXY_HOST_DEVICE
  virtual int getOrder() const = 0;

  /**
   * @brief Get the boundary type of a given node.
   * @param n Global node index
   * @return A combination of BoundaryFlag values
   */
  PROXY_HOST_DEVICE
  virtual BoundaryFlag boundaryType(ScalarType n) const = 0;

  /**
   * @brief Compute the outward unit normal vector of an element face.
   * @param e Element index
   * @param dir Axis direction (0 = x, 1 = y, 2 = z)
   * @param face 1 = negative side, 2 = positive side in that direction
   * @param[out] v Output array (size 3) holding the normal vector
   */
  PROXY_HOST_DEVICE
  virtual void faceNormal(ScalarType e, int dir, int face,
                          FloatType v[3]) const = 0;

  /**
   * @brief Get the size of the domain in the specified dimension.
   * @param dim The dimension index (0 for X, 1 for Y, 2 for Z).
   * @return The size of the domain along the specified dimension.
   */
  PROXY_HOST_DEVICE
  virtual FloatType domainSize(int dim) const = 0;

  /**
   * @brief Gets the minimum spacing between grid points or nodes
   * Returns the smallest distance between adjacent computational points
   * in the discretized domain.
   * @return Minimum spacing value between grid points
   */
  PROXY_HOST_DEVICE
  virtual FloatType getMinSpacing() const = 0;

  /**
   * @brief Gets the maximum wave speed in the computational domain
   * Returns the highest wave velocity present in the velocity model.
   * @return Maximum wave speed value
   */
  virtual FloatType getMaxSpeed() const = 0;

  /**
   * @brief Indicates if the model properties are defined at nodes
   * Returns true if model parameters (e.g., velocity, density) are
   * specified at mesh nodes, false if defined at element centers.
   * @return True if model is defined on nodes, false if on elements
   */
  PROXY_HOST_DEVICE
  virtual bool isModelOnNodes() const = 0;

  /**
   * @brief Indicates if the model is elastic.
   * Returns true if the model supports elastic wave propagation,
   * false if it is acoustic only.
   * @return True if the model is elastic, false if acoustic only
   */
  PROXY_HOST_DEVICE
  virtual bool isElastic() const = 0;
};

}  // namespace model
#endif  // SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_
