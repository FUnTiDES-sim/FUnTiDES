#ifndef SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_
#define SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_

#include <data_type.h>
#include <sem_macros.h>

namespace model
{

template <typename FloatType, typename ScalarType>
struct ModelDataBase
{
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
  InteriorNode = 0,
  Damping = 1 << 0,
  Sponge = 1 << 1,
  Surface = 1 << 2,
  Ghost = 1 << 3
};

/**
 * @enum CubicFace
 * @brief Local face identifiers for cubic elements
 * Convention: even=minus, odd=plus
 */
enum class CubicFace : int
{
  kXMinus = 0,
  kXPlus = 1,
  kYMinus = 2,
  kYPlus = 3,
  kZMinus = 4,
  kZPlus = 5
};

/**
 * @brief Abstract base class representing a 3D mesh.
 */
template <typename FloatType, typename ScalarType>
class ModelApi
{
 public:
  PROXY_HOST_DEVICE ModelApi() = default;
  PROXY_HOST_DEVICE ModelApi(const ModelDataBase<ScalarType, FloatType>& data)
  {
  }
  PROXY_HOST_DEVICE ModelApi(const ModelApi&) = default;
  PROXY_HOST_DEVICE ModelApi& operator=(const ModelApi&) = default;
  PROXY_HOST_DEVICE ~ModelApi() = default;

  PROXY_HOST_DEVICE
  virtual FloatType nodeCoord(ScalarType dofGlobal, int dim) const = 0;

  PROXY_HOST_DEVICE
  virtual ScalarType globalNodeIndex(ScalarType e, int i, int j,
                                     int k) const = 0;

  PROXY_HOST_DEVICE virtual FloatType getModelVpOnNodes(ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelVpOnElement(
      ScalarType e) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelRhoOnNodes(
      ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelRhoOnElement(
      ScalarType e) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelVsOnNodes(ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelVsOnElement(
      ScalarType e) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelDeltaOnNodes(
      ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelDeltaOnElement(
      ScalarType e) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelEpsilonOnNodes(
      ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelEpsilonOnElement(
      ScalarType e) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelGammaOnNodes(
      ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getModelGammaOnElement(
      ScalarType e) const = 0;
  PROXY_HOST_DEVICE virtual ScalarType getModelThetaOnNodes(
      ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual ScalarType getModelThetaOnElement(
      ScalarType e) const = 0;
  PROXY_HOST_DEVICE virtual ScalarType getModelPhiOnNodes(
      ScalarType n) const = 0;
  PROXY_HOST_DEVICE virtual ScalarType getModelPhiOnElement(
      ScalarType e) const = 0;

  PROXY_HOST_DEVICE
  virtual void getCTensorOnElement(ScalarType e,
                                   FloatType CTTI[6][6]) const = 0;

  virtual void initElasticityTensors() = 0;

  PROXY_HOST_DEVICE virtual ScalarType getNumberOfElements() const = 0;
  PROXY_HOST_DEVICE virtual ScalarType getNumberOfNodes() const = 0;
  PROXY_HOST_DEVICE virtual int getNumberOfPointsPerElement() const = 0;
  PROXY_HOST_DEVICE virtual int getOrder() const = 0;
  PROXY_HOST_DEVICE virtual BoundaryFlag boundaryType(ScalarType n) const = 0;

  /**
   * @brief Compute the outward unit normal vector of an element face.
   * @param e Element index
   * @param face Face identifier (using CubicFace enum)
   * @param[out] v Output array (size 3) holding the normal vector
   */
  PROXY_HOST_DEVICE
  virtual void faceNormal(ScalarType e, CubicFace face,
                          FloatType v[3]) const = 0;

  PROXY_HOST_DEVICE virtual FloatType domainSize(int dim) const = 0;
  PROXY_HOST_DEVICE virtual FloatType getMinSpacing() const = 0;
  virtual FloatType getMaxSpeed() const = 0;
  PROXY_HOST_DEVICE virtual bool isModelOnNodes() const = 0;
  PROXY_HOST_DEVICE virtual bool isElastic() const = 0;

  /**
   * @brief Build face connectivity for absorbing boundary conditions
   */
  virtual void buildFaceConnectivity() = 0;
};

}  // namespace model
#endif  // SRC_MODEL_MODELAPI_INCLUDE_MODEL_API_H_
