#ifndef FUNTIDES_MODEL_MESH_API_INCLUDE_MODEL_H_
#define FUNTIDES_MODEL_MESH_API_INCLUDE_MODEL_H_
#include "data_type.h"
#include "sem_macros.h"

/**
 * @namespace model
 * @brief Mesh geometry, connectivity and material model representations.
 */
namespace model {

/**
 * @brief Empty base of the initialization data passed to ModelApi implementations.
 *
 * @tparam FloatType Floating-point type of coordinates and material values.
 * @tparam ScalarType Integer type of node and element indices.
 */
template <typename FloatType, typename ScalarType>
struct ModelDataBase {
  PROXY_HOST_DEVICE ModelDataBase() = default;
  PROXY_HOST_DEVICE ~ModelDataBase() = default;
  /// Copyable on host and device.
  PROXY_HOST_DEVICE ModelDataBase(const ModelDataBase&) = default;
  /// Copyable on host and device.
  PROXY_HOST_DEVICE ModelDataBase& operator=(const ModelDataBase&) = default;
};

/**
 * @brief Classification of a mesh node with respect to the domain boundaries.
 *
 * Each node carries exactly one value; the values are not combinable bit flags.
 */
enum BoundaryFlag : int {
  InteriorNode = 0,  ///< Node inside the domain
  Damping = 1,       ///< Node on a domain boundary with a damping condition
  Sponge = 2,        ///< Never assigned by the current code
  Surface = 3,       ///< Node on a free surface
  Ghost = 4          ///< Never assigned by the current code
};

/**
 * @enum CubicFace
 * @brief Local face identifiers of a hexahedral element.
 *
 * value / 2 is the normal axis (0 = x, 1 = y, 2 = z) and value % 2 the side (0 = minus, 1 = plus),
 * so value ^ 1 is the opposite face. Normals point outward from the element.
 * @see docs/design.md, section "Hexahedron local numbering".
 */
enum class CubicFace : int {
  kXMinus = 0,  ///< Face at x = x_min (left face, normal = [-1, 0, 0])
  kXPlus = 1,   ///< Face at x = x_max (right face, normal = [+1, 0, 0])
  kYMinus = 2,  ///< Face at y = y_min (front face, normal = [0, -1, 0])
  kYPlus = 3,   ///< Face at y = y_max (back face, normal = [0, +1, 0])
  kZMinus = 4,  ///< Face at z = z_min (bottom face, normal = [0, 0, -1])
  kZPlus = 5    ///< Face at z = z_max (top face, normal = [0, 0, +1])
};

/**
 * @enum AnisotropyType
 * @brief Anisotropy class of the medium, which selects how elasticity tensors are obtained
 * (see ModelApi::initElasticityTensors()).
 */
enum AnisotropyType : uint8_t {
  kIso = 0,       ///< Isotropic medium
  kVTI = 1 << 0,  ///< Vertically Transverse Isotropic medium
  kTTI = 1 << 1   ///< Tilted Transverse Isotropic medium, symmetry axis given by theta and phi
};

/**
 * @brief Abstract interface to a 3D hexahedral spectral-element mesh and the material model
 * defined on it.
 *
 * Material properties are stored either per node or per element, see isModelOnNodes(). Methods marked PROXY_HOST_DEVICE
 * may be called from Kokkos kernels through the concrete type; the others are host only.
 * @see docs/design.md, sections "Hexahedron local numbering" and "Device calls on mesh objects".
 *
 * @todo VERIFY: are all lengths (coordinates, domainSize(), getMinSpacing()) in meters, as
 * implied by velocities in m/s?
 *
 * @tparam FloatType Floating-point type of coordinates and material values.
 * @tparam ScalarType Integer type of node, element and face indices.
 */
template <typename FloatType, typename ScalarType>
class ModelApi {
 public:
  PROXY_HOST_DEVICE ModelApi() = default;

  /**
   * @brief Construct from initialization data; the base class ignores the argument.
   */
  PROXY_HOST_DEVICE ModelApi(const ModelDataBase<ScalarType, FloatType>& data) {}

  /// Copyable on host and device.
  PROXY_HOST_DEVICE ModelApi(const ModelApi&) = default;
  /// Copyable on host and device.
  PROXY_HOST_DEVICE ModelApi& operator=(const ModelApi&) = default;
  PROXY_HOST_DEVICE ~ModelApi() = default;

  /**
   * @brief Get the coordinate of a node along one axis.
   * @param[in] dofGlobal Global node index.
   * @param[in] dim Axis index (0 = x, 1 = y, 2 = z).
   * @return Coordinate value. Whether the subdomain origin is included differs between
   * implementations, see docs/design-red-flags.md.
   */
  PROXY_HOST_DEVICE
  virtual FloatType nodeCoord(ScalarType dofGlobal, int dim) const = 0;

  /**
   * @brief Get the global node index of a node of an element.
   *
   * Nodes shared by adjacent elements have the same global index.
   * @param[in] e Element index.
   * @param[in] i Local index along x, in [0, order].
   * @param[in] j Local index along y, in [0, order].
   * @param[in] k Local index along z, in [0, order].
   * @return Global node index.
   */
  PROXY_HOST_DEVICE
  virtual ScalarType globalNodeIndex(ScalarType e, int i, int j, int k) const = 0;

  /**
   * @brief Get the boundary classification of a node.
   * @param[in] n Global node index.
   * @return BoundaryFlag of the node; InteriorNode when the mesh carries no node classification.
   */
  PROXY_HOST_DEVICE
  virtual BoundaryFlag boundaryType(ScalarType n) const = 0;

  /**
   * @brief Get the P-wave velocity at a node.
   * @param[in] n Global node index.
   * @return Vp in m/s.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVpOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the P-wave velocity of an element.
   * @param[in] e Element index.
   * @return Vp in m/s.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVpOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the density at a node.
   * @param[in] n Global node index.
   * @return Density in kg/m^3.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelRhoOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the density of an element.
   * @param[in] e Element index.
   * @return Density in kg/m^3.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelRhoOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the S-wave velocity at a node.
   * @param[in] n Global node index.
   * @return Vs in m/s.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVsOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the S-wave velocity of an element.
   * @param[in] e Element index.
   * @return Vs in m/s.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelVsOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the P-wave quality factor (attenuation) at a node.
   * @param[in] n Global node index.
   * @return Qp, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelQpOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the P-wave quality factor (attenuation) of an element.
   * @param[in] e Element index.
   * @return Qp, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelQpOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the S-wave quality factor (attenuation) at a node.
   * @param[in] n Global node index.
   * @return Qs, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelQsOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the S-wave quality factor (attenuation) of an element.
   * @param[in] e Element index.
   * @return Qs, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelQsOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the Thomsen delta anisotropy parameter at a node.
   * @param[in] n Global node index.
   * @return Delta, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelDeltaOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the Thomsen delta anisotropy parameter of an element.
   * @param[in] e Element index.
   * @return Delta, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelDeltaOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the Thomsen epsilon anisotropy parameter at a node.
   * @param[in] n Global node index.
   * @return Epsilon, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelEpsilonOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the Thomsen epsilon anisotropy parameter of an element.
   * @param[in] e Element index.
   * @return Epsilon, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelEpsilonOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the Thomsen gamma anisotropy parameter at a node.
   * @param[in] n Global node index.
   * @return Gamma, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelGammaOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the Thomsen gamma anisotropy parameter of an element.
   * @param[in] e Element index.
   * @return Gamma, dimensionless.
   */
  PROXY_HOST_DEVICE virtual FloatType getModelGammaOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the tilt angle theta of the TTI symmetry axis at a node.
   *
   * The ScalarType return type truncates the angle to an integer, see docs/design-red-flags.md.
   * @todo VERIFY: theta is measured from which axis (z?), and in degrees or radians?
   * @param[in] n Global node index.
   * @return Theta, truncated to an integer.
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelThetaOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the tilt angle theta of the TTI symmetry axis of an element.
   * @param[in] e Element index.
   * @return Theta, truncated to an integer, same convention as getModelThetaOnNodes().
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelThetaOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the azimuth angle phi of the TTI symmetry axis at a node.
   *
   * The ScalarType return type truncates the angle to an integer, see docs/design-red-flags.md.
   * @todo VERIFY: phi is measured from which axis (x?), in which plane, in degrees or radians?
   * @param[in] n Global node index.
   * @return Phi, truncated to an integer.
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelPhiOnNodes(ScalarType n) const = 0;

  /**
   * @brief Get the azimuth angle phi of the TTI symmetry axis of an element.
   * @param[in] e Element index.
   * @return Phi, truncated to an integer, same convention as getModelPhiOnNodes().
   */
  PROXY_HOST_DEVICE virtual ScalarType getModelPhiOnElement(ScalarType e) const = 0;

  /**
   * @brief Get the precomputed elasticity tensor of an element.
   *
   * Precondition: initElasticityTensors(kTTI) was called on an elastic model; otherwise no
   * tensor is stored.
   * @todo VERIFY: Voigt index order (xx, yy, zz, yz, xz, xy?) and unit (Pa, i.e. rho * Vp^2?).
   * @param[in] e Element index.
   * @param[out] CTTI 6x6 elasticity tensor in Voigt notation.
   */
  PROXY_HOST_DEVICE
  virtual void getCTensorOnElement(ScalarType e, FloatType CTTI[6][6]) const = 0;

  /**
   * @brief Precompute the per-element elasticity tensors read by getCTensorOnElement().
   *
   * Host only. Call once the material properties are final. Tensors are stored only for kTTI
   * on an elastic model; in every other case this call does nothing.
   */
  virtual void initElasticityTensors(AnisotropyType anisotropy_type) = 0;

  /**
   * @brief Get the number of elements held by this model.
   */
  PROXY_HOST_DEVICE virtual ScalarType getNumberOfElements() const = 0;

  /**
   * @brief Get the number of distinct nodes held by this model.
   */
  PROXY_HOST_DEVICE virtual ScalarType getNumberOfNodes() const = 0;

  /**
   * @brief Get the number of GLL nodes per element.
   * @return (order+1)^3.
   */
  PROXY_HOST_DEVICE virtual int getNumberOfPointsPerElement() const = 0;

  /**
   * @brief Get the polynomial order of the spectral element basis.
   */
  PROXY_HOST_DEVICE virtual int getOrder() const = 0;

  /**
   * @brief Compute the outward unit normal of an element face.
   * @param[in] e Element index.
   * @param[in] face Local face of element e.
   * @param[out] v Normal vector [nx, ny, nz], size 3.
   */
  PROXY_HOST_DEVICE
  virtual void faceNormal(ScalarType e, CubicFace face, FloatType v[3]) const = 0;

  /**
   * @brief Get the extent of the domain held by this model along one axis (the local subdomain
   * in distributed runs).
   * @param[in] dim Axis index (0 = x, 1 = y, 2 = z).
   */
  PROXY_HOST_DEVICE virtual FloatType domainSize(int dim) const = 0;

  /**
   * @brief Get the smallest distance between two adjacent nodes of the mesh.
   *
   * Not every implementation scans the whole mesh, see docs/design-red-flags.md.
   */
  PROXY_HOST_DEVICE virtual FloatType getMinSpacing() const = 0;

  /**
   * @brief Get the maximum P-wave velocity of the model, in m/s. Host only.
   *
   * Not every implementation reads the stored velocities, see docs/design-red-flags.md.
   */
  virtual FloatType getMaxSpeed() const = 0;

  /**
   * @brief Tell where material properties are stored.
   * @return true for one value per node (OnNodes getters), false for one value per element
   * (OnElement getters).
   */
  PROXY_HOST_DEVICE virtual bool isModelOnNodes() const = 0;

  /**
   * @brief Tell whether the model carries shear properties.
   * @return true for elastic models, false for acoustic models.
   */
  PROXY_HOST_DEVICE virtual bool isElastic() const = 0;

  /**
   * @brief Build the face tables used by getNumberOfFaces(), getGlobalFace(),
   * getGlobalNodeFromFace() and isBoundaryFace().
   *
   * Host only. Must be called before those queries; calling it again is harmless.
   */
  virtual void buildFaceConnectivity() = 0;

  /**
   * @brief Tell whether a node lies on a free surface.
   * @param[in] n Global node index.
   * @return true if boundaryType(n) is Surface.
   */
  PROXY_HOST_DEVICE
  virtual bool isFreeSurface(ScalarType n) const = 0;

  /**
   * @brief Set uniform quality factors on every element, replacing any per-element values.
   *
   * Host only. Per-node quality factors are left unchanged.
   * @param[in] qp P-wave quality factor.
   * @param[in] qs S-wave quality factor, meaningful for elastic models only.
   */
  virtual void setQualityFactors(FloatType qp, FloatType qs) = 0;

  /**
   * @brief Get the number of distinct faces of the mesh. Requires buildFaceConnectivity().
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getNumberOfFaces() const = 0;

  /**
   * @brief Tell whether a face lies on the boundary of the domain.
   *
   * When nodes carry a boundary classification, a face is a boundary face if none of its nodes
   * is InteriorNode; otherwise it is a boundary face if it has no neighbor element. This differs
   * from FaceConnectivityApi::isBoundaryFace(), which is purely topological.
   * @param[in] face_id Global face id.
   */
  PROXY_HOST_DEVICE
  virtual bool isBoundaryFace(ScalarType face_id) const = 0;

  /**
   * @brief Get the global face id of a local face of an element.
   * @param[in] elem Element index.
   * @param[in] local_face Local face of elem.
   * @return Global face id, in [0, getNumberOfFaces()).
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const = 0;

  /**
   * @brief Get the global node index of a node of a face.
   * @param[in] face_global Global face id.
   * @param[in] local_dof 2D face DOF index, in [0, (order+1)^2).
   * @return Global node index.
   * @see docs/design.md, section "Hexahedron local numbering".
   */
  PROXY_HOST_DEVICE
  virtual ScalarType getGlobalNodeFromFace(ScalarType face_global, int local_dof) const = 0;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_API_INCLUDE_MODEL_H_
