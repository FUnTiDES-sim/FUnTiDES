#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_MODEL_UNSTRUCT_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_MODEL_UNSTRUCT_H_
#include <elasticity_utils.h>
#include <model.h>

#include <algorithm>
#include <array>
#include <map>

#include "face_connectivity_unstruct.h"

namespace model {

/**
 * @brief Plain data holder used to build a ModelUnstruct.
 *
 * Per-node arrays have one entry per global node, per-element arrays one entry per element.
 * For each material property only the node or the element variant is expected to be filled,
 * according to isModelOnNodes_; the other one may be empty.
 * @tparam FloatType Floating-point type of the geometry and material scalars.
 * @tparam ScalarType Integer type of counts and indices.
 */
template <typename FloatType, typename ScalarType>
struct ModelUnstructData : public ModelDataBase<FloatType, ScalarType> {
  PROXY_HOST_DEVICE ModelUnstructData() = default;
  PROXY_HOST_DEVICE ~ModelUnstructData() = default;
  PROXY_HOST_DEVICE ModelUnstructData(const ModelUnstructData&) = default;
  PROXY_HOST_DEVICE ModelUnstructData& operator=(const ModelUnstructData&) = default;

  /**
   * @brief Full constructor. Every array is stored as given (shallow copy of the views).
   * @param order Polynomial order of the elements.
   * @param n_element Number of local elements.
   * @param n_node Number of local nodes.
   * @param lx Local domain length along x.
   * @param ly Local domain length along y.
   * @param lz Local domain length along z.
   * @param isModelOnNodes True if material properties are given per node, false if per element.
   * @param isElastic True for elastic propagation, false for acoustic.
   * @param global_node_index Node table, shape (n_element, (order+1)^3), indexed by (e, i + j*(order+1) + k*(order+1)^2).
   * @param nodes_coords_x X coordinate of each node (size n_node).
   * @param nodes_coords_y Y coordinate of each node (size n_node).
   * @param nodes_coords_z Z coordinate of each node (size n_node).
   * @param boundaries_t BoundaryFlag value of each node (size n_node), may be empty.
   * @param model_C_tensor_element Per-element 6x6 Voigt elasticity tensor, shape (n_element, 6, 6).
   * @param face_connectivity Optional precomputed face connectivity.
   * @todo VERIFY: units of the angles model_theta_* and model_phi_* (radians or degrees).
   */
  PROXY_HOST_DEVICE
  ModelUnstructData(ScalarType order, ScalarType n_element, ScalarType n_node, FloatType lx, FloatType ly, FloatType lz,
                    bool isModelOnNodes, bool isElastic, arrayInt global_node_index, vectorReal nodes_coords_x,
                    vectorReal nodes_coords_y, vectorReal nodes_coords_z, vectorReal model_vp_node,
                    vectorReal model_vp_element, vectorReal model_rho_node, vectorReal model_rho_element,
                    vectorReal model_vs_node, vectorReal model_vs_element, vectorReal model_delta_node,
                    vectorReal model_delta_element, vectorReal model_epsilon_node, vectorReal model_epsilon_element,
                    vectorReal model_gamma_node, vectorReal model_gamma_element, vectorReal model_theta_node,
                    vectorReal model_theta_element, vectorReal model_phi_node, vectorReal model_phi_element,
                    array3DReal model_C_tensor_element, vectorInt boundaries_t, vectorReal model_qp_node = vectorReal(),
                    vectorReal model_qp_element = vectorReal(), vectorReal model_qs_node = vectorReal(),
                    vectorReal model_qs_element = vectorReal(),
                    FaceConnectivityUnstructData<FloatType, ScalarType> face_connectivity = {})
      : order_(order),
        n_element_(n_element),
        n_node_(n_node),
        lx_(lx),
        ly_(ly),
        lz_(lz),
        isModelOnNodes_(isModelOnNodes),
        isElastic_(isElastic),
        global_node_index_(global_node_index),
        nodes_coords_x_(nodes_coords_x),
        nodes_coords_y_(nodes_coords_y),
        nodes_coords_z_(nodes_coords_z),
        model_vp_node_(model_vp_node),
        model_vp_element_(model_vp_element),
        model_rho_node_(model_rho_node),
        model_rho_element_(model_rho_element),
        model_vs_node_(model_vs_node),
        model_vs_element_(model_vs_element),
        model_qp_node_(model_qp_node),
        model_qp_element_(model_qp_element),
        model_qs_node_(model_qs_node),
        model_qs_element_(model_qs_element),
        model_delta_node_(model_delta_node),
        model_delta_element_(model_delta_element),
        model_epsilon_node_(model_epsilon_node),
        model_epsilon_element_(model_epsilon_element),
        model_gamma_node_(model_gamma_node),
        model_gamma_element_(model_gamma_element),
        model_theta_node_(model_theta_node),
        model_theta_element_(model_theta_element),
        model_phi_node_(model_phi_node),
        model_phi_element_(model_phi_element),
        model_C_tensor_element_(model_C_tensor_element),
        boundaries_t_(boundaries_t),
        face_connectivity_(face_connectivity) {}

  /// @todo VERIFY: is origin_x_/y_/z_ read anywhere, and how does it relate to ox_/oy_/oz_?
  FloatType origin_x_{0}, origin_y_{0}, origin_z_{0};
  FloatType ox_, oy_, oz_;  ///< Origin of the local subdomain.
  ScalarType order_;        ///< Polynomial order of the elements.
  ScalarType n_element_;    ///< Number of local elements.
  ScalarType n_node_;       ///< Number of local nodes.
  FloatType lx_, ly_, lz_;  ///< Local domain lengths.

  bool isModelOnNodes_;  ///< True if material properties are stored per node, false if per element.
  bool isElastic_;       ///< True for elastic propagation, false for acoustic.

  arrayInt global_node_index_;  ///< Shape (n_element, (order+1)^3), indexed by (e, i + j*(order+1) + k*(order+1)^2).
  vectorReal nodes_coords_x_;   ///< X coordinate of each node, size n_node.
  vectorReal nodes_coords_y_;   ///< Y coordinate of each node, size n_node.
  vectorReal nodes_coords_z_;   ///< Z coordinate of each node, size n_node.

  vectorReal model_vp_node_;
  vectorReal model_vp_element_;
  vectorReal model_rho_node_;
  vectorReal model_rho_element_;
  vectorReal model_vs_node_;
  vectorReal model_vs_element_;
  vectorReal model_qp_node_;
  vectorReal model_qp_element_;
  vectorReal model_qs_node_;
  vectorReal model_qs_element_;
  vectorReal model_delta_node_;
  vectorReal model_delta_element_;
  vectorReal model_epsilon_node_;
  vectorReal model_epsilon_element_;
  vectorReal model_gamma_node_;
  vectorReal model_gamma_element_;
  vectorReal model_theta_node_;
  vectorReal model_theta_element_;
  vectorReal model_phi_node_;
  vectorReal model_phi_element_;
  array3DReal model_C_tensor_element_;  ///< Shape (n_element, 6, 6), Voigt notation.
  vectorInt boundaries_t_;              ///< BoundaryFlag value of each node, size n_node; may be empty.
  FaceConnectivityUnstructData<FloatType, ScalarType> face_connectivity_;
};

/**
 * @brief Unstructured mesh of hexahedral elements with an explicit node table and node coordinates.
 *
 * Holds the geometry, the material properties (per node or per element) and the face
 * connectivity of the local subdomain. Copies share the underlying Kokkos views.
 * @tparam FloatType Floating-point type of the geometry and material scalars.
 * @tparam ScalarType Integer type of counts and indices.
 */
template <typename FloatType, typename ScalarType>
class ModelUnstruct final : public ModelApi<FloatType, ScalarType> {
 public:
  using IndexType = int;

  PROXY_HOST_DEVICE ModelUnstruct() = default;

  /**
   * @brief Construct from a data structure (shallow copy of the views).
   * @param data Mesh and material data.
   */
  PROXY_HOST_DEVICE ModelUnstruct(const ModelUnstructData<FloatType, ScalarType>& data)
      : order_(data.order_),
        n_element_(data.n_element_),
        n_node_(data.n_node_),
        lx_(data.lx_),
        ly_(data.ly_),
        lz_(data.lz_),
        ox_(data.ox_),
        oy_(data.oy_),
        oz_(data.oz_),
        isModelOnNodes_(data.isModelOnNodes_),
        isElastic_(data.isElastic_),
        global_node_index_(data.global_node_index_),
        nodes_coords_x_(data.nodes_coords_x_),
        nodes_coords_y_(data.nodes_coords_y_),
        nodes_coords_z_(data.nodes_coords_z_),
        model_vp_node_(data.model_vp_node_),
        model_vp_element_(data.model_vp_element_),
        model_rho_node_(data.model_rho_node_),
        model_rho_element_(data.model_rho_element_),
        model_vs_node_(data.model_vs_node_),
        model_vs_element_(data.model_vs_element_),
        model_qp_node_(data.model_qp_node_),
        model_qp_element_(data.model_qp_element_),
        model_qs_node_(data.model_qs_node_),
        model_qs_element_(data.model_qs_element_),
        model_delta_node_(data.model_delta_node_),
        model_delta_element_(data.model_delta_element_),
        model_epsilon_node_(data.model_epsilon_node_),
        model_epsilon_element_(data.model_epsilon_element_),
        model_gamma_node_(data.model_gamma_node_),
        model_gamma_element_(data.model_gamma_element_),
        model_phi_node_(data.model_phi_node_),
        model_phi_element_(data.model_phi_element_),
        model_theta_node_(data.model_theta_node_),
        model_theta_element_(data.model_theta_element_),
        model_C_tensor_element_(data.model_C_tensor_element_),
        boundaries_t_(data.boundaries_t_),
        face_connectivity_(data.face_connectivity_),
        n_points_per_element_((order_ + 1) * (order_ + 1) * (order_ + 1)) {}

  PROXY_HOST_DEVICE ModelUnstruct& operator=(const ModelUnstruct&) = default;
  PROXY_HOST_DEVICE ~ModelUnstruct() = default;

  /**
   * @brief Convert a linear element index to an element index.
   * @param linearIndex Linear element index.
   * @return linearIndex unchanged.
   */
  PROXY_HOST_DEVICE
  IndexType elementIndex(const int linearIndex) const { return linearIndex; }

  /**
   * @brief Get the global node index of a vertex of an element.
   * @param e Element index.
   * @param i Vertex coordinate along the first local axis, 0 or 1.
   * @param j Vertex coordinate along the second local axis, 0 or 1.
   * @param k Vertex coordinate along the third local axis, 0 or 1.
   * @return Global node index of the corner node.
   */
  PROXY_HOST_DEVICE
  IndexType globalVertexIndex(IndexType e, int const i, int const j, int const k) const {
    int local_i = i * order_;
    int local_j = j * order_;
    int local_k = k * order_;
    const auto localDofIndex = local_i + local_j * (order_ + 1) + local_k * (order_ + 1) * (order_ + 1);
    return global_node_index_(e, localDofIndex);
  }

  /**
   * @brief Get the coordinates of a node.
   * @param dofGlobal Global node index.
   * @param[out] coords Array of at least 3 entries, receives x, y, z.
   */
  PROXY_HOST_DEVICE
  void vertexCoords(IndexType dofGlobal, FloatType* const coords) const {
    coords[0] = nodes_coords_x_[dofGlobal];
    coords[1] = nodes_coords_y_[dofGlobal];
    coords[2] = nodes_coords_z_[dofGlobal];
  }

  /**
   * @brief Get one coordinate of a node.
   * @param dofGlobal Global node index.
   * @param dim Axis: 0 = x, 1 = y, 2 = z.
   * @return Coordinate value, or -1 if dim is not in 0..2.
   */
  PROXY_HOST_DEVICE
  FloatType nodeCoord(ScalarType dofGlobal, int dim) const final {
    switch (dim) {
      case 0: {
        return nodes_coords_x_[dofGlobal];
      }
      case 1: {
        return nodes_coords_y_[dofGlobal];
      }
      case 2: {
        return nodes_coords_z_[dofGlobal];
      }
      default:
        return FloatType(-1);
    }
  }

  /**
   * @brief Get the global node index of a local node of an element.
   * @param e Element index.
   * @param i Local index along the first axis, in [0, order].
   * @param j Local index along the second axis, in [0, order].
   * @param k Local index along the third axis, in [0, order].
   * @return Global node index. The local node number is i + j*(order+1) + k*(order+1)^2.
   */
  PROXY_HOST_DEVICE
  ScalarType globalNodeIndex(ScalarType e, int i, int j, int k) const final {
    const auto localDofIndex = i + j * (order_ + 1) + k * (order_ + 1) * (order_ + 1);
    return global_node_index_(e, localDofIndex);
  }

  /**
   * @brief Get the P-wave velocity at a node.
   * @param n Node index.
   * @return P-wave velocity (m/s).
   */
  PROXY_HOST_DEVICE FloatType getModelVpOnNodes(ScalarType n) const final { return model_vp_node_[n]; }

  /**
   * @brief Get the P-wave velocity of an element.
   * @param e Element index.
   * @return P-wave velocity (m/s).
   */
  PROXY_HOST_DEVICE FloatType getModelVpOnElement(ScalarType e) const final { return model_vp_element_[e]; }

  /**
   * @brief Get the density at a node.
   * @param n Node index.
   * @return Density (kg/m^3).
   */
  PROXY_HOST_DEVICE FloatType getModelRhoOnNodes(ScalarType n) const final { return model_rho_node_[n]; }

  /**
   * @brief Get the density of an element.
   * @param e Element index.
   * @return Density (kg/m^3).
   */
  PROXY_HOST_DEVICE FloatType getModelRhoOnElement(ScalarType e) const final { return model_rho_element_[e]; }

  /**
   * @brief Get the S-wave velocity at a node.
   * @param n Node index.
   * @return S-wave velocity (m/s).
   */
  PROXY_HOST_DEVICE FloatType getModelVsOnNodes(ScalarType n) const final { return model_vs_node_[n]; }

  /**
   * @brief Get the S-wave velocity of an element.
   * @param e Element index.
   * @return S-wave velocity (m/s).
   */
  PROXY_HOST_DEVICE FloatType getModelVsOnElement(ScalarType e) const final { return model_vs_element_[e]; }

  /**
   * @brief Overwrite the per-node vp, vs and rho at a node.
   *
   * Host only. The per-node arrays must be allocated.
   * @param n Global node index.
   * @param vp P-wave velocity (m/s).
   * @param vs S-wave velocity (m/s).
   * @param rho Density (kg/m^3).
   */
  void setModelNodeProps(ScalarType n, FloatType vp, FloatType vs, FloatType rho) {
    model_vp_node_[n] = vp;
    model_vs_node_[n] = vs;
    model_rho_node_[n] = rho;
  }

  /**
   * @brief Get the P-wave quality factor at a node.
   * @param n Node index.
   * @return Quality factor Qp (dimensionless), or 1.0e9 if no per-node array is stored.
   */
  PROXY_HOST_DEVICE FloatType getModelQpOnNodes(ScalarType n) const final {
    if (model_qp_node_.extent(0) > 0) return model_qp_node_[n];
    return static_cast<FloatType>(1.0e9);
  }

  /**
   * @brief Get the P-wave quality factor of an element.
   * @param e Element index.
   * @return Quality factor Qp (dimensionless), or 1.0e9 if no per-element array is stored.
   */
  PROXY_HOST_DEVICE FloatType getModelQpOnElement(ScalarType e) const final {
    if (model_qp_element_.extent(0) > 0) return model_qp_element_[e];
    return static_cast<FloatType>(1.0e9);
  }

  /**
   * @brief Get the S-wave quality factor at a node.
   * @param n Node index.
   * @return Quality factor Qs (dimensionless), or 1.0e9 if no per-node array is stored.
   */
  PROXY_HOST_DEVICE FloatType getModelQsOnNodes(ScalarType n) const final {
    if (model_qs_node_.extent(0) > 0) return model_qs_node_[n];
    return static_cast<FloatType>(1.0e9);
  }

  /**
   * @brief Get the S-wave quality factor of an element.
   * @param e Element index.
   * @return Quality factor Qs (dimensionless), or 1.0e9 if no per-element array is stored.
   */
  PROXY_HOST_DEVICE FloatType getModelQsOnElement(ScalarType e) const final {
    if (model_qs_element_.extent(0) > 0) return model_qs_element_[e];
    return static_cast<FloatType>(1.0e9);
  }

  /**
   * @brief Get the Thomsen delta parameter at a node.
   * @param n Node index.
   * @return Thomsen delta (dimensionless).
   */
  PROXY_HOST_DEVICE FloatType getModelDeltaOnNodes(ScalarType n) const final { return model_delta_node_[n]; }

  /**
   * @brief Get the Thomsen delta parameter of an element.
   * @param e Element index.
   * @return Thomsen delta (dimensionless).
   */
  PROXY_HOST_DEVICE FloatType getModelDeltaOnElement(ScalarType e) const final { return model_delta_element_[e]; }

  /**
   * @brief Get the Thomsen epsilon parameter at a node.
   * @param n Node index.
   * @return Thomsen epsilon (dimensionless).
   */
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnNodes(ScalarType n) const final { return model_epsilon_node_[n]; }

  /**
   * @brief Get the Thomsen epsilon parameter of an element.
   * @param e Element index.
   * @return Thomsen epsilon (dimensionless).
   */
  PROXY_HOST_DEVICE FloatType getModelEpsilonOnElement(ScalarType e) const final { return model_epsilon_element_[e]; }

  /**
   * @brief Get the Thomsen gamma parameter at a node.
   * @param n Node index.
   * @return Thomsen gamma (dimensionless).
   */
  PROXY_HOST_DEVICE FloatType getModelGammaOnNodes(ScalarType n) const final { return model_gamma_node_[n]; }

  /**
   * @brief Get the Thomsen gamma parameter of an element.
   * @param e Element index.
   * @return Thomsen gamma (dimensionless).
   */
  PROXY_HOST_DEVICE FloatType getModelGammaOnElement(ScalarType e) const final { return model_gamma_element_[e]; }

  /**
   * @brief Get the azimuth angle phi at a node.
   * @param n Node index.
   * @return Azimuth angle phi.
   * @todo VERIFY: unit of the stored angle (radians here, degrees in computeCTensor).
   */
  PROXY_HOST_DEVICE ScalarType getModelPhiOnNodes(ScalarType n) const final { return model_phi_node_[n]; }

  /**
   * @brief Get the azimuth angle phi of an element.
   * @param e Element index.
   * @return Azimuth angle phi.
   * @todo VERIFY: unit of the stored angle (radians here, degrees in computeCTensor).
   */
  PROXY_HOST_DEVICE ScalarType getModelPhiOnElement(ScalarType e) const final { return model_phi_element_[e]; }

  /**
   * @brief Get the tilt angle theta at a node.
   * @param n Node index.
   * @return Tilt angle theta.
   * @todo VERIFY: unit of the stored angle (radians here, degrees in computeCTensor).
   */
  PROXY_HOST_DEVICE ScalarType getModelThetaOnNodes(ScalarType n) const final { return model_theta_node_[n]; }

  /**
   * @brief Get the tilt angle theta of an element.
   * @param e Element index.
   * @return Tilt angle theta.
   * @todo VERIFY: unit of the stored angle (radians here, degrees in computeCTensor).
   */
  PROXY_HOST_DEVICE ScalarType getModelThetaOnElement(ScalarType e) const final { return model_theta_element_[e]; }

  /**
   * @brief Precompute the per-element 6x6 Voigt elasticity tensors.
   *
   * Does nothing if the mesh is not elastic. Only the TTI case allocates and fills
   * model_C_tensor_element_ from the per-element vp, vs, rho and Thomsen parameters.
   * @param anisotropy_type Anisotropy of the medium.
   */
  void initElasticityTensors(AnisotropyType anisotropy_type) final {
    if (!isElastic_) return;

    if (anisotropy_type == AnisotropyType::kIso || anisotropy_type == AnisotropyType::kVTI) {
      // Isotropic and VTI tensors are not stored: the solver builds them on the fly.
      return;
    }

    if (anisotropy_type == AnisotropyType::kTTI) {
      model_C_tensor_element_ = allocateArray3D<array3DReal>(n_element_, 6, 6);

      auto C_tensor = model_C_tensor_element_;
      auto vp = model_vp_element_;
      auto vs = model_vs_element_;
      auto rho = model_rho_element_;
      auto delta = model_delta_element_;
      auto epsilon = model_epsilon_element_;
      auto gamma = model_gamma_element_;
      auto theta = model_theta_element_;
      auto phi = model_phi_element_;

      Kokkos::parallel_for(
          "Model init ElasticTensors",
          Kokkos::RangePolicy<Kokkos::LaunchBounds<LaunchMaxThreadsPerBlock, LaunchMinBlocksPerSM>>(0, n_element_),
          KOKKOS_LAMBDA(const int i) {
            FloatType CTTI[6][6];
            FloatType vp_val = static_cast<FloatType>(vp[i]);
            FloatType vs_val = static_cast<FloatType>(vs[i]);
            FloatType rho_val = static_cast<FloatType>(rho[i]);
            FloatType delta_val = static_cast<FloatType>(delta[i]);
            FloatType epsilon_val = static_cast<FloatType>(epsilon[i]);
            FloatType gamma_val = static_cast<FloatType>(gamma[i]);
            FloatType theta_val = static_cast<FloatType>(theta[i]);
            FloatType phi_val = static_cast<FloatType>(phi[i]);

            computeCTensor(vp_val, vs_val, rho_val, delta_val, epsilon_val, gamma_val, theta_val, phi_val, CTTI);

            for (int k = 0; k < 6; k++)
              for (int l = 0; l < 6; l++) C_tensor(i, k, l) = CTTI[k][l];
          });
    }
  }

  /**
   * @brief Get the stored elasticity tensor of an element.
   * @param e Element index.
   * @param[out] CTTI 6x6 elasticity tensor in Voigt notation.
   */
  PROXY_HOST_DEVICE
  void getCTensorOnElement(ScalarType e, FloatType CTTI[6][6]) const final {
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) CTTI[i][j] = model_C_tensor_element_(e, i, j);
  }

  /**
   * @brief Tell where the material properties are stored.
   * @return True if per node, false if per element.
   */
  PROXY_HOST_DEVICE bool isModelOnNodes() const final { return isModelOnNodes_; }

  /**
   * @brief Tell the wave equation the mesh is set up for.
   * @return True if elastic, false if acoustic.
   */
  PROXY_HOST_DEVICE bool isElastic() const final { return isElastic_; }

  /**
   * @brief Get the number of local elements.
   * @return Number of elements.
   */
  PROXY_HOST_DEVICE ScalarType getNumberOfElements() const final { return n_element_; }

  /**
   * @brief Get the number of local nodes.
   * @return Number of nodes.
   */
  PROXY_HOST_DEVICE ScalarType getNumberOfNodes() const final { return n_node_; }

  /**
   * @brief Get the number of nodes per element.
   * @return (order+1)^3.
   */
  PROXY_HOST_DEVICE int getNumberOfPointsPerElement() const final { return n_points_per_element_; }

  /**
   * @brief Get the polynomial order of the elements.
   * @return Element order.
   */
  PROXY_HOST_DEVICE int getOrder() const final { return static_cast<int>(order_); }

  /**
   * @brief Compute the unit outward normal of an element face.
   *
   * The normal is built from three corner nodes of the face. It is left unnormalized
   * if its norm is below 1e-12.
   * @param e Element index.
   * @param local_face Local face of the element.
   * @param[out] v Array of 3 entries, receives the normal (x, y, z).
   */
  PROXY_HOST_DEVICE
  void faceNormal(ScalarType e, CubicFace local_face, FloatType v[3]) const final {
    ScalarType n0, n1, n2;
    const int o = order_;

    switch (local_face) {
      case CubicFace::kXMinus:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, 0, o, 0);
        n2 = globalNodeIndex(e, 0, 0, o);
        break;
      case CubicFace::kXPlus:
        n0 = globalNodeIndex(e, o, 0, 0);
        n1 = globalNodeIndex(e, o, 0, o);
        n2 = globalNodeIndex(e, o, o, 0);
        break;
      case CubicFace::kYMinus:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, 0, 0, o);
        n2 = globalNodeIndex(e, o, 0, 0);
        break;
      case CubicFace::kYPlus:
        n0 = globalNodeIndex(e, 0, o, 0);
        n1 = globalNodeIndex(e, o, o, 0);
        n2 = globalNodeIndex(e, 0, o, o);
        break;
      case CubicFace::kZMinus:
        n0 = globalNodeIndex(e, 0, 0, 0);
        n1 = globalNodeIndex(e, o, 0, 0);
        n2 = globalNodeIndex(e, 0, o, 0);
        break;
      case CubicFace::kZPlus:
        n0 = globalNodeIndex(e, 0, 0, o);
        n1 = globalNodeIndex(e, 0, o, o);
        n2 = globalNodeIndex(e, o, 0, o);
        break;
    }

    FloatType p0[3], p1[3], p2[3];
    for (int d = 0; d < 3; ++d) {
      p0[d] = nodeCoord(n0, d);
      p1[d] = nodeCoord(n1, d);
      p2[d] = nodeCoord(n2, d);
    }

    FloatType t1[3], t2[3];
    t1[0] = p1[0] - p0[0];
    t1[1] = p1[1] - p0[1];
    t1[2] = p1[2] - p0[2];
    t2[0] = p2[0] - p0[0];
    t2[1] = p2[1] - p0[1];
    t2[2] = p2[2] - p0[2];

    // t1 x t2 points inward for this vertex ordering on all six faces; negate to face outward.
    v[0] = -(t1[1] * t2[2] - t1[2] * t2[1]);
    v[1] = -(t1[2] * t2[0] - t1[0] * t2[2]);
    v[2] = -(t1[0] * t2[1] - t1[1] * t2[0]);

    FloatType norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (norm > 1e-12) {
      v[0] /= norm;
      v[1] /= norm;
      v[2] /= norm;
    }
  }

  /**
   * @brief Get the boundary flag of a node.
   * @param n Node index.
   * @return BoundaryFlag of the node; InteriorNode if no boundary flags are stored.
   */
  PROXY_HOST_DEVICE
  BoundaryFlag boundaryType(ScalarType n) const override {
    if (boundaries_t_.extent(0) == 0) return BoundaryFlag::InteriorNode;
    return static_cast<BoundaryFlag>(boundaries_t_[n]);
  }

  /**
   * @brief Get the local domain length along an axis.
   * @param dim Axis: 0 = x, 1 = y, 2 = z.
   * @return Domain length (m), or -1 if dim is not in 0..2.
   */
  PROXY_HOST_DEVICE FloatType domainSize(int dim) const final {
    switch (dim) {
      case 0:
        return lx_;
      case 1:
        return ly_;
      case 2:
        return lz_;
      default:
        return FloatType(-1);
    }
  }

  /**
   * @brief Get the minimum distance between adjacent nodes.
   * @return Minimum node spacing (m).
   */
  PROXY_HOST_DEVICE FloatType getMinSpacing() const final {
    FloatType minSpacing = std::numeric_limits<FloatType>::max();
    constexpr ScalarType e = 0;

    for (int k = 0; k <= order_; ++k)
      for (int j = 0; j <= order_; ++j)
        for (int i = 0; i < order_; ++i) {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i + 1, j, k);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          minSpacing = fmin(minSpacing, sqrt(dx * dx + dy * dy + dz * dz));
        }

    for (int k = 0; k <= order_; ++k)
      for (int i = 0; i <= order_; ++i)
        for (int j = 0; j < order_; ++j) {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i, j + 1, k);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          minSpacing = fmin(minSpacing, sqrt(dx * dx + dy * dy + dz * dz));
        }

    for (int j = 0; j <= order_; ++j)
      for (int i = 0; i <= order_; ++i)
        for (int k = 0; k < order_; ++k) {
          ScalarType node1 = globalNodeIndex(e, i, j, k);
          ScalarType node2 = globalNodeIndex(e, i, j, k + 1);
          FloatType dx = nodeCoord(node2, 0) - nodeCoord(node1, 0);
          FloatType dy = nodeCoord(node2, 1) - nodeCoord(node1, 1);
          FloatType dz = nodeCoord(node2, 2) - nodeCoord(node1, 2);
          minSpacing = fmin(minSpacing, sqrt(dx * dx + dy * dy + dz * dz));
        }

    return minSpacing;
  }

  /**
   * @brief Get the maximum P-wave velocity of the local mesh.
   * @return Maximum P-wave velocity (m/s).
   * @throws std::runtime_error if neither per-node nor per-element vp is stored.
   */
  FloatType getMaxSpeed() const final {
    FloatType maxSpeedNode = std::numeric_limits<FloatType>::lowest();
    FloatType maxSpeedElem = std::numeric_limits<FloatType>::lowest();

    if (model_vp_node_.extent(0) > 0) {
      FIND_MAX_1D(model_vp_node_, n_node_, maxSpeedNode);
    } else if (model_vp_element_.extent(0) > 0) {
      FIND_MAX_1D(model_vp_element_, n_element_, maxSpeedElem);
    } else {
      throw std::runtime_error("No model initialized (model unstruct getMaxSpeed).");
    }
    return max(maxSpeedElem, maxSpeedNode);
  }

  /**
   * @brief Build the face connectivity tables.
   *
   * Identifies the unique faces, the element pair sharing each internal face, and the
   * nodes of each face. Does nothing if the tables already exist. Must be called before
   * any face query.
   */
  void buildFaceConnectivity() override {
    if (face_connectivity_.getNumberOfFaces() > 0) return;

    face_connectivity_.build(*this);
  }
  /**
   * @brief Get the global face index of a local face of an element.
   * @param elem Element index.
   * @param local_face Local face of the element.
   * @return Global face index.
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const override {
    return face_connectivity_.getGlobalFace(elem, local_face);
  }

  /**
   * @brief Get the global node index of a node lying on a face.
   * @param face_global Global face index.
   * @param local_dof Node number on the face, in [0, (order+1)^2).
   * @return Global node index.
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalNodeFromFace(ScalarType face_global, int local_dof) const override {
    return face_connectivity_.getGlobalNodeFromFace(face_global, local_dof);
  }

  /**
   * @brief Tell whether a face lies on the domain boundary.
   *
   * If node boundary flags are stored, a face is a boundary face when none of its nodes
   * is InteriorNode; otherwise the face connectivity decides (no neighbor element).
   * @param face_global Global face index.
   * @return True if boundary face.
   */
  PROXY_HOST_DEVICE
  bool isBoundaryFace(ScalarType face_global) const override {
    if (boundaries_t_.extent(0) == 0) return face_connectivity_.isBoundaryFace(face_global);
    int const n_dofs = face_connectivity_.getDofsPerFace();
    for (int q = 0; q < n_dofs; ++q) {
      if (boundaries_t_[getGlobalNodeFromFace(face_global, q)] == static_cast<ScalarType>(BoundaryFlag::InteriorNode))
        return false;
    }
    return true;
  }

  /**
   * @brief Get the number of unique faces.
   * @return Number of faces (0 until buildFaceConnectivity() has run).
   */
  PROXY_HOST_DEVICE
  ScalarType getNumberOfFaces() const override { return face_connectivity_.getNumberOfFaces(); }

  /**
   * @brief Tell whether a node is on the free surface.
   * @param n Node index.
   * @return True if the node flag is Surface; false if no boundary flags are stored.
   */
  PROXY_HOST_DEVICE
  bool isFreeSurface(ScalarType n) const override {
    if (boundaries_t_.extent(0) == 0) return false;
    return boundaries_t_[n] == static_cast<ScalarType>(BoundaryFlag::Surface);
  }

  /**
   * @brief Set uniform per-element quality factors.
   *
   * Allocates and fills the per-element Qp and Qs arrays; host only.
   * @param qp P-wave quality factor (dimensionless).
   * @param qs S-wave quality factor (dimensionless).
   */
  void setQualityFactors(FloatType qp, FloatType qs) override {
    ScalarType nElem = getNumberOfElements();
    model_qp_element_ = allocateVector<vectorReal>(nElem, "model_qp_element");
    model_qs_element_ = allocateVector<vectorReal>(nElem, "model_qs_element");
    for (ScalarType e = 0; e < nElem; ++e) {
      model_qp_element_(e) = qp;
      model_qs_element_(e) = qs;
    }
  }

 private:
  ScalarType order_;
  ScalarType n_element_;
  ScalarType n_node_;
  FloatType lx_, ly_, lz_;
  FloatType ox_, oy_, oz_;

  int n_points_per_element_;
  bool isModelOnNodes_;
  bool isElastic_;

  arrayInt global_node_index_;
  vectorReal nodes_coords_x_;
  vectorReal nodes_coords_y_;
  vectorReal nodes_coords_z_;

  vectorReal model_vp_node_;
  vectorReal model_vp_element_;
  vectorReal model_rho_node_;
  vectorReal model_rho_element_;
  vectorReal model_vs_node_;
  vectorReal model_vs_element_;
  vectorReal model_qp_node_;
  vectorReal model_qp_element_;
  vectorReal model_qs_node_;
  vectorReal model_qs_element_;
  vectorReal model_delta_node_;
  vectorReal model_delta_element_;
  vectorReal model_epsilon_node_;
  vectorReal model_epsilon_element_;
  vectorReal model_theta_node_;
  vectorReal model_theta_element_;
  vectorReal model_gamma_node_;
  vectorReal model_gamma_element_;
  vectorReal model_phi_node_;
  vectorReal model_phi_element_;
  array3DReal model_C_tensor_element_;
  vectorInt boundaries_t_;

  FaceConnectivityUnstruct<FloatType, ScalarType> face_connectivity_;
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_MODEL_UNSTRUCT_INCLUDE_MODEL_UNSTRUCT_H_
