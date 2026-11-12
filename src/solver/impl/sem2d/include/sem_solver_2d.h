#ifndef FUNTIDES_SOLVER_IMPL_SEM2D_INCLUDE_SEM_SOLVER_2D_H_
#define FUNTIDES_SOLVER_IMPL_SEM2D_INCLUDE_SEM_SOLVER_2D_H_
/**
 * @file sem_solver_2d.h
 * @brief 2-D acoustic spectral-element solver on quadrilaterals in the (x, z) plane.
 */

#include <data_type.h>
#include <sem_macros.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Qk_Quadrilateral_Lagrange_GaussLobatto.h"
#include "model.h"
#include "model_unstruct.h"
#include "physics_traits.h"
#include "physics_traits_acoustic.h"
#include "sem_enums.h"
#include "sem_solver_data.h"
#include "solver.h"
#include "sponge_profile.h"

namespace solver {

/**
 * @brief Acoustic SEM solver in 2-D: quadrilateral elements in the (x, z) plane, explicit Verlet
 * time scheme.
 *
 * Same equation and time scheme as the 3-D acoustic SEMsolver, on a ModelUnstruct of
 * dimension 2 (see ModelUnstructData::dim_), with Qk_Quadrilateral_Lagrange_GaussLobatto
 * elements:
 *  - mass: |det J| w / (vp^2 rho), lumped (diagonal);
 *  - stiffness: sum-factorized, coefficient 1 / rho;
 *  - absorbing boundary (first-order paraxial): w * edge length / 2 / (rho vp) on the nodes
 *    flagged Damping of the boundary edges (edges of a single element) that have no InteriorNode;
 *  - free surface: nodes flagged Surface are kept at 0;
 *  - boundary nodes without a flag (InteriorNode) get the natural condition (rigid wall);
 *  - optional sponge layers along x and z (sponge_size[0] and sponge_size[2]): both time levels
 *    multiplied by exp(-sigma dt) at each step, sigma given by sponge::dampingRate();
 *  - optional interface operator (setInterfaceOperator()): the mesh may be made of blocks that
 *    share no node, glued by dense segment matrices (SIPG) added to the stiffness term. Nodes of
 *    such interfaces must not be flagged Damping or Surface.
 * Each step adds the stiffness term and the sources into the force vector, then updates the
 * wavefield: the new value is written into the previous-field buffer, which the caller swaps.
 *
 * Not available in 2-D yet: adjoint update (updateSolutionBackward), attenuation, anisotropy,
 * distributed (MPI) runs.
 *
 * The member functions that launch Kokkos kernels are public: nvcc does not accept extended
 * lambdas in private or protected member functions.
 *
 * @tparam ORDER Polynomial order of the elements.
 * @tparam IS_MODEL_ON_NODES True if vp and rho are stored per node, false if per element.
 */
template <int ORDER, bool IS_MODEL_ON_NODES>
class SEMsolver2D : public Solver {
 public:
  using MeshType = model::ModelUnstruct<float, int>;
  using Element = typename Qk_Quadrilateral_Lagrange_GaussLobatto_Selector<ORDER>::type;
  using DataType = SEMsolverData<utils::enums::physicType::kAcoustic>;

  static constexpr int kN1d = ORDER + 1;                 ///< Nodes along one element axis.
  static constexpr int kPointsPerElement = kN1d * kN1d;  ///< Nodes of an element.
  static constexpr int kMaxInterfacePoints = 16;         ///< Quadrature points per interface segment, at most.

  SEMsolver2D() = default;
  ~SEMsolver2D() = default;

  // ---------------------------------------------------------------------------------------
  // Setup
  // ---------------------------------------------------------------------------------------

  /**
   * @brief Take the mesh, then build the mass and damping matrices and the sponge taper.
   * @param[in] mesh ModelUnstruct<float, int> of dimension 2 and order ORDER.
   * @param[in] sponge_size Sponge widths along x, (unused), z.
   * @param[in] surface_sponge When true, no sponge along the top (z max) edge (free surface).
   * @param[in] taper_delta Attenuation R of the sponge: amplitude ratio of a wave crossing it at
   * normal incidence and coming back, at the largest velocity of the model (see
   * sponge::dampingRate(); 1e-2 if not in (0, 1)).
   * @throws std::runtime_error If the mesh is not a 2-D ModelUnstruct of order ORDER.
   */
  void computeFEInit(model::ModelApi<float, int>& mesh, const std::array<float, 3>& sponge_size,
                     const bool surface_sponge, const float taper_delta) override {
    auto* typed = dynamic_cast<MeshType*>(&mesh);
    if (typed == nullptr) throw std::runtime_error("SEMsolver2D: the mesh must be a ModelUnstruct<float, int>");
    if (typed->getDimension() != 2) throw std::runtime_error("SEMsolver2D: the mesh must be 2-D");
    if (typed->getOrder() != ORDER) throw std::runtime_error("SEMsolver2D: mesh order differs from the solver order");
    m_mesh = *typed;
    sponge_size_ = sponge_size;
    surface_sponge_ = surface_sponge;
    taper_delta_ = taper_delta;

    allocateFEarrays();
    findBoundaryEdges();
    initFEarrays();
    computeGlobalMassMatrix();
    computeDampingMatrix();
    FENCE
  }

  /// Allocates the mass, damping, force and sponge vectors (zero-filled).
  void allocateFEarrays() override {
    int const n = m_mesh.getNumberOfNodes();
    massMatrixGlobal_ = allocateVector<vectorReal>(n, "massMatrix2D");
    dampingMatrixGlobal_ = allocateVector<vectorReal>(n, "dampingMatrix2D");
    workVectorGlobal_ = allocateVector<vectorReal>(n, "workVector2D");
    spongeTaperCoeff_ = allocateVector<dualVectorReal>(n, "spongeTaper2D");
    spongeRate_ = allocateVector<dualVectorReal>(n, "spongeRate2D");
  }

  /// Fills the sponge damping rates.
  void initFEarrays() override { initSpongeValues(); }

  /**
   * @brief Sponge damping rate per node (1/s), from sponge::dampingRate() along x and z (the
   * rates add up on the corners); d is the distance to the edges of the domain bounding box and
   * vmax the largest velocity of the model. The taper exp(-rate dt) is computed at the first
   * update with a new dt.
   */
  void initSpongeValues() override {
    int const n = m_mesh.getNumberOfNodes();
    auto rate = spongeRate_.view_host();
    double xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
    double zmin = xmin, zmax = xmax;
    for (int i = 0; i < n; ++i) {
      xmin = std::min<double>(xmin, m_mesh.nodeCoord(i, 0));
      xmax = std::max<double>(xmax, m_mesh.nodeCoord(i, 0));
      zmin = std::min<double>(zmin, m_mesh.nodeCoord(i, 2));
      zmax = std::max<double>(zmax, m_mesh.nodeCoord(i, 2));
    }
    double vmax = 0.0;
    if constexpr (IS_MODEL_ON_NODES) {
      for (int i = 0; i < n; ++i) vmax = std::max<double>(vmax, m_mesh.getModelVpOnNodes(i));
    } else {
      for (int e = 0; e < m_mesh.getNumberOfElements(); ++e)
        vmax = std::max<double>(vmax, m_mesh.getModelVpOnElement(e));
    }
    for (int i = 0; i < n; ++i) {
      double const x = m_mesh.nodeCoord(i, 0), z = m_mesh.nodeCoord(i, 2);
      double const dx = std::min(x - xmin, xmax - x);
      double const dz = surface_sponge_ ? (z - zmin) : std::min(z - zmin, zmax - z);
      rate(i) = static_cast<float>(sponge::dampingRate(dx, sponge_size_[0], vmax, taper_delta_) +
                                   sponge::dampingRate(dz, sponge_size_[2], vmax, taper_delta_));
    }
    spongeRate_.modify_host();
    spongeRate_.sync_device();
    Kokkos::deep_copy(spongeTaperCoeff_.view_host(), 1.0f);  // 1 until the first update with a dt
    spongeTaperCoeff_.modify_host();
    spongeTaperCoeff_.sync_device();
    spongeTaperDt_ = -1.0f;
  }

  /// Sets the taper to exp(-rate dt) if dt differs from the one of the last call.
  void updateSpongeTaper(float dt);

  /**
   * @brief Lists the boundary edges (edges of a single element) as (element, local edge) pairs.
   *
   * Host only. An edge is identified by its two vertices; local edges are numbered as in
   * Qk_Quadrilateral_Lagrange_GaussLobatto::edgeNodeIndex().
   */
  void findBoundaryEdges() {
    int const n_elem = m_mesh.getNumberOfElements();
    long long const n_node = m_mesh.getNumberOfNodes();
    std::unordered_map<long long, std::pair<int, int>> owner;  // edge key -> (count, element*4+edge)
    owner.reserve(static_cast<size_t>(n_elem) * 4);
    for (int e = 0; e < n_elem; ++e)
      for (int ed = 0; ed < 4; ++ed) {
        long long a = edgeVertex(e, ed, 0), b = edgeVertex(e, ed, 1);
        if (a > b) std::swap(a, b);
        auto& entry = owner[a * n_node + b];
        if (entry.first++ == 0) entry.second = e * 4 + ed;
      }
    std::vector<int> edges;
    for (auto const& kv : owner)
      if (kv.second.first == 1) edges.push_back(kv.second.second);
    std::sort(edges.begin(), edges.end());

    boundaryEdges_ = allocateVector<vectorInt>(static_cast<int>(edges.size()), "boundaryEdges2D");
    auto host = Kokkos::create_mirror_view(boundaryEdges_);
    for (size_t i = 0; i < edges.size(); ++i) host(i) = edges[i];
    Kokkos::deep_copy(boundaryEdges_, host);
  }

  /// Global node of vertex v (0: first, 1: last) of local edge ed of element e. Host only.
  int edgeVertex(int e, int ed, int v) const {
    int const l = Element::edgeNodeIndex(ed, v == 0 ? 0 : ORDER);
    return m_mesh.globalNodeIndex(e, l % kN1d, 0, l / kN1d);
  }

  /// Assembles the lumped mass matrix: |det J| w / (vp^2 rho).
  void computeGlobalMassMatrix() override;

  /// Assembles the absorbing-boundary term on the Damping nodes of the boundary edges.
  void computeDampingMatrix() override;

  // ---------------------------------------------------------------------------------------
  // Time stepping
  // ---------------------------------------------------------------------------------------

  /// Zeroes the force vector.
  void resetGlobalVectors(int) override { Kokkos::deep_copy(workVectorGlobal_, 0.0f); }

  /**
   * @brief Force vector of the time sample: stiffness term minus the sources.
   * @throws std::bad_cast If @p data is not an acoustic SEMsolverData.
   */
  void computeForces(const float& dt, const int& timeSample, DataStruct& data) override {
    auto& myData = dynamic_cast<DataType&>(data);
    resetGlobalVectors(m_mesh.getNumberOfNodes());
    FENCE
    applyRHSTerm(timeSample, dt, myData);
    FENCE
    computeElementContributions(myData);
    FENCE
    if (nInterfaceSegments_ > 0) {
      applyInterfaceOperator(myData);
      FENCE
    }
  }

  /**
   * @brief Store the interface operator (see Solver::setInterfaceOperator()).
   * @throws std::invalid_argument If the sizes do not match or a node index is out of range.
   */
  void setInterfaceOperator(const vectorInt& face, const vectorInt& line, const vectorReal& dn, const vectorReal& lag,
                            const vectorReal& half_inv_rho, const vectorReal& weight, const vectorReal& gamma,
                            int points_per_segment) override {
    long long const nf = face.extent(0);
    if (nf == 0) {
      nInterfaceSegments_ = 0;
      return;
    }
    int const nq = points_per_segment;
    if (nq <= 0 || nq > kMaxInterfacePoints || nf % (2 * kN1d) != 0)
      throw std::invalid_argument("SEMsolver2D::setInterfaceOperator: face must hold n_seg * 2 * (order + 1) nodes");
    long long const n_seg = nf / (2 * kN1d);
    if (static_cast<long long>(line.extent(0)) != n_seg * 2 * kN1d * kN1d ||
        static_cast<long long>(dn.extent(0)) != n_seg * 2 * kN1d ||
        static_cast<long long>(lag.extent(0)) != n_seg * 2 * nq * kN1d ||
        static_cast<long long>(half_inv_rho.extent(0)) != n_seg * 2 * nq ||
        static_cast<long long>(weight.extent(0)) != n_seg * nq || static_cast<long long>(gamma.extent(0)) != n_seg)
      throw std::invalid_argument("SEMsolver2D::setInterfaceOperator: array sizes do not match n_seg and n_q");
    int bad = 0;
    int const n_node = m_mesh.getNumberOfNodes();
    Kokkos::parallel_reduce(
        "SEM2D interface nodes check", Kokkos::RangePolicy<>(0, static_cast<int>(line.extent(0))),
        KOKKOS_LAMBDA(const int i, int& b) { b += (line(i) < 0 || line(i) >= n_node) ? 1 : 0; }, bad);
    Kokkos::parallel_reduce(
        "SEM2D interface faces check", Kokkos::RangePolicy<>(0, static_cast<int>(nf)),
        KOKKOS_LAMBDA(const int i, int& b) { b += (face(i) < 0 || face(i) >= n_node) ? 1 : 0; }, bad);
    if (bad > 0) throw std::invalid_argument("SEMsolver2D::setInterfaceOperator: node index out of range");
    ifFace_ = face;
    ifLine_ = line;
    ifDn_ = dn;
    ifLag_ = lag;
    ifHalfInvRho_ = half_inv_rho;
    ifWeight_ = weight;
    ifGamma_ = gamma;
    interfacePointsPerSegment_ = nq;
    nInterfaceSegments_ = static_cast<int>(n_seg);
  }

  /// Adds the interface operator applied to the current field to the force vector.
  void applyInterfaceOperator(const DataType& data);

  /// @return Number of interface segments (0: no interface operator).
  int getNumberOfInterfaceSegments() const { return nInterfaceSegments_; }

  /// Subtracts the weighted source terms of the time sample from the force vector.
  void applyRHSTerm(int timeSample, float dt, const DataType& data);

  /// Adds the stiffness operator applied to the current field to the force vector.
  void computeElementContributions(const DataType& data);

  /**
   * @brief Verlet update of every node from the force vector:
   * p_new = (2 M p - (M - dt D / 2) p_prev - dt^2 F) / (M + dt D / 2) * taper, written into the
   * previous-field buffer; the current field is multiplied by the taper; free-surface nodes are
   * set to 0.
   * @throws std::runtime_error If the wavefield has a third (adjoint) buffer.
   */
  void updateSolutionForward(const float& dt, DataStruct& data) override;

  /// computeForces() then updateSolutionForward(). Distributed runs are not supported.
  void computeOneStep(const float& dt, const int& timeSample, DataStruct& data) override {
    computeForces(dt, timeSample, data);
    updateSolutionForward(dt, data);
  }

  /// @throws std::runtime_error Always: the adjoint update is not available in 2-D yet.
  void updateSolutionBackward(const float&, DataStruct&) override {
    throw std::runtime_error("SEMsolver2D: updateSolutionBackward (adjoint) is not available in 2-D yet");
  }

  // ---------------------------------------------------------------------------------------
  // Accessors and unsupported options
  // ---------------------------------------------------------------------------------------

  int getNumComponents() const override { return 1; }
  vectorReal& getMassMatrixAcoustic() override { return massMatrixGlobal_; }
  vectorReal& getMassMatrixElastic() override { return massMatrixGlobal_; }
  vectorReal& getDampingMatrix(int) override { return dampingMatrixGlobal_; }
  vectorReal& getForceVector(int) override { return workVectorGlobal_; }
  /// @return Sponge taper exp(-rate dt) per node, for the dt of the last update (device).
  vectorReal getSpongeTaperCoeff() const { return spongeTaperCoeff_.view_device(); }
  /// @return Boundary edges, as element * 4 + local edge.
  vectorInt getBoundaryEdges() const { return boundaryEdges_; }

  /// Isotropic only: other anisotropy types are rejected.
  void setAnisotropyType(model::AnisotropyType type) override {
    if (type != model::AnisotropyType::kIso)
      throw std::invalid_argument("SEMsolver2D: only isotropic acoustic media are supported");
  }

  /// @throws std::runtime_error Always: no attenuation in 2-D yet.
  void setSLSAttenuation(const vectorReal&, const vectorReal&) override {
    throw std::runtime_error("SEMsolver2D: attenuation is not available in 2-D yet");
  }

  /// Prints the field value at the first node of element @p e.
  void outputSolutionValues(const int& t, int& e, const vectorReal& field, const char* fieldName) override {
    float value = 0.0f;
    Kokkos::deep_copy(value, Kokkos::subview(field, m_mesh.globalNodeIndex(e, 0, 0, 0)));
    std::cout << "TimeStep=" << t << ";  " << fieldName << " @ element " << e << " = " << value << std::endl;
  }

  /// Prints the first column of @p field at the first node of element @p e.
  void outputSolutionValues(const int& t, int& e, const arrayReal& field, const char* fieldName) override {
    float value = 0.0f;
    Kokkos::deep_copy(value, Kokkos::subview(field, m_mesh.globalNodeIndex(e, 0, 0, 0), 0));
    std::cout << "TimeStep=" << t << ";  " << fieldName << " @ element " << e << " = " << value << std::endl;
  }

 private:
  MeshType m_mesh;
  vectorReal massMatrixGlobal_;
  vectorReal dampingMatrixGlobal_;
  vectorReal workVectorGlobal_;
  dualVectorReal spongeTaperCoeff_;  ///< exp(-spongeRate_ dt), for dt = spongeTaperDt_
  dualVectorReal spongeRate_;        ///< sponge damping rate per node (1/s)
  float spongeTaperDt_{-1.0f};
  vectorInt boundaryEdges_;  ///< element * 4 + local edge, for each boundary edge
  std::array<float, 3> sponge_size_{0, 0, 0};
  bool surface_sponge_{false};
  float taper_delta_{0};
  // Interface operator (see Solver::setInterfaceOperator()), per segment, side, point, node
  vectorInt ifFace_;           ///< face nodes (n_seg * 2 * n1)
  vectorInt ifLine_;           ///< normal lines through the face nodes (n_seg * 2 * n1 * n1)
  vectorReal ifDn_;            ///< normal derivative coefficients on a line (n_seg * 2 * n1)
  vectorReal ifLag_;           ///< Lagrange basis of the face nodes at the points (n_seg * 2 * n_q * n1)
  vectorReal ifHalfInvRho_;    ///< (1 / rho) / 2 at the points (n_seg * 2 * n_q)
  vectorReal ifWeight_;        ///< Gauss weight x length / 2 (n_seg * n_q)
  vectorReal ifGamma_;         ///< penalty (n_seg)
  int interfacePointsPerSegment_{0}; ///< n_q
  int nInterfaceSegments_{0};
};

// Member functions launching Kokkos kernels are defined out of the class: they are then not
// inline, so that the extern template declarations below keep them out of other files.

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2D<ORDER, IS_MODEL_ON_NODES>::computeGlobalMassMatrix() {
  auto mesh = m_mesh;
  auto mass = massMatrixGlobal_;
  Kokkos::parallel_for(
      "SEM2D mass matrix", Kokkos::RangePolicy<>(0, mesh.getNumberOfElements()), KOKKOS_LAMBDA(const int e) {
        real_t X[4][2];
        for (int v = 0; v < 4; ++v) {
          int const node = mesh.globalNodeIndex(e, (v % 2) * ORDER, 0, (v / 2) * ORDER);
          X[v][0] = mesh.nodeCoord(node, 0);
          X[v][1] = mesh.nodeCoord(node, 2);
        }
        Element::computeMassTerm(X, [&](int const q, real_t const val) {
          int const node = mesh.globalNodeIndex(e, q % kN1d, 0, q / kN1d);
          real_t vp, rho;
          if constexpr (IS_MODEL_ON_NODES) {
            vp = mesh.getModelVpOnNodes(node);
            rho = mesh.getModelRhoOnNodes(node);
          } else {
            vp = mesh.getModelVpOnElement(e);
            rho = mesh.getModelRhoOnElement(e);
          }
          ATOMICADD(mass(node), val / (vp * vp * rho));
        });
      });
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2D<ORDER, IS_MODEL_ON_NODES>::computeDampingMatrix() {
  auto mesh = m_mesh;
  auto damping = dampingMatrixGlobal_;
  auto edges = boundaryEdges_;
  Kokkos::parallel_for(
      "SEM2D damping matrix", Kokkos::RangePolicy<>(0, edges.extent(0)), KOKKOS_LAMBDA(const int b) {
        int const e = edges(b) / 4;
        int const ed = edges(b) % 4;
        real_t E[2][2];
        for (int v = 0; v < 2; ++v) {
          int const l = Element::edgeNodeIndex(ed, v == 0 ? 0 : ORDER);
          int const node = mesh.globalNodeIndex(e, l % kN1d, 0, l / kN1d);
          E[v][0] = mesh.nodeCoord(node, 0);
          E[v][1] = mesh.nodeCoord(node, 2);
        }
        // An edge of a single element absorbs only if it lies on the outer boundary: none of its
        // nodes is an interior node. Edges of non-conforming block interfaces (glued by the
        // interface operator) have interior nodes and are skipped, even when one of their ends
        // is a flagged node of the outer boundary.
        for (int i = 0; i < kN1d; ++i) {
          int const l = Element::edgeNodeIndex(ed, i);
          if (mesh.boundaryType(mesh.globalNodeIndex(e, l % kN1d, 0, l / kN1d)) == model::BoundaryFlag::InteriorNode)
            return;
        }
        for (int i = 0; i < kN1d; ++i) {
          int const l = Element::edgeNodeIndex(ed, i);
          int const node = mesh.globalNodeIndex(e, l % kN1d, 0, l / kN1d);
          if (mesh.boundaryType(node) != model::BoundaryFlag::Damping) continue;
          real_t vp, rho;
          if constexpr (IS_MODEL_ON_NODES) {
            vp = mesh.getModelVpOnNodes(node);
            rho = mesh.getModelRhoOnNodes(node);
          } else {
            vp = mesh.getModelVpOnElement(e);
            rho = mesh.getModelRhoOnElement(e);
          }
          ATOMICADD(damping(node), Element::computeDampingTerm(i, E) / (rho * vp));
        }
      });
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2D<ORDER, IS_MODEL_ON_NODES>::applyRHSTerm(int timeSample, float, const DataType& data) {
  auto const rhs_element = data.getRhsElement();
  auto const rhs_term = data.getRhsTerm(0);
  auto const rhs_weights = data.getRhsWeights(0);
  auto mesh = m_mesh;
  auto work = workVectorGlobal_;
  int const n = rhs_element.extent(0) * kPointsPerElement;
  Kokkos::parallel_for(
      "SEM2D apply RHS", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const int t) {
        int const i = t / kPointsPerElement;
        int const l = t - i * kPointsPerElement;
        int const node = mesh.globalNodeIndex(rhs_element(i), l % kN1d, 0, l / kN1d);
        ATOMICADD(work(node), -rhs_term(i, timeSample) * rhs_weights(i, l));
      });
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2D<ORDER, IS_MODEL_ON_NODES>::computeElementContributions(const DataType& data) {
  auto mesh = m_mesh;
  auto work = workVectorGlobal_;
  auto field = data.getCurrentField(0);
  Kokkos::parallel_for(
      "SEM2D stiffness", Kokkos::RangePolicy<>(0, mesh.getNumberOfElements()), KOKKOS_LAMBDA(const int e) {
        real_t u[kPointsPerElement];
        real_t v[kPointsPerElement];
        int nodes[kPointsPerElement];
        for (int l = 0; l < kPointsPerElement; ++l) {
          nodes[l] = mesh.globalNodeIndex(e, l % kN1d, 0, l / kN1d);
          u[l] = field(nodes[l]);
          v[l] = 0;
        }
        real_t X[4][2];
        for (int c = 0; c < 4; ++c) {
          int const node = nodes[Element::meshIndexToLinearIndex2D(c)];
          X[c][0] = mesh.nodeCoord(node, 0);
          X[c][1] = mesh.nodeCoord(node, 2);
        }
        real_t inv_rho = 0;
        if constexpr (!IS_MODEL_ON_NODES) inv_rho = 1.0f / mesh.getModelRhoOnElement(e);
        Element::computeStiffnessTermSumFact(X, u, v, [&](int const qa, int const qb) -> real_t {
          if constexpr (IS_MODEL_ON_NODES) {
            return 1.0f / mesh.getModelRhoOnNodes(nodes[qa + qb * kN1d]);
          } else {
            return inv_rho;
          }
        });
        for (int l = 0; l < kPointsPerElement; ++l) ATOMICADD(work(nodes[l]), v[l]);
      });
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2D<ORDER, IS_MODEL_ON_NODES>::applyInterfaceOperator(const DataType& data) {
  auto work = workVectorGlobal_;
  auto u = data.getCurrentField(0);
  auto face = ifFace_;
  auto line = ifLine_;
  auto dn = ifDn_;
  auto lag = ifLag_;
  auto hir = ifHalfInvRho_;
  auto wgt = ifWeight_;
  auto gamma = ifGamma_;
  int const nq = interfacePointsPerSegment_;
  constexpr int n1 = kN1d;
  // One thread per segment. Traces and normal derivatives are computed once per face node (the normal
  // derivative on an element line), then interpolated at the Gauss points: about four times less data
  // than full rows over the element nodes. Two flat kernels with full rows, and a team per segment,
  // were measured 2 to 5 times slower on GPU.
  Kokkos::parallel_for(
      "SEM2D interface operator", Kokkos::RangePolicy<>(0, nInterfaceSegments_), KOKKOS_LAMBDA(const int s) {
        real_t t[2][n1], d[2][n1], sf[2][n1], el[2][n1];
        for (int c = 0; c < 2; ++c)
          for (int b = 0; b < n1; ++b) {
            long long const fb = (static_cast<long long>(s) * 2 + c) * n1 + b;
            t[c][b] = u(face(fb));
            real_t acc = 0;
            for (int a = 0; a < n1; ++a) acc += dn((static_cast<long long>(s) * 2 + c) * n1 + a) * u(line(fb * n1 + a));
            d[c][b] = acc;
            sf[c][b] = 0;
            el[c][b] = 0;
          }
        real_t const g = gamma(s);
        for (int q = 0; q < nq; ++q) {
          real_t jump = 0, flux = 0;
          for (int c = 0; c < 2; ++c) {
            long long const lq = ((static_cast<long long>(s) * 2 + c) * nq + q) * n1;
            real_t tr = 0, dr = 0;
            for (int b = 0; b < n1; ++b) {
              tr += lag(lq + b) * t[c][b];
              dr += lag(lq + b) * d[c][b];
            }
            jump += (c == 0 ? tr : -tr);
            flux += hir((static_cast<long long>(s) * 2 + c) * nq + q) * dr;
          }
          real_t const w = wgt(static_cast<long long>(s) * nq + q);
          real_t const c1 = w * (g * jump - flux);
          real_t const c2 = -w * jump;
          for (int c = 0; c < 2; ++c) {
            long long const lq = ((static_cast<long long>(s) * 2 + c) * nq + q) * n1;
            real_t const h = hir((static_cast<long long>(s) * 2 + c) * nq + q) * c2;
            real_t const sgn = (c == 0 ? c1 : -c1);
            for (int b = 0; b < n1; ++b) {
              sf[c][b] += lag(lq + b) * sgn;
              el[c][b] += lag(lq + b) * h;
            }
          }
        }
        for (int c = 0; c < 2; ++c)
          for (int b = 0; b < n1; ++b) {
            long long const fb = (static_cast<long long>(s) * 2 + c) * n1 + b;
            ATOMICADD(work(face(fb)), sf[c][b]);
            for (int a = 0; a < n1; ++a)
              ATOMICADD(work(line(fb * n1 + a)), dn((static_cast<long long>(s) * 2 + c) * n1 + a) * el[c][b]);
          }
      });
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2D<ORDER, IS_MODEL_ON_NODES>::updateSpongeTaper(float dt) {
  if (dt == spongeTaperDt_) return;
  auto rate = spongeRate_.view_device();
  auto taper = spongeTaperCoeff_.view_device();
  Kokkos::parallel_for(
      "SEM2D sponge taper", Kokkos::RangePolicy<>(0, m_mesh.getNumberOfNodes()),
      KOKKOS_LAMBDA(const int i) { taper(i) = Kokkos::exp(-rate(i) * dt); });
  spongeTaperCoeff_.modify_device();
  FENCE
  spongeTaperDt_ = dt;
}

template <int ORDER, bool IS_MODEL_ON_NODES>
void SEMsolver2D<ORDER, IS_MODEL_ON_NODES>::updateSolutionForward(const float& dt, DataStruct& data) {
  auto& myData = dynamic_cast<DataType&>(data);
  if (myData.getPrevPrevField(0).extent(0) > 0)
    throw std::runtime_error("SEMsolver2D: updateSolutionForward called with a 3-buffer wavefield");
  auto mesh = m_mesh;
  auto mass = massMatrixGlobal_;
  auto damping = dampingMatrixGlobal_;
  auto work = workVectorGlobal_;
  updateSpongeTaper(dt);
  auto taper = spongeTaperCoeff_.view_device();
  auto cur = myData.getCurrentField(0);
  auto prev = myData.getPreviousField(0);
  float const dt2 = dt * dt;
  float const half_dt = 0.5f * dt;
  Kokkos::parallel_for(
      "SEM2D update", Kokkos::RangePolicy<>(0, mesh.getNumberOfNodes()), KOKKOS_LAMBDA(const int i) {
        float const m = mass(i);
        if (m <= 0.0f) return;
        if (mesh.isFreeSurface(i)) {
          cur(i) = 0.0f;
          prev(i) = 0.0f;
          return;
        }
        float const d = damping(i);
        float const next = 2.0f * m * cur(i) - (m - half_dt * d) * prev(i) - dt2 * work(i);
        prev(i) = next / (m + half_dt * d) * taper(i);
        cur(i) = cur(i) * taper(i);
      });
  FENCE
}

// Each SEMsolver2D is instantiated in a single file (generated SEM2DQ<order>_* sources): the
// other files must not compile its Kokkos lambdas again (nvcc __host__ __device__ lambda
// wrappers break when the same lambda is compiled in several files).
#ifndef FUNTIDES_SEM2D_INSTANTIATION
#define FUNTIDES_SEM2D_EXTERN(O)                 \
  extern template class SEMsolver2D<O, false>; \
  extern template class SEMsolver2D<O, true>;
#if MAX_SOLVER_ACOUSTIC_ORDER >= 1
FUNTIDES_SEM2D_EXTERN(1)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 2
FUNTIDES_SEM2D_EXTERN(2)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 3
FUNTIDES_SEM2D_EXTERN(3)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 4
FUNTIDES_SEM2D_EXTERN(4)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 5
FUNTIDES_SEM2D_EXTERN(5)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 6
FUNTIDES_SEM2D_EXTERN(6)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 7
FUNTIDES_SEM2D_EXTERN(7)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 8
FUNTIDES_SEM2D_EXTERN(8)
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 9
FUNTIDES_SEM2D_EXTERN(9)
#endif
#undef FUNTIDES_SEM2D_EXTERN
#endif

}  // namespace solver

#endif  // FUNTIDES_SOLVER_IMPL_SEM2D_INCLUDE_SEM_SOLVER_2D_H_
