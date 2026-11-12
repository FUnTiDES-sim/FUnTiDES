// Tests of the 2-D acoustic SEM solver (COMPILE_SEM_2D).

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "cartesian_unstruct_2d_builder.h"
#include "rhs_acoustic.h"
#include "sem_solver_2d.h"
#include "solver_factory.h"
#include "wavefield_acoustic.h"

namespace solver {
namespace test {

constexpr float kVp = 1500.f;  // default vp of CartesianUnstruct2DBuilder, rho = 1

// ---------------------------------------------------------------------------
// Meshes
// ---------------------------------------------------------------------------

/// Cartesian 2-D model from the builder (boundary flags: Damping, Surface on top if asked).
std::shared_ptr<model::ModelApi<float, int>> builderModel(int order, int ex, int ez, float lx, float lz,
                                                          bool freeSurface) {
  model::CartesianParams<float, int> p(order, ex, 1, ez, lx, 1.f, lz, false, false);
  return model::CartesianUnstruct2DBuilder<float, int>(p).getModel(freeSurface);
}

/// Same grid without boundary flags: natural (rigid wall) condition on every edge.
std::shared_ptr<model::ModelUnstruct<float, int>> rigidModel(int order, int ex, int ez, float lx, float lz) {
  int const n1d = order + 1, nx = ex * order + 1, nz = ez * order + 1;
  model::ModelUnstructData<float, int> d;
  d.dim_ = 2;
  d.order_ = order;
  d.n_element_ = ex * ez;
  d.n_node_ = nx * nz;
  d.lx_ = lx;
  d.ly_ = 0;
  d.lz_ = lz;
  d.ox_ = d.oy_ = d.oz_ = 0;
  d.isModelOnNodes_ = false;
  d.isElastic_ = false;
  d.global_node_index_ = allocateArray2D<dualArrayInt>(ex * ez, n1d * n1d, "gni");
  d.nodes_coords_x_ = allocateVector<dualVectorReal>(nx * nz, "x");
  d.nodes_coords_y_ = allocateVector<dualVectorReal>(nx * nz, "y");
  d.nodes_coords_z_ = allocateVector<dualVectorReal>(nx * nz, "z");
  d.model_vp_element_ = allocateVector<dualVectorReal>(ex * ez, "vp");
  d.model_rho_element_ = allocateVector<dualVectorReal>(ex * ez, "rho");
  auto gni = d.global_node_index_.view_host();
  for (int k = 0; k < ez; ++k)
    for (int i = 0; i < ex; ++i)
      for (int m = 0; m < n1d; ++m)
        for (int l = 0; l < n1d; ++l) gni(i + k * ex, l + m * n1d) = (i * order + l) + (k * order + m) * nx;
  float const hx = lx / ex, hz = lz / ez;
  auto coord = [&](int g, float h) { return (g / order) * h + (GLLPoints::get(order, g % order) + 1) * h / 2; };
  for (int gk = 0; gk < nz; ++gk)
    for (int gi = 0; gi < nx; ++gi) {
      d.nodes_coords_x_.view_host()(gi + gk * nx) = coord(gi, hx);
      d.nodes_coords_z_.view_host()(gi + gk * nx) = coord(gk, hz);
    }
  Kokkos::deep_copy(d.model_vp_element_.view_host(), kVp);
  Kokkos::deep_copy(d.model_rho_element_.view_host(), 1.f);
  for (auto* v : {&d.nodes_coords_x_, &d.nodes_coords_y_, &d.nodes_coords_z_, &d.model_vp_element_,
                  &d.model_rho_element_})
    v->modify_host();
  d.global_node_index_.modify_host();
  return std::make_shared<model::ModelUnstruct<float, int>>(d);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

template <typename T>
std::vector<float> toHost(T const& view) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
  return std::vector<float>(h.data(), h.data() + h.extent(0));
}

inline void fromHost(vectorReal& view, std::vector<float> const& values) {
  auto h = Kokkos::create_mirror_view(view);
  for (size_t i = 0; i < values.size(); ++i) h(i) = values[i];
  Kokkos::deep_copy(view, h);
}

/// Wavefield, sources and data of one run. Sources: one element, uniform weights, Ricker 2.
struct Run {
  vectorReal prev, cur;
  std::shared_ptr<SEMsolverDataAcoustic> data;

  Run(int n_node, int n_points, int n_steps, int source_element = -1, float dt = 0, float f0 = 5) {
    prev = allocateVector<vectorReal>(n_node, "prev");
    cur = allocateVector<vectorReal>(n_node, "cur");
    int const n_src = source_element < 0 ? 0 : 1;
    arrayReal term = allocateArray2D<arrayReal>(n_src, n_steps, "term");
    vectorInt elem = allocateVector<vectorInt>(n_src, "elem");
    arrayReal weights = allocateArray2D<arrayReal>(n_src, n_points, "weights");
    if (n_src) {
      auto t = Kokkos::create_mirror_view(term);
      auto e = Kokkos::create_mirror_view(elem);
      auto w = Kokkos::create_mirror_view(weights);
      e(0) = source_element;
      for (int q = 0; q < n_points; ++q) w(0, q) = 1.f / n_points;
      double const a = M_PI * M_PI * f0 * f0;
      for (int s = 0; s < n_steps; ++s) {
        double const tt = s * dt - 1.0 / f0;
        t(0, s) = static_cast<float>((1 - 2 * a * tt * tt) * std::exp(-a * tt * tt));
      }
      Kokkos::deep_copy(term, t);
      Kokkos::deep_copy(elem, e);
      Kokkos::deep_copy(weights, w);
    }
    data = std::make_shared<SEMsolverDataAcoustic>(WavefieldAcoustic(prev, cur), RhsAcoustic(term, elem, weights));
  }

  template <typename SOLVER>
  void steps(SOLVER& solver, float dt, int first, int count) {
    for (int s = first; s < first + count; ++s) {
      solver.computeOneStep(dt, s, *data);
      data->swapWavefields();
    }
  }

  std::vector<float> current() const { return toHost(data->getCurrentField(0)); }
};

// ---------------------------------------------------------------------------
// Typed tests over the compiled orders
// ---------------------------------------------------------------------------

template <typename T>
class SemSolver2DTest : public ::testing::Test {
 public:
  static constexpr int kOrder = T::value;
  using Solver2D = SEMsolver2D<kOrder, false>;
  using Element = typename Solver2D::Element;
};

using Orders = ::testing::Types<std::integral_constant<int, 1>
#if MAX_SOLVER_ACOUSTIC_ORDER >= 2
                                ,
                                std::integral_constant<int, 2>
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 3
                                ,
                                std::integral_constant<int, 3>
#endif
#if MAX_SOLVER_ACOUSTIC_ORDER >= 4
                                ,
                                std::integral_constant<int, 4>
#endif
                                >;
TYPED_TEST_SUITE(SemSolver2DTest, Orders);

TYPED_TEST(SemSolver2DTest, MassSumsToAreaOverVp2Rho) {
  constexpr int O = TestFixture::kOrder;
  auto model = builderModel(O, 5, 3, 1000.f, 600.f, true);
  typename TestFixture::Solver2D solver;
  solver.computeFEInit(*model, {0, 0, 0}, false, 0);
  double total = 0;
  for (float m : toHost(solver.getMassMatrixAcoustic())) {
    EXPECT_GT(m, 0.f);
    total += m;
  }
  double const expected = 1000.0 * 600.0 / (kVp * kVp);
  EXPECT_NEAR(total, expected, 1e-4 * expected);
}

TYPED_TEST(SemSolver2DTest, BoundaryEdgesAndDampingLength) {
  constexpr int O = TestFixture::kOrder;
  int const ex = 5, ez = 3;
  float const lx = 1000.f, lz = 600.f, hz = lz / ez;
  for (bool freeSurface : {false, true}) {
    auto model = builderModel(O, ex, ez, lx, lz, freeSurface);
    typename TestFixture::Solver2D solver;
    solver.computeFEInit(*model, {0, 0, 0}, false, 0);
    EXPECT_EQ(static_cast<int>(solver.getBoundaryEdges().extent(0)), 2 * (ex + ez));
    double total = 0;
    for (float d : toHost(solver.getDampingMatrix(0))) total += d;
    // Absorbing length / (rho vp). With a free surface the top edge is not absorbing, nor the
    // top node of the two side edges (weight w0 * hz / 2 each).
    double const w0 = TestFixture::Element::quadratureWeight(0);
    double const length = freeSurface ? lx + 2 * lz - w0 * hz : 2 * (lx + lz);
    EXPECT_NEAR(total, length / kVp, 1e-4 * length / kVp) << "free surface " << freeSurface;
  }
}

TYPED_TEST(SemSolver2DTest, ConstantFieldHasNoStiffnessForce) {
  constexpr int O = TestFixture::kOrder;
  auto model = rigidModel(O, 4, 4, 800.f, 800.f);
  typename TestFixture::Solver2D solver;
  solver.computeFEInit(*model, {0, 0, 0}, false, 0);
  int const n = model->getNumberOfNodes();
  Run run(n, (O + 1) * (O + 1), 1);
  // Reference scale: forces of a non-constant field
  std::vector<float> f(n);
  for (int i = 0; i < n; ++i) f[i] = model->nodeCoord(i, 0) * model->nodeCoord(i, 2);
  fromHost(run.cur, f);
  solver.computeForces(1e-3f, 0, *run.data);
  double scale = 0;
  for (float w : toHost(solver.getForceVector(0))) scale = std::max(scale, (double)std::abs(w));
  ASSERT_GT(scale, 0);
  // Constant field: zero forces
  std::fill(f.begin(), f.end(), 3.f);
  fromHost(run.cur, f);
  solver.computeForces(1e-3f, 0, *run.data);
  for (float w : toHost(solver.getForceVector(0))) EXPECT_NEAR(w, 0.f, 1e-4 * scale * 3);
}

TYPED_TEST(SemSolver2DTest, InterfaceOperatorAddsSegmentProducts) {
  // Two segments with arbitrary (structured) data: the force vector gains, for each segment,
  // sum_q w_q ((gamma [u]_q - {u}_q) [v]_q - [u]_q {v}_q) (see Solver::setInterfaceOperator()).
  constexpr int O = TestFixture::kOrder;
  constexpr int n1 = O + 1;
  auto model = rigidModel(O, 3, 3, 300.f, 300.f);
  typename TestFixture::Solver2D solver;
  solver.computeFEInit(*model, {0, 0, 0}, false, 0);
  int const n = model->getNumberOfNodes();
  Run run(n, n1 * n1, 1);
  std::vector<float> u(n);
  for (int i = 0; i < n; ++i) u[i] = std::sin(0.7f * i) + 0.1f * i;
  fromHost(run.cur, u);
  solver.computeForces(1e-3f, 0, *run.data);
  std::vector<float> const base = toHost(solver.getForceVector(0));

  int const n_seg = 2, nq = 3;
  std::vector<int> face(n_seg * 2 * n1), line(n_seg * 2 * n1 * n1);
  std::vector<float> dn(n_seg * 2 * n1), lag(n_seg * 2 * nq * n1), hir(n_seg * 2 * nq), w(n_seg * nq),
      gam = {3.f, 0.5f};
  for (size_t i = 0; i < face.size(); ++i) face[i] = static_cast<int>((7 * i + 3) % n);
  for (size_t i = 0; i < line.size(); ++i) line[i] = static_cast<int>((11 * i + 1) % n);
  for (size_t i = 0; i < dn.size(); ++i) dn[i] = 0.3f * (i % 4) - 0.5f;
  for (size_t i = 0; i < lag.size(); ++i) lag[i] = 0.2f * (i % 5) - 0.1f;
  for (size_t i = 0; i < hir.size(); ++i) hir[i] = 0.5f + 0.1f * (i % 3);
  for (size_t i = 0; i < w.size(); ++i) w[i] = 0.7f + 0.2f * (i % 2);
  auto toInt = [](std::vector<int> const& v, const char* name) {
    vectorInt out = allocateVector<vectorInt>(v.size(), name);
    auto h = Kokkos::create_mirror_view(out);
    for (size_t i = 0; i < v.size(); ++i) h(i) = v[i];
    Kokkos::deep_copy(out, h);
    return out;
  };
  auto toReal = [](std::vector<float> const& v, const char* name) {
    vectorReal out = allocateVector<vectorReal>(v.size(), name);
    fromHost(out, v);
    return out;
  };
  vectorInt vface = toInt(face, "face"), vline = toInt(line, "line");
  vectorReal vdn = toReal(dn, "dn"), vlag = toReal(lag, "lag"), vhir = toReal(hir, "hir"), vw = toReal(w, "w"),
             vgam = toReal(gam, "gamma");
  solver.setInterfaceOperator(vface, vline, vdn, vlag, vhir, vw, vgam, nq);
  EXPECT_EQ(solver.getNumberOfInterfaceSegments(), n_seg);
  solver.computeForces(1e-3f, 0, *run.data);
  std::vector<float> const with = toHost(solver.getForceVector(0));

  std::vector<double> expected(n, 0.0);
  for (int s = 0; s < n_seg; ++s) {
    double t[2][n1], d[2][n1];
    for (int c = 0; c < 2; ++c)
      for (int b = 0; b < n1; ++b) {
        int const fb = (s * 2 + c) * n1 + b;
        t[c][b] = u[face[fb]];
        d[c][b] = 0;
        for (int a = 0; a < n1; ++a) d[c][b] += dn[(s * 2 + c) * n1 + a] * u[line[fb * n1 + a]];
      }
    for (int q = 0; q < nq; ++q) {
      double jump = 0, flux = 0;
      for (int c = 0; c < 2; ++c) {
        double tr = 0, dr = 0;
        for (int b = 0; b < n1; ++b) {
          tr += lag[((s * 2 + c) * nq + q) * n1 + b] * t[c][b];
          dr += lag[((s * 2 + c) * nq + q) * n1 + b] * d[c][b];
        }
        jump += c == 0 ? tr : -tr;
        flux += hir[(s * 2 + c) * nq + q] * dr;
      }
      double const c1 = w[s * nq + q] * (gam[s] * jump - flux), c2 = -w[s * nq + q] * jump;
      for (int c = 0; c < 2; ++c)
        for (int b = 0; b < n1; ++b) {
          double const l = lag[((s * 2 + c) * nq + q) * n1 + b];
          int const fb = (s * 2 + c) * n1 + b;
          expected[face[fb]] += l * (c == 0 ? c1 : -c1);
          for (int a = 0; a < n1; ++a)
            expected[line[fb * n1 + a]] += dn[(s * 2 + c) * n1 + a] * hir[(s * 2 + c) * nq + q] * l * c2;
        }
    }
  }
  for (int i = 0; i < n; ++i) EXPECT_NEAR(with[i] - base[i], expected[i], 1e-4 * (1 + std::abs(expected[i])));

  // Wrong sizes or node indices are rejected; an empty operator removes it.
  vectorReal short_lag = toReal(std::vector<float>(lag.size() - 1, 0.f), "short");
  EXPECT_THROW(solver.setInterfaceOperator(vface, vline, vdn, short_lag, vhir, vw, vgam, nq), std::invalid_argument);
  EXPECT_THROW(solver.setInterfaceOperator(vface, vline, vdn, vlag, vhir, vw, vgam, nq + 1), std::invalid_argument);
  std::vector<int> bad_face = face;
  bad_face[0] = n;
  EXPECT_THROW(solver.setInterfaceOperator(toInt(bad_face, "bad"), vline, vdn, vlag, vhir, vw, vgam, nq),
               std::invalid_argument);
  solver.setInterfaceOperator(vectorInt(), vectorInt(), vectorReal(), vectorReal(), vectorReal(), vectorReal(),
                              vectorReal(), nq);
  EXPECT_EQ(solver.getNumberOfInterfaceSegments(), 0);
}

TYPED_TEST(SemSolver2DTest, RigidBoxEigenmodeIsReproduced) {
  // p = cos(pi x / L) cos(pi z / L) cos(w t), w = vp pi sqrt(2) / L, in a rigid square box.
  constexpr int O = TestFixture::kOrder;
  float const L = 1000.f;
  auto model = rigidModel(O, 6, 6, L, L);
  typename TestFixture::Solver2D solver;
  solver.computeFEInit(*model, {0, 0, 0}, false, 0);
  int const n = model->getNumberOfNodes();
  double const w = kVp * M_PI * std::sqrt(2.0) / L;
  float const dt = 0.2f * model->getMinSpacing() / kVp;
  int const n_steps = static_cast<int>(0.6 / dt);
  std::vector<float> mode(n), prev(n);
  for (int i = 0; i < n; ++i) {
    mode[i] = std::cos(M_PI * model->nodeCoord(i, 0) / L) * std::cos(M_PI * model->nodeCoord(i, 2) / L);
    prev[i] = mode[i] * std::cos(w * dt);  // p(-dt)
  }
  Run run(n, (O + 1) * (O + 1), n_steps);
  fromHost(run.cur, mode);
  fromHost(run.prev, prev);
  run.steps(solver, dt, 0, n_steps);
  auto p = run.current();
  double const c = std::cos(w * n_steps * dt);
  double err = 0, ref = 0;
  for (int i = 0; i < n; ++i) {
    err = std::max(err, std::abs(p[i] - c * mode[i]));
    ref = std::max(ref, std::abs(c * mode[i]));
  }
  // 6 x 6 elements, 0.6 s (about 0.6 period)
  std::cout << "[ eigenmode ] order " << O << ": relative max error " << err / ref << std::endl;
  double const tol[] = {0, 1e-1, 1e-2, 1e-3, 1e-3};
  EXPECT_LT(err / ref, tol[O]) << "order " << O;
}

TYPED_TEST(SemSolver2DTest, PointSourceFieldIsSymmetric) {
  constexpr int O = TestFixture::kOrder;
  int const ex = 7, ez = 7;  // odd: the source element is centred
  auto model = builderModel(O, ex, ez, 700.f, 700.f, false);
  typename TestFixture::Solver2D solver;
  solver.computeFEInit(*model, {0, 0, 0}, false, 0);
  int const n = model->getNumberOfNodes(), nx = ex * O + 1;
  float const dt = 0.2f * model->getMinSpacing() / kVp;
  int const n_steps = static_cast<int>(0.3 / dt);
  Run run(n, (O + 1) * (O + 1), n_steps, ex / 2 + (ez / 2) * ex, dt, 10);
  run.steps(solver, dt, 0, n_steps);
  auto p = run.current();
  double pmax = 0, asym = 0;
  for (int i = 0; i < n; ++i) {
    int const gi = i % nx, gk = i / nx;
    pmax = std::max(pmax, (double)std::abs(p[i]));
    asym = std::max(asym, (double)std::abs(p[i] - p[(nx - 1 - gi) + gk * nx]));  // x mirror
    asym = std::max(asym, (double)std::abs(p[i] - p[gk + gi * nx]));             // x <-> z
  }
  ASSERT_GT(pmax, 0);
  EXPECT_LT(asym / pmax, 1e-4);
}

TYPED_TEST(SemSolver2DTest, FreeSurfaceNodesStayAtZero) {
  constexpr int O = TestFixture::kOrder;
  int const ex = 6, ez = 4;
  auto model = builderModel(O, ex, ez, 600.f, 400.f, true);
  typename TestFixture::Solver2D solver;
  solver.computeFEInit(*model, {0, 0, 0}, false, 0);
  int const n = model->getNumberOfNodes();
  float const dt = 0.2f * model->getMinSpacing() / kVp;
  int const n_steps = static_cast<int>(0.3 / dt);
  Run run(n, (O + 1) * (O + 1), n_steps, ex / 2 + (ez - 1) * ex, dt, 10);  // top row
  run.steps(solver, dt, 0, n_steps);
  auto p = run.current();
  double pmax = 0;
  for (int i = 0; i < n; ++i) {
    pmax = std::max(pmax, (double)std::abs(p[i]));
    if (model->isFreeSurface(i)) EXPECT_EQ(p[i], 0.f) << "node " << i;
  }
  EXPECT_GT(pmax, 0);
}

TYPED_TEST(SemSolver2DTest, AbsorbingEdgesRemoveTheEnergy) {
  // A Gaussian pulse leaves the box. Measure: sum of M p^2, averaged over the last quarter of
  // the run, relative to its initial value. It stays of order 1 in the rigid box (the energy
  // only moves between pressure and velocity) and must drop with absorbing edges.
  constexpr int O = TestFixture::kOrder;
  float const L = 800.f;
  int const ne = 16;
  auto absorbing = builderModel(O, ne, ne, L, L, false);
  auto rigid = rigidModel(O, ne, ne, L, L);
  float const dt = 0.2f * rigid->getMinSpacing() / kVp;
  int const n_steps = static_cast<int>(1.5 * L / kVp / dt);
  auto remaining = [&](model::ModelApi<float, int>& m) {
    typename TestFixture::Solver2D solver;
    solver.computeFEInit(m, {0, 0, 0}, false, 0);
    int const n = m.getNumberOfNodes();
    std::vector<float> p0(n);
    for (int i = 0; i < n; ++i) {
      double const x = m.nodeCoord(i, 0) - L / 2, z = m.nodeCoord(i, 2) - L / 2;
      p0[i] = std::exp(-(x * x + z * z) / (2 * 60.0 * 60.0));
    }
    auto mass = toHost(solver.getMassMatrixAcoustic());
    auto mp2 = [&](std::vector<float> const& p) {
      double e = 0;
      for (int i = 0; i < n; ++i) e += mass[i] * p[i] * p[i];
      return e;
    };
    Run run(n, (O + 1) * (O + 1), n_steps);
    fromHost(run.cur, p0);
    fromHost(run.prev, p0);
    int const late = n_steps - n_steps / 4;
    run.steps(solver, dt, 0, late);
    double avg = 0;
    for (int s = late; s < n_steps; ++s) {
      run.steps(solver, dt, s, 1);
      avg += mp2(run.current());
    }
    return avg / (n_steps - late) / mp2(p0);
  };
  double const r_abs = remaining(*absorbing), r_rig = remaining(*rigid);
  std::cout << "[ absorbing ] order " << O << ": remaining " << r_abs << " (absorbing), " << r_rig << " (rigid)"
            << std::endl;
  EXPECT_GT(r_rig, 0.1);
  EXPECT_LT(r_abs, 0.05 * r_rig);
}

TYPED_TEST(SemSolver2DTest, RejectsUnsupportedInput) {
  constexpr int O = TestFixture::kOrder;
  typename TestFixture::Solver2D solver;
  auto wrongOrder = builderModel(O == 1 ? 2 : 1, 2, 2, 100.f, 100.f, false);
  EXPECT_THROW(solver.computeFEInit(*wrongOrder, {0, 0, 0}, false, 0), std::runtime_error);
  auto model = builderModel(O, 2, 2, 100.f, 100.f, false);
  solver.computeFEInit(*model, {0, 0, 0}, false, 0);
  Run run(model->getNumberOfNodes(), (O + 1) * (O + 1), 1);
  EXPECT_THROW(solver.updateSolutionBackward(1e-3f, *run.data), std::runtime_error);
  EXPECT_THROW(solver.setAnisotropyType(model::AnisotropyType::kVTI), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Factory
// ---------------------------------------------------------------------------

namespace e = utils::enums;

TEST(SemSolver2DFactory, CreatesThe2DSolver) {
  for (auto loc : {e::modelLocationType::kOnElements, e::modelLocationType::kOnNodes}) {
    auto s = solver_factory::createSolver(e::methodType::kSem, e::implemType::kMakutu, e::meshType::kUnstruct, loc,
                                          e::physicType::kAcoustic, 2, 0, 2);
    ASSERT_NE(s, nullptr);
    if (loc == e::modelLocationType::kOnNodes)
      EXPECT_NE((dynamic_cast<SEMsolver2D<2, true>*>(s.get())), nullptr);
    else
      EXPECT_NE((dynamic_cast<SEMsolver2D<2, false>*>(s.get())), nullptr);
  }
}

TEST(SemSolver2DFactory, RejectsUnsupportedConfigurations) {
  auto create = [](e::methodType m, e::meshType mesh, e::physicType ph, int order, int dim) {
    return solver_factory::createSolver(m, e::implemType::kMakutu, mesh, e::modelLocationType::kOnElements, ph, order,
                                        0, dim);
  };
  EXPECT_THROW(create(e::methodType::kSem, e::meshType::kStruct, e::physicType::kAcoustic, 2, 2), std::runtime_error);
  EXPECT_THROW(create(e::methodType::kSem, e::meshType::kUnstruct, e::physicType::kElastic, 2, 2), std::runtime_error);
  EXPECT_THROW(create(e::methodType::kDg, e::meshType::kUnstruct, e::physicType::kAcoustic, 2, 2), std::runtime_error);
  EXPECT_THROW(create(e::methodType::kSem, e::meshType::kUnstruct, e::physicType::kAcoustic, 2, 4), std::runtime_error);
  EXPECT_THROW(create(e::methodType::kSem, e::meshType::kUnstruct, e::physicType::kAcoustic, MAX_SOLVER_ACOUSTIC_ORDER + 1, 2),
               std::runtime_error);
}

TEST(SemSolver2DFactory, ThreeDimensionalSolverRejectsA2DMesh) {
  auto model = builderModel(2, 3, 3, 300.f, 300.f, false);
  auto s3d = solver_factory::createSolver(e::methodType::kSem, e::implemType::kMakutu, e::meshType::kUnstruct,
                                          e::modelLocationType::kOnElements, e::physicType::kAcoustic, 2);
  EXPECT_THROW(s3d->computeFEInit(*model, {0, 0, 0}, false, 0), std::runtime_error);
}

}  // namespace test
}  // namespace solver
