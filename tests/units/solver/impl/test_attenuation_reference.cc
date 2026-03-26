/**
 * @file test_attenuation_reference.cc
 * @brief Quantitative reference test for acoustic SLS wave attenuation
 *        at order 2.
 *
 * Validates the SLS attenuation implementation by checking that the
 * energy decay rate scales correctly with the quality factor Q.
 *
 * Physics background
 * ------------------
 * For a viscoelastic medium modeled with Standard Linear Solid (SLS)
 * relaxation mechanisms, the total wavefield energy decays as:
 *
 *     E(t) ~ exp(-omega_eff * t / Q)
 *
 * where omega_eff is the effective angular frequency of the wavefield
 * and Q is the quality factor.  For two simulations with different Q
 * values but the same source and domain:
 *
 *     ln[E_Q1(t) / E_ref(t)]     Q2
 *     ────────────────────── ≈   ────
 *     ln[E_Q2(t) / E_ref(t)]     Q1
 *
 * This ratio is independent of omega_eff and wave speed, making it a
 * robust verification of the Q-factor scaling.  The derivation follows
 * directly from the exponential decay law: since ln(E_Q/E_ref) ~ -1/Q,
 * the ratio of the log-decays for two different Q values equals
 * Q2/Q1 regardless of the (unknown) dominant frequency.
 *
 * Setup
 * -----
 *   - Domain: 1000 m cube, 5x5x5 elements, order 2 (1331 nodes)
 *   - Source: impulse (p=1) at center node
 *   - Material: Vp = 1500 m/s, rho = 1 kg/m^3 (model default)
 *   - SLS: 1 mechanism at f0 = 10 Hz
 *   - Two Q values: 20 (strong attenuation) and 60 (moderate)
 *   - dt = 0.0003 s, 3000 steps (T = 0.9 s)
 *   - No sponge layer, but the solver applies absorbing (Clayton-
 *     Engquist) boundary conditions via the damping matrix, so the
 *     reference energy is NOT exactly conserved.  The energy ratio
 *     E_att/E_ref isolates the attenuation effect because both runs
 *     share the same boundary treatment.
 *
 * Reference
 * ---------
 *   Blanch, Robertsson, Symes (1995). "Modeling of a constant Q:
 *   Methodology and algorithm for an efficient and optimally inexpensive
 *   viscoelastic technique."  Geophysics, 60(1), 176-184.
 *
 *   Carcione J.M. (2007). "Wave Fields in Real Media: Wave Propagation
 *   in Anisotropic, Anelastic, Porous and Electromagnetic Media."
 *   Elsevier, 2nd edition.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "cartesian_struct_builder.h"
#include "common_macros.h"
#include "data_type.h"
#include "rhs_acoustic.h"
#include "sem_solver.h"
#include "sem_solver_data.h"
#include "solver_factory.h"
#include "wavefield_acoustic.h"

namespace solver
{
namespace fe
{
namespace test
{
static VECTOR_REAL_VIEW toView(const std::vector<float>& v, const char* name)
{
  if (v.empty()) return VECTOR_REAL_VIEW();
  auto view = allocateVector<VECTOR_REAL_VIEW>(v.size(), name);
  for (size_t i = 0; i < v.size(); ++i) view[i] = v[i];
  return view;
}

// ======================================================================
// Helper: build a 5x5x5 order-2 mesh with given Q factors
// ======================================================================
static std::shared_ptr<model::ModelApi<float, int>> buildRefMesh(
    float qp = 1e9f, float qs = 1e9f)
{
  constexpr int EX = 5, EY = 5, EZ = 5;
  constexpr float LX = 1000.0f, LY = 1000.0f, LZ = 1000.0f;
  model::CartesianStructBuilder<float, int, 2> builder(EX, LX, EY, LY, EZ, LZ,
                                                       false, false);
  auto mesh = builder.getModel(false);
  mesh->setQualityFactors(qp, qs);
  return mesh;
}

// ======================================================================
// Helper: total wavefield energy  E = sum_i p_i^2
// ======================================================================
static float totalEnergy(const VECTOR_REAL_VIEW& field, int numNodes)
{
  float sum = 0.0f;
  for (int i = 0; i < numNodes; ++i) sum += field(i) * field(i);
  return sum;
}

// ======================================================================
// Helper: run acoustic simulation, recording energy at checkpoints
// ======================================================================
static std::vector<float> runAndRecordEnergy(
    std::shared_ptr<model::ModelApi<float, int>> mesh,
    const std::vector<float>& slsFreqs, int numTimeSteps, float dt,
    const std::vector<int>& checkpoints)
{
  constexpr int order = 2;
  int numNodes = mesh->getNumberOfNodes();
  int npp = (order + 1) * (order + 1) * (order + 1);

  auto solver = solver_factory::createSolver(
      utils::enums::methodType::kSem, utils::enums::implemType::kMakutu,
      utils::enums::meshType::kStruct,
      utils::enums::modelLocationType::kOnElements,
      utils::enums::physicType::kAcoustic, order);
  solver->setAnisotropyType(model::AnisotropyType::kIso);

  if (!slsFreqs.empty()) solver->setSLSAttenuation(toView(slsFreqs, "f"));

  solver->computeFEInit(*mesh, {0, 0, 0}, false, 0.0f);

  auto pPrev = allocateVector<VECTOR_REAL_VIEW>(numNodes, "pP_ref");
  auto pCurr = allocateVector<VECTOR_REAL_VIEW>(numNodes, "pC_ref");
  for (int i = 0; i < numNodes; ++i)
  {
    pPrev(i) = 0.0f;
    pCurr(i) = 0.0f;
  }
  // Impulse at center node
  pCurr(numNodes / 2) = 1.0f;

  auto rhsTerm = allocateArray2D<ARRAY_REAL_VIEW>(1, numTimeSteps, "rt_ref");
  auto rhsElem = allocateVector<VECTOR_INT_VIEW>(1, "re_ref");
  auto rhsWeights = allocateArray2D<ARRAY_REAL_VIEW>(1, npp, "rw_ref");
  rhsElem(0) = 0;
  for (int j = 0; j < npp; ++j)
  {
    rhsTerm(0, j) = 0.0f;
    rhsWeights(0, j) = 0.0f;
  }

  WavefieldAcoustic wf(pPrev, pCurr);
  RhsAcoustic rhs(rhsTerm, rhsElem, rhsWeights);
  SEMsolverDataAcoustic data(wf, rhs);

  std::vector<float> energies;
  int checkIdx = 0;

  for (int t = 0; t < numTimeSteps; ++t)
  {
    solver->computeForces(dt, t, data);
    solver->updateSolution(dt, data);
    data.swapWavefields();

    if (checkIdx < static_cast<int>(checkpoints.size()) &&
        t + 1 == checkpoints[checkIdx])
    {
      energies.push_back(totalEnergy(data.getCurrentField(0), numNodes));
      ++checkIdx;
    }
  }
  return energies;
}

// ======================================================================
// TEST: Energy decay rate scales correctly with Q
//
//   Three runs (same source, same mesh, same dt):
//     - No attenuation (reference)
//     - Q = 20  (strong)
//     - Q = 60  (moderate)
//
//   At the early-time checkpoint (t = 0.3 s, before the wave has
//   fully reflected off all boundaries), verify:
//     1. Q=20 loses more energy than Q=60
//     2. Both lose more energy than the reference
//     3. The log-decay ratio  ln(R20)/ln(R60) ~ Q60/Q20 = 3
//        (within a tolerance to account for the broadband source
//        and absorbing-boundary interplay)
// ======================================================================
TEST(AttenuationReferenceAcoustic, Order2EnergyDecayScalesWithQ)
{
  const int numTimeSteps = 1000;
  const float dt = 0.0003f;

  // Checkpoints at t = 0.1 s, 0.2 s, 0.3 s
  // At t = 0.3 s the wave has just reached the far boundaries
  // (Vp * t = 1500 * 0.3 = 450 m in a 1000 m domain), so the
  // attenuation effect is cleanest before boundary reflections
  // wash out the signal.
  std::vector<int> checkpoints = {333, 667, 1000};

  // 1 SLS mechanism at f0 = 10 Hz
  std::vector<float> freqs = {2.0f * static_cast<float>(M_PI) * 10.0f};

  // ---- Run 1: no attenuation (reference) ----
  auto mesh_ref = buildRefMesh();
  auto E_ref = runAndRecordEnergy(mesh_ref, {}, numTimeSteps, dt, checkpoints);

  // ---- Run 2: Q = 20 ----
  auto mesh_q20 = buildRefMesh(20.0f, 20.0f);
  auto E_q20 =
      runAndRecordEnergy(mesh_q20, freqs, numTimeSteps, dt, checkpoints);

  // ---- Run 3: Q = 60 ----
  auto mesh_q60 = buildRefMesh(60.0f, 60.0f);
  auto E_q60 =
      runAndRecordEnergy(mesh_q60, freqs, numTimeSteps, dt, checkpoints);

  // ==== Assertions ====

  // 1. All energies must be finite and positive
  for (size_t i = 0; i < checkpoints.size(); ++i)
  {
    ASSERT_TRUE(std::isfinite(E_ref[i]))
        << "Reference energy not finite at checkpoint " << i;
    ASSERT_TRUE(std::isfinite(E_q20[i]))
        << "Q=20 energy not finite at checkpoint " << i;
    ASSERT_TRUE(std::isfinite(E_q60[i]))
        << "Q=60 energy not finite at checkpoint " << i;
    ASSERT_GT(E_ref[i], 0.0f);
    ASSERT_GT(E_q20[i], 0.0f);
    ASSERT_GT(E_q60[i], 0.0f);
  }

  // 2. Both attenuated runs should have less energy than reference
  //    at the last checkpoint (earliest checkpoints may not yet show
  //    a clear difference if the wave hasn't propagated far enough).
  EXPECT_LT(E_q20.back(), E_ref.back())
      << "Q=20 should have less energy than reference";
  EXPECT_LT(E_q60.back(), E_ref.back())
      << "Q=60 should have less energy than reference";

  // 3. Q = 20 should lose more energy than Q = 60 at last checkpoint
  EXPECT_LT(E_q20.back(), E_q60.back())
      << "Q=20 should have less energy than Q=60"
      << " (E_Q20=" << E_q20.back() << ", E_Q60=" << E_q60.back() << ")";

  // 4. Quantitative scaling check at the last checkpoint.
  //
  //    Since E_Q(t)/E_ref(t) ~ exp(-omega_eff * t / Q), we have
  //
  //      ln(E_Q20/E_ref) / ln(E_Q60/E_ref) = Q60/Q20 = 3.
  //
  //    We use a tolerance of [2.0, 5.0] because the broadband
  //    impulse excites frequencies where Q_eff != Q_input, and
  //    the absorbing boundary removes some energy non-uniformly.
  float lnR20 = std::log(E_q20.back() / E_ref.back());
  float lnR60 = std::log(E_q60.back() / E_ref.back());

  ASSERT_LT(lnR60, 0.0f)
      << "Q=60 log-ratio should be negative (energy decreasing)";
  ASSERT_LT(lnR20, lnR60)
      << "Q=20 should have a more negative log-ratio than Q=60";

  float scalingRatio = lnR20 / lnR60;

  EXPECT_GT(scalingRatio, 2.0f) << "Q-scaling ratio too low: " << scalingRatio
                                << " (expected ~3.0 = Q60/Q20)";
  EXPECT_LT(scalingRatio, 5.0f) << "Q-scaling ratio too high: " << scalingRatio
                                << " (expected ~3.0 = Q60/Q20)";

  // ---- Diagnostic output ----
  std::cout << "\n=== Acoustic Attenuation Reference Test (Order 2) ===\n"
            << "Domain: 1000 m cube, 5x5x5 elements, order 2\n"
            << "Material: Vp = 1500 m/s, rho = 1 kg/m^3\n"
            << "SLS: 1 mechanism at f0 = 10 Hz\n"
            << "dt = " << dt << " s, steps = " << numTimeSteps
            << ", T = " << dt * numTimeSteps << " s\n\n";
  for (size_t i = 0; i < checkpoints.size(); ++i)
  {
    float t = checkpoints[i] * dt;
    std::cout << "  t = " << t << " s : "
              << "E_ref = " << E_ref[i] << "  E_Q20 = " << E_q20[i]
              << "  E_Q60 = " << E_q60[i] << "  R20 = " << E_q20[i] / E_ref[i]
              << "  R60 = " << E_q60[i] / E_ref[i] << "\n";
  }
  std::cout << "\n  Q-scaling ratio: ln(R20)/ln(R60) = " << scalingRatio
            << " (theoretical: " << 60.0 / 20.0 << ")\n"
            << "=====================================================\n\n";
}

}  // namespace test
}  // namespace fe
}  // namespace solver
