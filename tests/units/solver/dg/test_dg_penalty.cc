/**
 * @file test_dg_penalty.cc
 * @brief Unit tests for dg_penalty.h: computeFaceArea, computeHexVolume,
 *        computeSIPGPenalty.
 *
 * All three functions are PROXY_HOST_DEVICE free functions operating on plain
 * float arrays, so no mesh or Kokkos views are needed. Expected values are
 * derived analytically from the formulas documented in dg_penalty.h:
 *
 *   computeHexVolume returns 8 * |det(J_center)|, NOT the geometric volume.
 *   For a cube of side h: returns 8 * h^3.
 *   computeFaceArea computes the area of a bilinear quad by summing triangle
 *   areas from vertex 0.
 *   computeSIPGPenalty = penalty_factor * (ORDER+1)^2 / (vol / area).
 */
#include <gtest/gtest.h>

#include <cmath>

#include "dg_penalty.h"

namespace solver {
namespace fe {
namespace test {

// ============================================================
// computeFaceArea
// ============================================================

TEST(ComputeFaceArea, UnitSquareXYPlane) {
  float X[4][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  EXPECT_NEAR(computeFaceArea(X), 1.0f, 1e-6f);
}

TEST(ComputeFaceArea, UnitSquareXZPlane) {
  float X[4][3] = {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}};
  EXPECT_NEAR(computeFaceArea(X), 1.0f, 1e-6f);
}

TEST(ComputeFaceArea, Rectangle2x3InXYPlane) {
  float X[4][3] = {{0, 0, 0}, {2, 0, 0}, {2, 3, 0}, {0, 3, 0}};
  EXPECT_NEAR(computeFaceArea(X), 6.0f, 1e-5f);
}

TEST(ComputeFaceArea, DegenerateAllSamePoint) {
  float X[4][3] = {{1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {1, 2, 3}};
  EXPECT_NEAR(computeFaceArea(X), 0.0f, 1e-6f);
}

TEST(ComputeFaceArea, PositiveForNonDegenerateQuad) {
  float X[4][3] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}};
  EXPECT_GT(computeFaceArea(X), 0.0f);
}

// ============================================================
// computeHexVolume
// Corner ordering: X[iv + 2*jv + 4*kv] for iv,jv,kv in {0,1}.
// The function returns 8 * |det(J_center)|.
// For a uniform cube of side h: returns 8 * h^3.
// ============================================================

TEST(ComputeHexVolume, UnitCubeReturns8) {
  // X[iv+2*jv+4*kv]
  float X[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  EXPECT_NEAR(computeHexVolume(X), 8.0f, 1e-5f);
}

TEST(ComputeHexVolume, CubeOfSide2Returns64) {
  float X[8][3] = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {2, 2, 0}, {0, 0, 2}, {2, 0, 2}, {0, 2, 2}, {2, 2, 2}};
  EXPECT_NEAR(computeHexVolume(X), 64.0f, 1e-4f);
}

TEST(ComputeHexVolume, ScalesAsCubicPowerOfSide) {
  // For side h the function returns 8*h^3.
  for (float h : {0.5f, 1.0f, 2.0f, 3.0f}) {
    float X[8][3];
    for (int kv = 0; kv < 2; ++kv)
      for (int jv = 0; jv < 2; ++jv)
        for (int iv = 0; iv < 2; ++iv) {
          int idx = iv + 2 * jv + 4 * kv;
          X[idx][0] = iv * h;
          X[idx][1] = jv * h;
          X[idx][2] = kv * h;
        }
    EXPECT_NEAR(computeHexVolume(X), 8.0f * h * h * h, 1e-3f) << "Failed for h=" << h;
  }
}

TEST(ComputeHexVolume, PositiveForValidHex) {
  float X[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  EXPECT_GT(computeHexVolume(X), 0.0f);
}

// ============================================================
// computeSIPGPenalty
// gamma = penalty_factor * (ORDER+1)^2 / h_f
// h_f  = vol / area
// Unit cube: area=1.0, vol=8.0, h_f=8.0.
//   ORDER=1: gamma = 10 * 4 / 8 = 5.0
//   ORDER=2: gamma = 10 * 9 / 8 = 11.25
// ============================================================

namespace {
// Unit square face and unit cube element for reuse.
float kUnitFace[4][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
float kUnitHex[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
}  // namespace

TEST(ComputeSIPGPenalty, Order1UnitCubePenalty10) {
  float gamma = computeSIPGPenalty<1>(kUnitFace, kUnitHex, 10.0f);
  EXPECT_NEAR(gamma, 5.0f, 1e-5f);
}

TEST(ComputeSIPGPenalty, Order2UnitCubePenalty10) {
  float gamma = computeSIPGPenalty<2>(kUnitFace, kUnitHex, 10.0f);
  EXPECT_NEAR(gamma, 11.25f, 1e-5f);
}

TEST(ComputeSIPGPenalty, LinearInPenaltyFactor) {
  float g1 = computeSIPGPenalty<1>(kUnitFace, kUnitHex, 10.0f);
  float g2 = computeSIPGPenalty<1>(kUnitFace, kUnitHex, 20.0f);
  EXPECT_NEAR(g2, 2.0f * g1, 1e-5f);
}

TEST(ComputeSIPGPenalty, AlwaysPositive) {
  float gamma = computeSIPGPenalty<1>(kUnitFace, kUnitHex, 1.0f);
  EXPECT_GT(gamma, 0.0f);
}

TEST(ComputeSIPGPenalty, InverseLengthScaling) {
  // Doubling all lengths doubles h_f = vol/area (vol ~ h^3, area ~ h^2),
  // so gamma should halve.
  float face2[4][3] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}};
  float hex2[8][3] = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {2, 2, 0}, {0, 0, 2}, {2, 0, 2}, {0, 2, 2}, {2, 2, 2}};
  float g1 = computeSIPGPenalty<1>(kUnitFace, kUnitHex, 10.0f);
  float g2 = computeSIPGPenalty<1>(face2, hex2, 10.0f);
  EXPECT_NEAR(g2, g1 * 0.5f, 1e-5f);
}

}  // namespace test
}  // namespace fe
}  // namespace solver
