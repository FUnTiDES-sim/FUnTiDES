#pragma once
#include <gtest/gtest.h>

#include <set>

#include "face_connectivity.h"

namespace model {

template <typename FloatType, typename ScalarType, typename FCType>
void testSingleCubeFaceCount(const FCType& fc) {
  EXPECT_EQ(fc.getNumberOfFaces(), 6);
  EXPECT_EQ(fc.getDofsPerFace(), 4);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testTwoAdjacentCubesFaceCount(const FCType& fc) {
  EXPECT_EQ(fc.getNumberOfFaces(), 11);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testGlobalFaceIds(const FCType& fc) {
  for (int lf = 0; lf < 6; ++lf) {
    auto face_id = fc.getGlobalFace(0, static_cast<CubicFace>(lf));
    EXPECT_GE(face_id, 0);
    EXPECT_LT(face_id, 6);
  }
}

template <typename FloatType, typename ScalarType, typename FCType>
void testGlobalFaceUniqueness(const FCType& fc) {
  std::set<ScalarType> face_ids;
  for (int lf = 0; lf < 6; ++lf) face_ids.insert(fc.getGlobalFace(0, static_cast<CubicFace>(lf)));
  EXPECT_EQ(face_ids.size(), 6);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testSingleCubeAllBoundary(const FCType& fc) {
  for (ScalarType f = 0; f < fc.getNumberOfFaces(); ++f) EXPECT_TRUE(fc.isBoundaryFace(f));
}

template <typename FloatType, typename ScalarType, typename FCType>
void testSharedFaceNotBoundary(const FCType& fc) {
  auto face_elem0_xplus = fc.getGlobalFace(0, CubicFace::kXPlus);
  auto face_elem1_xminus = fc.getGlobalFace(1, CubicFace::kXMinus);
  EXPECT_EQ(face_elem0_xplus, face_elem1_xminus);
  EXPECT_FALSE(fc.isBoundaryFace(face_elem0_xplus));
}

template <typename FloatType, typename ScalarType, typename FCType>
void testBoundaryFaceCount(const FCType& fc) {
  int boundary_count = 0;
  for (ScalarType f = 0; f < fc.getNumberOfFaces(); ++f)
    if (fc.isBoundaryFace(f)) boundary_count++;
  EXPECT_EQ(boundary_count, 10);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testFaceNodes(const FCType& fc, int n_nodes) {
  auto face0 = fc.getGlobalFace(0, CubicFace::kXMinus);
  for (int dof = 0; dof < 4; ++dof) {
    auto node = fc.getGlobalNodeFromFace(face0, dof);
    EXPECT_GE(node, 0);
    EXPECT_LT(node, n_nodes);
  }
}

template <typename FloatType, typename ScalarType, typename FCType>
void testFaceNodesUnique(const FCType& fc) {
  auto face0 = fc.getGlobalFace(0, CubicFace::kXMinus);
  std::set<ScalarType> nodes;
  for (int dof = 0; dof < 4; ++dof) nodes.insert(fc.getGlobalNodeFromFace(face0, dof));
  EXPECT_EQ(nodes.size(), 4);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testInternalFaceOwners(const FCType& fc) {
  auto shared_face = fc.getGlobalFace(0, CubicFace::kXPlus);
  EXPECT_GE(fc.elemOwner(shared_face), 0);
  EXPECT_GE(fc.elemNeighbor(shared_face), 0);
  EXPECT_NE(fc.elemOwner(shared_face), fc.elemNeighbor(shared_face));
}

template <typename FloatType, typename ScalarType, typename FCType>
void testBoundaryFaceNoNeighbor(const FCType& fc) {
  auto boundary_face = fc.getGlobalFace(0, CubicFace::kXMinus);
  EXPECT_GE(fc.elemOwner(boundary_face), 0);
  EXPECT_EQ(fc.elemNeighbor(boundary_face), -1);
}

/// Reconstruct the face DOFs for a given element and local face, following
/// the same loop convention as FaceConnectivityUnstruct::build().
template <typename FloatType, typename ScalarType>
std::vector<ScalarType> computeFaceDofs(const ModelApi<FloatType, ScalarType>& mesh, ScalarType elem, CubicFace face,
                                        int order) {
  const int ndofs = (order + 1) * (order + 1);
  std::vector<ScalarType> dofs(ndofs);
  int idx = 0;
  switch (face) {
    case CubicFace::kXMinus:
      for (int k = 0; k <= order; ++k)
        for (int j = 0; j <= order; ++j) dofs[idx++] = mesh.globalNodeIndex(elem, 0, j, k);
      break;
    case CubicFace::kXPlus:
      for (int k = 0; k <= order; ++k)
        for (int j = 0; j <= order; ++j) dofs[idx++] = mesh.globalNodeIndex(elem, order, j, k);
      break;
    case CubicFace::kYMinus:
      for (int k = 0; k <= order; ++k)
        for (int i = 0; i <= order; ++i) dofs[idx++] = mesh.globalNodeIndex(elem, i, 0, k);
      break;
    case CubicFace::kYPlus:
      for (int k = 0; k <= order; ++k)
        for (int i = 0; i <= order; ++i) dofs[idx++] = mesh.globalNodeIndex(elem, i, order, k);
      break;
    case CubicFace::kZMinus:
      for (int j = 0; j <= order; ++j)
        for (int i = 0; i <= order; ++i) dofs[idx++] = mesh.globalNodeIndex(elem, i, j, 0);
      break;
    case CubicFace::kZPlus:
      for (int j = 0; j <= order; ++j)
        for (int i = 0; i <= order; ++i) dofs[idx++] = mesh.globalNodeIndex(elem, i, j, order);
      break;
  }
  return dofs;
}

/// Verify that getNeighborFaceDof produces a correct permutation: for each
/// owner DOF i, the global node at owner's i equals the global node at
/// the neighbor's perm(i).
template <typename FloatType, typename ScalarType, typename FCType>
void testNeighborFaceDofCorrectness(const FCType& fc, const ModelApi<FloatType, ScalarType>& mesh,
                                    ScalarType shared_face) {
  const int ndofs = fc.getDofsPerFace();
  const int order = mesh.getOrder();
  const ScalarType e1 = fc.elemNeighbor(shared_face);
  const CubicFace lf1 = static_cast<CubicFace>(fc.localFaceNeighbor(shared_face));

  auto neigh_dofs = computeFaceDofs<FloatType, ScalarType>(mesh, e1, lf1, order);

  for (int i = 0; i < ndofs; ++i) {
    ScalarType owner_node = fc.getGlobalNodeFromFace(shared_face, i);
    int perm_i = fc.getNeighborFaceDof(shared_face, i);
    EXPECT_GE(perm_i, 0);
    EXPECT_LT(perm_i, ndofs);
    EXPECT_EQ(owner_node, neigh_dofs[perm_i]) << "owner DOF " << i << " (node " << owner_node << ") does not match"
                                              << " neighbor DOF " << perm_i << " (node " << neigh_dofs[perm_i] << ")";
  }

  // Permutation must be a bijection
  std::set<int> seen;
  for (int i = 0; i < ndofs; ++i) seen.insert(fc.getNeighborFaceDof(shared_face, i));
  EXPECT_EQ(static_cast<int>(seen.size()), ndofs) << "getNeighborFaceDof is not a bijection";
}

}  // namespace model