#pragma once
#include <gtest/gtest.h>

#include <set>

#include "face_connectivity.h"

namespace model
{

template <typename FloatType, typename ScalarType, typename FCType>
void testSingleCubeFaceCount(const FCType& fc)
{
  EXPECT_EQ(fc.getNumberOfFaces(), 6);
  EXPECT_EQ(fc.getDofsPerFace(), 4);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testTwoAdjacentCubesFaceCount(const FCType& fc)
{
  EXPECT_EQ(fc.getNumberOfFaces(), 11);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testGlobalFaceIds(const FCType& fc)
{
  for (int lf = 0; lf < 6; ++lf)
  {
    auto face_id = fc.getGlobalFace(0, static_cast<CubicFace>(lf));
    EXPECT_GE(face_id, 0);
    EXPECT_LT(face_id, 6);
  }
}

template <typename FloatType, typename ScalarType, typename FCType>
void testGlobalFaceUniqueness(const FCType& fc)
{
  std::set<ScalarType> face_ids;
  for (int lf = 0; lf < 6; ++lf)
    face_ids.insert(fc.getGlobalFace(0, static_cast<CubicFace>(lf)));
  EXPECT_EQ(face_ids.size(), 6);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testSingleCubeAllBoundary(const FCType& fc)
{
  for (ScalarType f = 0; f < fc.getNumberOfFaces(); ++f)
    EXPECT_TRUE(fc.isBoundaryFace(f));
}

template <typename FloatType, typename ScalarType, typename FCType>
void testSharedFaceNotBoundary(const FCType& fc)
{
  auto face_elem0_xplus = fc.getGlobalFace(0, CubicFace::kXPlus);
  auto face_elem1_xminus = fc.getGlobalFace(1, CubicFace::kXMinus);
  EXPECT_EQ(face_elem0_xplus, face_elem1_xminus);
  EXPECT_FALSE(fc.isBoundaryFace(face_elem0_xplus));
}

template <typename FloatType, typename ScalarType, typename FCType>
void testBoundaryFaceCount(const FCType& fc)
{
  int boundary_count = 0;
  for (ScalarType f = 0; f < fc.getNumberOfFaces(); ++f)
    if (fc.isBoundaryFace(f)) boundary_count++;
  EXPECT_EQ(boundary_count, 10);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testFaceNodes(const FCType& fc, int n_nodes)
{
  auto face0 = fc.getGlobalFace(0, CubicFace::kXMinus);
  for (int dof = 0; dof < 4; ++dof)
  {
    auto node = fc.getGlobalNodeFromFace(face0, dof);
    EXPECT_GE(node, 0);
    EXPECT_LT(node, n_nodes);
  }
}

template <typename FloatType, typename ScalarType, typename FCType>
void testFaceNodesUnique(const FCType& fc)
{
  auto face0 = fc.getGlobalFace(0, CubicFace::kXMinus);
  std::set<ScalarType> nodes;
  for (int dof = 0; dof < 4; ++dof)
    nodes.insert(fc.getGlobalNodeFromFace(face0, dof));
  EXPECT_EQ(nodes.size(), 4);
}

template <typename FloatType, typename ScalarType, typename FCType>
void testInternalFaceOwners(const FCType& fc)
{
  auto shared_face = fc.getGlobalFace(0, CubicFace::kXPlus);
  EXPECT_GE(fc.elemOwner(shared_face), 0);
  EXPECT_GE(fc.elemNeighbor(shared_face), 0);
  EXPECT_NE(fc.elemOwner(shared_face), fc.elemNeighbor(shared_face));
}

template <typename FloatType, typename ScalarType, typename FCType>
void testBoundaryFaceNoNeighbor(const FCType& fc)
{
  auto boundary_face = fc.getGlobalFace(0, CubicFace::kXMinus);
  EXPECT_GE(fc.elemOwner(boundary_face), 0);
  EXPECT_EQ(fc.elemNeighbor(boundary_face), -1);
}

}  // namespace model