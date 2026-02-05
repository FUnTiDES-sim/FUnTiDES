#ifndef SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_H_
#define SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_H_

#include <model.h>

#include <algorithm>
#include <array>
#include <map>

namespace model
{

/**
 * @brief GPU-compatible face connectivity structure
 *
 * Independent topology manager for hexahedral meshes.
 * Works with any mesh implementing the ModelApi interface.
 */
template <typename FloatType, typename ScalarType>
class FaceConnectivity
{
 public:
  /**
   * @brief Build face connectivity from any mesh
   *
   * @param mesh Mesh object (ModelStruct or ModelUnstruct)
   * @return FaceConnectivity object with all face tables constructed
   */
  static FaceConnectivity build(const ModelApi<FloatType, ScalarType>& mesh);

  *@brief Get total number of unique faces in
          the mesh* @ return Number of faces *
      / PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const
  {
    return n_faces_;
  }

  /**
   * @brief Get number of degrees of freedom (nodes) per face
   * @return Number of DOFs per face: (order+1)²
   */
  PROXY_HOST_DEVICE ScalarType getDofsPerFace() const
  {
    return ndofs_per_face_;
  }

  /**
   * @brief Map element and local face to global face ID
   * @param elem Element index
   * @param local_face Local face identifier (kXMinus, kXPlus, etc.)
   * @return Global face ID [0, n_faces)
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const
  {
    return elem_to_faces_(elem, static_cast<int>(local_face));
  }

  /**
   * @brief Get global node index from face and local DOF
   * @param face_id Global face ID
   * @param local_dof Local node index on face [0, ndofs_per_face)
   * @return Global node index
   */
  PROXY_HOST_DEVICE
  ScalarType getGlobalNodeFromFace(ScalarType face_id, int local_dof) const
  {
    return face_dofs_(face_id, local_dof);
  }

  /**
   * @brief Check if face is on domain boundary
   * @param face_id Global face ID
   * @return True if boundary face (no neighbor element), false if internal
   */
  PROXY_HOST_DEVICE
  bool isBoundaryFace(ScalarType face_id) const
  {
    return face_elem_neighbor_(face_id) == -1;
  }

  /**
   * @brief Get element that owns this face
   * @param face_id Global face ID
   * @return Element index of owner
   * @note Every face has exactly one owner
   */
  PROXY_HOST_DEVICE ScalarType elemOwner(ScalarType face_id) const
  {
    return face_elem_owner_(face_id);
  }

  /**
   * @brief Get neighboring element across this face
   * @param face_id Global face ID
   * @return Element index of neighbor, or -1 if boundary face
   */
  PROXY_HOST_DEVICE ScalarType elemNeighbor(ScalarType face_id) const
  {
    return face_elem_neighbor_(face_id);
  }

  /**
   * @brief Get local face index in owner element
   * @param face_id Global face ID
   * @return Local face index [0-5] in owner element
   * @see CubicFace enum for face numbering convention
   */
  PROXY_HOST_DEVICE int localFaceOwner(ScalarType face_id) const
  {
    return face_local_owner_(face_id);
  }

  /**
   * @brief Get local face index in neighbor element
   * @param face_id Global face ID
   * @return Local face index [0-5] in neighbor element, or -1 if boundary
   * @see CubicFace enum for face numbering convention
   */
  PROXY_HOST_DEVICE int localFaceNeighbor(ScalarType face_id) const
  {
    return face_local_neighbor_(face_id);
  }

 private:
  // Storage
  ARRAY_INT_VIEW elem_to_faces_;
  ARRAY_INT_VIEW face_dofs_;
  VECTOR_INT_VIEW face_elem_owner_;
  VECTOR_INT_VIEW face_elem_neighbor_;
  VECTOR_INT_VIEW face_local_owner_;
  VECTOR_INT_VIEW face_local_neighbor_;

  ScalarType n_faces_ = 0;
  int ndofs_per_face_ = 0;

  // Construction helpers
  static std::array<ScalarType, 4> extractFaceCorners(
      const ModelApi<FloatType, ScalarType>& mesh, ScalarType elem,
      CubicFace local_face);

  static std::array<ScalarType, 4> makeFaceKey(
      std::array<ScalarType, 4> corners);
};

// ============================================================================
// IMPLEMENTATION
// ============================================================================

template <typename FloatType, typename ScalarType>
FaceConnectivity<FloatType, ScalarType>
FaceConnectivity<FloatType, ScalarType>::build(
    const ModelApi<FloatType, ScalarType>& mesh)
{
  FaceConnectivity result;

  const ScalarType n_element = mesh.getNumberOfElements();
  const int order = mesh.getOrder();
  const ScalarType max_faces = n_element * 6;
  const int ndofs_per_face = (order + 1) * (order + 1);

  result.ndofs_per_face_ = ndofs_per_face;

  // Temporary arrays for construction
  auto elem_to_faces_temp = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
  auto face_dofs_temp =
      allocateArray2D<ARRAY_INT_VIEW>(max_faces, ndofs_per_face);
  auto face_elem_owner_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
  auto face_elem_neighbor_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
  auto face_local_owner_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
  auto face_local_neighbor_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);

  // Initialize neighbors to -1 (boundary marker)
  for (ScalarType i = 0; i < max_faces; ++i) face_elem_neighbor_temp(i) = -1;

  // Map for face uniqueness
  using FaceKey = std::array<ScalarType, 4>;
  std::map<FaceKey, ScalarType> face_map;

  ScalarType face_count = 0;

  // Build face connectivity
  for (ScalarType elem = 0; elem < n_element; ++elem)
  {
    for (int lf = 0; lf < 6; ++lf)
    {
      CubicFace local_face = static_cast<CubicFace>(lf);
      auto corners = extractFaceCorners(mesh, elem, local_face);
      auto face_key = makeFaceKey(corners);

      auto it = face_map.find(face_key);
      if (it == face_map.end())
      {
        // New face
        ScalarType face_id = face_count++;
        face_map[face_key] = face_id;

        // Fill face DOFs
        int idx = 0;
        switch (local_face)
        {
          case CubicFace::kXMinus:
            for (int k = 0; k <= order; ++k)
              for (int j = 0; j <= order; ++j)
                face_dofs_temp(face_id, idx++) =
                    mesh.globalNodeIndex(elem, 0, j, k);
            break;
          case CubicFace::kXPlus:
            for (int k = 0; k <= order; ++k)
              for (int j = 0; j <= order; ++j)
                face_dofs_temp(face_id, idx++) =
                    mesh.globalNodeIndex(elem, order, j, k);
            break;
          case CubicFace::kYMinus:
            for (int k = 0; k <= order; ++k)
              for (int i = 0; i <= order; ++i)
                face_dofs_temp(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, 0, k);
            break;
          case CubicFace::kYPlus:
            for (int k = 0; k <= order; ++k)
              for (int i = 0; i <= order; ++i)
                face_dofs_temp(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, order, k);
            break;
          case CubicFace::kZMinus:
            for (int j = 0; j <= order; ++j)
              for (int i = 0; i <= order; ++i)
                face_dofs_temp(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, j, 0);
            break;
          case CubicFace::kZPlus:
            for (int j = 0; j <= order; ++j)
              for (int i = 0; i <= order; ++i)
                face_dofs_temp(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, j, order);
            break;
        }

        face_elem_owner_temp(face_id) = elem;
        face_local_owner_temp(face_id) = lf;
        elem_to_faces_temp(elem, lf) = face_id;
      }
      else
      {
        // Face already seen (internal face)
        ScalarType face_id = it->second;
        face_elem_neighbor_temp(face_id) = elem;
        face_local_neighbor_temp(face_id) = lf;
        elem_to_faces_temp(elem, lf) = face_id;
      }
    }
  }

  // Allocate final arrays with exact size
  result.n_faces_ = face_count;

  result.elem_to_faces_ = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
  result.face_dofs_ =
      allocateArray2D<ARRAY_INT_VIEW>(face_count, ndofs_per_face);
  result.face_elem_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
  result.face_elem_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);
  result.face_local_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
  result.face_local_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);

  // Copy to final arrays
  for (ScalarType elem = 0; elem < n_element; ++elem)
    for (int lf = 0; lf < 6; ++lf)
      result.elem_to_faces_(elem, lf) = elem_to_faces_temp(elem, lf);

  for (ScalarType face_id = 0; face_id < face_count; ++face_id)
  {
    result.face_elem_owner_(face_id) = face_elem_owner_temp(face_id);
    result.face_elem_neighbor_(face_id) = face_elem_neighbor_temp(face_id);
    result.face_local_owner_(face_id) = face_local_owner_temp(face_id);
    result.face_local_neighbor_(face_id) = face_local_neighbor_temp(face_id);

    for (int dof = 0; dof < ndofs_per_face; ++dof)
      result.face_dofs_(face_id, dof) = face_dofs_temp(face_id, dof);
  }

  return result;
}

template <typename FloatType, typename ScalarType>
std::array<ScalarType, 4>
FaceConnectivity<FloatType, ScalarType>::extractFaceCorners(
    const ModelApi<FloatType, ScalarType>& mesh, ScalarType elem,
    CubicFace local_face)
{
  const int o = mesh.getOrder();

  switch (local_face)
  {
    case CubicFace::kXMinus:
      return {mesh.globalNodeIndex(elem, 0, 0, 0),
              mesh.globalNodeIndex(elem, 0, o, 0),
              mesh.globalNodeIndex(elem, 0, o, o),
              mesh.globalNodeIndex(elem, 0, 0, o)};
    case CubicFace::kXPlus:
      return {mesh.globalNodeIndex(elem, o, 0, 0),
              mesh.globalNodeIndex(elem, o, o, 0),
              mesh.globalNodeIndex(elem, o, o, o),
              mesh.globalNodeIndex(elem, o, 0, o)};
    case CubicFace::kYMinus:
      return {mesh.globalNodeIndex(elem, 0, 0, 0),
              mesh.globalNodeIndex(elem, o, 0, 0),
              mesh.globalNodeIndex(elem, o, 0, o),
              mesh.globalNodeIndex(elem, 0, 0, o)};
    case CubicFace::kYPlus:
      return {mesh.globalNodeIndex(elem, 0, o, 0),
              mesh.globalNodeIndex(elem, o, o, 0),
              mesh.globalNodeIndex(elem, o, o, o),
              mesh.globalNodeIndex(elem, 0, o, o)};
    case CubicFace::kZMinus:
      return {mesh.globalNodeIndex(elem, 0, 0, 0),
              mesh.globalNodeIndex(elem, o, 0, 0),
              mesh.globalNodeIndex(elem, o, o, 0),
              mesh.globalNodeIndex(elem, 0, o, 0)};
    case CubicFace::kZPlus:
      return {mesh.globalNodeIndex(elem, 0, 0, o),
              mesh.globalNodeIndex(elem, o, 0, o),
              mesh.globalNodeIndex(elem, o, o, o),
              mesh.globalNodeIndex(elem, 0, o, o)};
  }
  return {};
}

template <typename FloatType, typename ScalarType>
std::array<ScalarType, 4> FaceConnectivity<FloatType, ScalarType>::makeFaceKey(
    std::array<ScalarType, 4> corners)
{
  std::sort(corners.begin(), corners.end());
  return corners;
}

}  // namespace model

#endif  // SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_H_