#ifndef SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_UNSTRUCT_H_
#define SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_UNSTRUCT_H_

#include <algorithm>
#include <array>
#include <map>

#include "face_connectivity.h"

namespace model
{

/**
 * @brief Face connectivity for unstructured meshes
 * @tparam FloatType Floating point type
 * @tparam ScalarType Integer type for indexing
 */
template <typename FloatType, typename ScalarType>
class FaceConnectivityUnstruct
    : public FaceConnectivityApi<FloatType, ScalarType>
{
 public:
  FaceConnectivityUnstruct() = default;

  /**
   * @brief Build face connectivity from the mesh
   * - Extracts faces from elements, identifies unique faces, and fills
   * connectivity tables.
   * - Uses a map to identify unique faces based on their corner nodes.
   * - Stores: * - elem_to_faces_(elem, local_face) → global_face_id *
   * - face_dofs_(face_id, local_dof) → global_node_id *
   * - face_elem_owner_(face_id) → owner element index *
   * - face_elem_neighbor_(face_id) → neighbor element index or -1 if boundary *
   * - face_local_owner_(face_id) → local face index in owner element *
   * - face_local_neighbor_(face_id) → local face index in neighbor element or
   * -1 if boundary
   */
  void build(const ModelApi<FloatType, ScalarType>& mesh)
  {
    const ScalarType n_element = mesh.getNumberOfElements();
    const int order = mesh.getOrder();
    const ScalarType max_faces = n_element * 6;
    ndofs_per_face_ = (order + 1) * (order + 1);

    auto elem_to_faces_temp = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
    auto face_dofs_temp =
        allocateArray2D<ARRAY_INT_VIEW>(max_faces, ndofs_per_face_);
    auto face_elem_owner_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
    auto face_elem_neighbor_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
    auto face_local_owner_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);
    auto face_local_neighbor_temp = allocateVector<VECTOR_INT_VIEW>(max_faces);

    for (ScalarType i = 0; i < max_faces; ++i) face_elem_neighbor_temp(i) = -1;

    using FaceKey = std::array<ScalarType, 4>;
    std::map<FaceKey, ScalarType> face_map;
    ScalarType face_count = 0;

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
          ScalarType face_id = face_count++;
          face_map[face_key] = face_id;

          // Remplissage face_dofs
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
          ScalarType face_id = it->second;
          face_elem_neighbor_temp(face_id) = elem;
          face_local_neighbor_temp(face_id) = lf;
          elem_to_faces_temp(elem, lf) = face_id;
        }
      }
    }

    n_faces_ = face_count;
    elem_to_faces_ = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
    face_dofs_ = allocateArray2D<ARRAY_INT_VIEW>(face_count, ndofs_per_face_);
    face_elem_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    face_elem_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    face_local_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    face_local_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);

    for (ScalarType elem = 0; elem < n_element; ++elem)
      for (int lf = 0; lf < 6; ++lf)
        elem_to_faces_(elem, lf) = elem_to_faces_temp(elem, lf);

    for (ScalarType f = 0; f < face_count; ++f)
    {
      face_elem_owner_(f) = face_elem_owner_temp(f);
      face_elem_neighbor_(f) = face_elem_neighbor_temp(f);
      face_local_owner_(f) = face_local_owner_temp(f);
      face_local_neighbor_(f) = face_local_neighbor_temp(f);
      for (int dof = 0; dof < ndofs_per_face_; ++dof)
        face_dofs_(f, dof) = face_dofs_temp(f, dof);
    }
  }

  /**
   * @brief Get total number of unique faces in the mesh
   * @return Number of faces
   */
  PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const override
  {
    return n_faces_;
  }

  /**
   * @brief Get number of DOFs (nodes) per face
   * @return Number of DOFs per face
   */
  PROXY_HOST_DEVICE int getDofsPerFace() const override
  {
    return ndofs_per_face_;
  }

  /**
   * @brief Get global face ID from element and local face
   * @param elem Element index
   * @param local_face Local face identifier
   * @return Global face ID
   */
  PROXY_HOST_DEVICE ScalarType
  getGlobalFace(ScalarType elem, CubicFace local_face) const override
  {
    return elem_to_faces_(elem, static_cast<int>(local_face));
  }

  /**
   * @brief Get global node index from face and local DOF
   * @param face_id Global face ID
   * @param local_dof Local DOF index on face [0, (order+1)²)
   * @return Global node index
   */
  PROXY_HOST_DEVICE ScalarType
  getGlobalNodeFromFace(ScalarType face_id, int local_dof) const override
  {
    return face_dofs_(face_id, local_dof);
  }

  /**
   * @brief Check if face is on domain boundary (no neighbor)
   * @param face_id Global face ID
   * @return True if boundary face
   */
  PROXY_HOST_DEVICE bool isBoundaryFace(ScalarType face_id) const override
  {
    return face_elem_neighbor_(face_id) == -1;
  }

  /**
   * @brief Get owner element of a face
   * @param face_id Global face ID
   * @return Owner element index
   */
  PROXY_HOST_DEVICE ScalarType elemOwner(ScalarType face_id) const override
  {
    return face_elem_owner_(face_id);
  }

  /**
   * @brief Get neighbor element of a face (-1 if boundary)
   * @param face_id Global face ID
   * @return Neighbor element index, or -1 if boundary
   */
  PROXY_HOST_DEVICE ScalarType elemNeighbor(ScalarType face_id) const override
  {
    return face_elem_neighbor_(face_id);
  }

  /**
   * @brief Get local face index of the owner element
   * @param face_id Global face ID
   * @return Local face index (0-5) in owner element
   */
  PROXY_HOST_DEVICE int localFaceOwner(ScalarType face_id) const override
  {
    return face_local_owner_(face_id);
  }

  /**
   * @brief Get local face index of the neighbor element (-1 if boundary)
   * @param face_id Global face ID
   * @return Local face index (0-5) in neighbor element, or -1 if boundary
   */
  PROXY_HOST_DEVICE int localFaceNeighbor(ScalarType face_id) const override
  {
    return face_local_neighbor_(face_id);
  }

 private:
  ScalarType n_faces_ = 0;
  int ndofs_per_face_ = 0;

  ARRAY_INT_VIEW elem_to_faces_;
  ARRAY_INT_VIEW face_dofs_;
  VECTOR_INT_VIEW face_elem_owner_;
  VECTOR_INT_VIEW face_elem_neighbor_;
  VECTOR_INT_VIEW face_local_owner_;
  VECTOR_INT_VIEW face_local_neighbor_;

  // ------------------------------------------------------------------
  // Helpers for the build process
  // ------------------------------------------------------------------

  static std::array<ScalarType, 4> extractFaceCorners(
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

  static std::array<ScalarType, 4> makeFaceKey(
      std::array<ScalarType, 4> corners)
  {
    std::sort(corners.begin(), corners.end());
    return corners;
  }
};

}  // namespace model

#endif  // SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_UNSTRUCT_H_