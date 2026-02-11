#ifndef SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_UNSTRUCT_H_
#define SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_UNSTRUCT_H_

#include <algorithm>
#include <array>
#include <map>

#include "face_connectivity.h"

namespace model
{

template <typename FloatType, typename ScalarType>
class FaceConnectivityUnstruct
    : public FaceConnectivityApi<FloatType, ScalarType>
{
 public:
  FaceConnectivity<FloatType, ScalarType> build(
      const ModelApi<FloatType, ScalarType>& mesh) const override
  {
    FaceConnectivity<FloatType, ScalarType> result;

    const ScalarType n_element = mesh.getNumberOfElements();
    const int order = mesh.getOrder();
    const ScalarType max_faces = n_element * 6;
    const int ndofs_per_face = (order + 1) * (order + 1);

    result.ndofs_per_face_ = ndofs_per_face;

    auto elem_to_faces_temp = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
    auto face_dofs_temp =
        allocateArray2D<ARRAY_INT_VIEW>(max_faces, ndofs_per_face);
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

    // Allocate final arrays with exact size
    result.n_faces_ = face_count;
    result.elem_to_faces_ = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
    result.face_dofs_ =
        allocateArray2D<ARRAY_INT_VIEW>(face_count, ndofs_per_face);
    result.face_elem_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    result.face_elem_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    result.face_local_owner_ = allocateVector<VECTOR_INT_VIEW>(face_count);
    result.face_local_neighbor_ = allocateVector<VECTOR_INT_VIEW>(face_count);

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

 private:
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