#ifndef SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_STRUCT_H_
#define SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_STRUCT_H_

#include "face_connectivity.h"

namespace model
{

/**
 * @brief Face connectivity implementation for structured Cartesian meshes
 *
 * Computes all face data using index formulas instead of a hash map.
 * Owner/neighbor relationships are derived deterministically from
 * the Cartesian topology.
 *
 * @tparam FloatType Floating point type
 * @tparam ScalarType Integer type for indexing
 * @tparam Order Polynomial order of spectral elements
 */
template <typename FloatType, typename ScalarType, int Order>
class FaceConnectivityStruct : public FaceConnectivityApi<FloatType, ScalarType>
{
 public:
  FaceConnectivityStruct() = default;

  FaceConnectivityStruct(ScalarType ex, ScalarType ey, ScalarType ez)
      : ex_(ex), ey_(ey), ez_(ez)
  {
    num_faces_x_ = (ex_ + 1) * ey_ * ez_;
    num_faces_y_ = ex_ * (ey_ + 1) * ez_;
    offset_y_ = num_faces_x_;
    offset_z_ = num_faces_x_ + num_faces_y_;
  }

  FaceConnectivity<FloatType, ScalarType> build(
      const ModelApi<FloatType, ScalarType>& mesh) const override
  {
    const ScalarType n_element = mesh.getNumberOfElements();
    const int ndofs_per_face = (Order + 1) * (Order + 1);
    const ScalarType n_faces = offset_z_ + ex_ * ey_ * (ez_ + 1);

    FaceConnectivity<FloatType, ScalarType> result;
    result.ndofs_per_face_ = ndofs_per_face;
    result.n_faces_ = n_faces;

    result.elem_to_faces_ = allocateArray2D<ARRAY_INT_VIEW>(n_element, 6);
    result.face_dofs_ =
        allocateArray2D<ARRAY_INT_VIEW>(n_faces, ndofs_per_face);
    result.face_elem_owner_ = allocateVector<VECTOR_INT_VIEW>(n_faces);
    result.face_elem_neighbor_ = allocateVector<VECTOR_INT_VIEW>(n_faces);
    result.face_local_owner_ = allocateVector<VECTOR_INT_VIEW>(n_faces);
    result.face_local_neighbor_ = allocateVector<VECTOR_INT_VIEW>(n_faces);

    // Initialize all neighbors to -1 (boundary marker)
    for (ScalarType f = 0; f < n_faces; ++f) result.face_elem_neighbor_(f) = -1;

    // -----------------------------------------------------------------------
    // Boucle 1 : elem_to_faces + face_dofs
    // -----------------------------------------------------------------------
    for (ScalarType elem = 0; elem < n_element; ++elem)
    {
      ScalarType elem_k = elem / (ex_ * ey_);
      ScalarType tmp = elem % (ex_ * ey_);
      ScalarType elem_j = tmp / ex_;
      ScalarType elem_i = tmp % ex_;

      // elem_to_faces via formules cartésiennes
      result.elem_to_faces_(elem, static_cast<int>(CubicFace::kXMinus)) =
          elem_i + elem_j * (ex_ + 1) + elem_k * (ex_ + 1) * ey_;
      result.elem_to_faces_(elem, static_cast<int>(CubicFace::kXPlus)) =
          (elem_i + 1) + elem_j * (ex_ + 1) + elem_k * (ex_ + 1) * ey_;
      result.elem_to_faces_(elem, static_cast<int>(CubicFace::kYMinus)) =
          offset_y_ + elem_i + elem_j * ex_ + elem_k * ex_ * (ey_ + 1);
      result.elem_to_faces_(elem, static_cast<int>(CubicFace::kYPlus)) =
          offset_y_ + elem_i + (elem_j + 1) * ex_ + elem_k * ex_ * (ey_ + 1);
      result.elem_to_faces_(elem, static_cast<int>(CubicFace::kZMinus)) =
          offset_z_ + elem_i + elem_j * ex_ + elem_k * ex_ * ey_;
      result.elem_to_faces_(elem, static_cast<int>(CubicFace::kZPlus)) =
          offset_z_ + elem_i + elem_j * ex_ + (elem_k + 1) * ex_ * ey_;

      // face_dofs pour chaque face locale
      for (int lf = 0; lf < 6; ++lf)
      {
        CubicFace local_face = static_cast<CubicFace>(lf);
        ScalarType face_id = result.elem_to_faces_(elem, lf);
        int idx = 0;

        switch (local_face)
        {
          case CubicFace::kXMinus:
            for (int k = 0; k <= Order; ++k)
              for (int j = 0; j <= Order; ++j)
                result.face_dofs_(face_id, idx++) =
                    mesh.globalNodeIndex(elem, 0, j, k);
            break;
          case CubicFace::kXPlus:
            for (int k = 0; k <= Order; ++k)
              for (int j = 0; j <= Order; ++j)
                result.face_dofs_(face_id, idx++) =
                    mesh.globalNodeIndex(elem, Order, j, k);
            break;
          case CubicFace::kYMinus:
            for (int k = 0; k <= Order; ++k)
              for (int i = 0; i <= Order; ++i)
                result.face_dofs_(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, 0, k);
            break;
          case CubicFace::kYPlus:
            for (int k = 0; k <= Order; ++k)
              for (int i = 0; i <= Order; ++i)
                result.face_dofs_(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, Order, k);
            break;
          case CubicFace::kZMinus:
            for (int j = 0; j <= Order; ++j)
              for (int i = 0; i <= Order; ++i)
                result.face_dofs_(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, j, 0);
            break;
          case CubicFace::kZPlus:
            for (int j = 0; j <= Order; ++j)
              for (int i = 0; i <= Order; ++i)
                result.face_dofs_(face_id, idx++) =
                    mesh.globalNodeIndex(elem, i, j, Order);
            break;
        }
      }
    }

    // -----------------------------------------------------------------------
    // Boucle 2 : owner/neighbor de façon déterministe
    //
    // Principe : on ne traite que les faces kXMinus/kYMinus/kZMinus
    // car chaque face interne est partagée par exactement 2 éléments :
    //   kXPlus(elem_i)  == kXMinus(elem_i+1)  → traité via elem_i+1
    //   kYPlus(elem_j)  == kYMinus(elem_j+1)  → traité via elem_j+1
    //   kZPlus(elem_k)  == kZMinus(elem_k+1)  → traité via elem_k+1
    // Les faces Plus en bord de domaine restent à neighbor=-1 (boundary)
    // -----------------------------------------------------------------------
    for (ScalarType elem = 0; elem < n_element; ++elem)
    {
      ScalarType elem_k = elem / (ex_ * ey_);
      ScalarType tmp = elem % (ex_ * ey_);
      ScalarType elem_j = tmp / ex_;
      ScalarType elem_i = tmp % ex_;

      // kXMinus : owner=elem, neighbor=elem(i-1) si existe
      {
        ScalarType face =
            result.elem_to_faces_(elem, static_cast<int>(CubicFace::kXMinus));
        result.face_elem_owner_(face) = elem;
        result.face_local_owner_(face) = static_cast<int>(CubicFace::kXMinus);
        if (elem_i > 0)
        {
          ScalarType neighbor =
              (elem_i - 1) + elem_j * ex_ + elem_k * ex_ * ey_;
          result.face_elem_neighbor_(face) = neighbor;
          result.face_local_neighbor_(face) =
              static_cast<int>(CubicFace::kXPlus);
        }
      }

      // kYMinus : owner=elem, neighbor=elem(j-1) si existe
      {
        ScalarType face =
            result.elem_to_faces_(elem, static_cast<int>(CubicFace::kYMinus));
        result.face_elem_owner_(face) = elem;
        result.face_local_owner_(face) = static_cast<int>(CubicFace::kYMinus);
        if (elem_j > 0)
        {
          ScalarType neighbor =
              elem_i + (elem_j - 1) * ex_ + elem_k * ex_ * ey_;
          result.face_elem_neighbor_(face) = neighbor;
          result.face_local_neighbor_(face) =
              static_cast<int>(CubicFace::kYPlus);
        }
      }

      // kZMinus : owner=elem, neighbor=elem(k-1) si existe
      {
        ScalarType face =
            result.elem_to_faces_(elem, static_cast<int>(CubicFace::kZMinus));
        result.face_elem_owner_(face) = elem;
        result.face_local_owner_(face) = static_cast<int>(CubicFace::kZMinus);
        if (elem_k > 0)
        {
          ScalarType neighbor =
              elem_i + elem_j * ex_ + (elem_k - 1) * ex_ * ey_;
          result.face_elem_neighbor_(face) = neighbor;
          result.face_local_neighbor_(face) =
              static_cast<int>(CubicFace::kZPlus);
        }
      }
    }

    return result;
  }

 private:
  ScalarType ex_{0}, ey_{0}, ez_{0};

  // Pre-computed face offsets
  ScalarType num_faces_x_{0};
  ScalarType num_faces_y_{0};
  ScalarType offset_y_{0};
  ScalarType offset_z_{0};
};

}  // namespace model

#endif  // SRC_MODEL_MESH_IMPL_FACE_CONNECTIVITY_STRUCT_H_