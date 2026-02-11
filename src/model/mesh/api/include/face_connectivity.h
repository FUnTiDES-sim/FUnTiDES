#ifndef SRC_MODEL_MESH_API_FACE_CONNECTIVITY_H_
#define SRC_MODEL_MESH_API_FACE_CONNECTIVITY_H_

#include <model.h>

namespace model
{

// ============================================================================
// FACE CONNECTIVITY DATA (GPU-compatible)
// ============================================================================

template <typename FloatType, typename ScalarType>
class FaceConnectivity
{
 public:
  PROXY_HOST_DEVICE ScalarType getNumberOfFaces() const { return n_faces_; }
  PROXY_HOST_DEVICE ScalarType getDofsPerFace() const
  {
    return ndofs_per_face_;
  }

  PROXY_HOST_DEVICE
  ScalarType getGlobalFace(ScalarType elem, CubicFace local_face) const
  {
    return elem_to_faces_(elem, static_cast<int>(local_face));
  }

  PROXY_HOST_DEVICE
  ScalarType getGlobalNodeFromFace(ScalarType face_id, int local_dof) const
  {
    return face_dofs_(face_id, local_dof);
  }

  PROXY_HOST_DEVICE
  bool isBoundaryFace(ScalarType face_id) const
  {
    return face_elem_neighbor_(face_id) == -1;
  }

  PROXY_HOST_DEVICE ScalarType elemOwner(ScalarType face_id) const
  {
    return face_elem_owner_(face_id);
  }
  PROXY_HOST_DEVICE ScalarType elemNeighbor(ScalarType face_id) const
  {
    return face_elem_neighbor_(face_id);
  }
  PROXY_HOST_DEVICE int localFaceOwner(ScalarType face_id) const
  {
    return face_local_owner_(face_id);
  }
  PROXY_HOST_DEVICE int localFaceNeighbor(ScalarType face_id) const
  {
    return face_local_neighbor_(face_id);
  }

  ARRAY_INT_VIEW elem_to_faces_;
  ARRAY_INT_VIEW face_dofs_;
  VECTOR_INT_VIEW face_elem_owner_;
  VECTOR_INT_VIEW face_elem_neighbor_;
  VECTOR_INT_VIEW face_local_owner_;
  VECTOR_INT_VIEW face_local_neighbor_;
  ScalarType n_faces_ = 0;
  int ndofs_per_face_ = 0;
};

template <typename FloatType, typename ScalarType>
class FaceConnectivityApi
{
 public:
  virtual ~FaceConnectivityApi() = default;

  virtual FaceConnectivity<FloatType, ScalarType> build(
      const ModelApi<FloatType, ScalarType>& mesh) const = 0;
};

}  // namespace model

#endif  // SRC_MODEL_MESH_API_FACE_CONNECTIVITY_H_