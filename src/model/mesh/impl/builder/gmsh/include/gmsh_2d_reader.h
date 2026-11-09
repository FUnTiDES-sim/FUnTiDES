#ifndef FUNTIDES_MODEL_GMSH_2D_READER_H_
#define FUNTIDES_MODEL_GMSH_2D_READER_H_

#include <filesystem>
#include <map>

#include "model_unstruct.h"

namespace model {
struct GmshAcousticMaterial {
  float vp;
  float rho{1};
};

/// Read a planar XZ mesh from Gmsh MSH 2.2 ASCII. Only linear quadrilaterals
/// (type 3), boundary lines (type 1) and point records (type 15) are accepted.
/// Vertex tags can be sparse. Rebuilds conforming tensor-product GLL nodes at
/// the requested SEM order; Gmsh high-order nodes are NOT used as GLL nodes.
/// Physical surface tags select materials; physical line tags select boundary
/// conditions. Every external quad edge must have a boundary line, and internal
/// edges must be shared by exactly two elements (no hanging-node interfaces).
/// CW elements are reoriented; folded, nonconvex, nonplanar and unsupported
/// elements are rejected. Pressure-free boundary tags take priority at corners.
ModelUnstructData<float, int> readGmsh2D(const std::filesystem::path& file, int order,
                                         const std::map<int, GmshAcousticMaterial>& materials,
                                         const std::map<int, BoundaryFlag>& boundaries);
}  // namespace model
#endif
