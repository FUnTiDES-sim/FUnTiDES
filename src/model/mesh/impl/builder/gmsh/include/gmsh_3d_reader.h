#ifndef FUNTIDES_MODEL_GMSH_3D_READER_H_
#define FUNTIDES_MODEL_GMSH_3D_READER_H_
#include "gmsh_2d_reader.h"
namespace model {
/// Import conforming straight-sided Hex8 volumes (MSH 2.2 ASCII, type 5).
/// Physical volumes select acoustic materials; physical Q4 surfaces select boundaries.
/// Lines and points are ignored. All external faces must be marked. Internal
/// physical surfaces, mixed cells, high-order geometry and inverted cells are rejected.
/// Reconstructs tensor-product GLL nodes, sharing vertices, edges and faces by topology.
/// Positive Jacobians are checked on a (2*order+1)^3 reference grid (not a proof
/// of positivity everywhere for a general warped trilinear cell).
ModelUnstructData<float, int> readGmsh3D(const std::filesystem::path& file, int order,
                                         const std::map<int, GmshAcousticMaterial>& materials,
                                         const std::map<int, BoundaryFlag>& boundaries);
}  // namespace model
#endif
