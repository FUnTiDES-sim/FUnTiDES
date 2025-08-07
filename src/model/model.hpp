#ifndef MODEL_H_
#define MODEL_H_

#ifdef USE_SIMPLE_MESH
#include <cartesianSEMmesh.hpp>
using Mesh = CartesianSEMmesh<float, float, int, int, 2>;
#endif

#endif // MODEL_H_
