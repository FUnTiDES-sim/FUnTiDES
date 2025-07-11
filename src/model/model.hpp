#ifndef MODEL_H_
#define MODEL_H_

#ifdef USE_SIMPLE_MESH
#include <SEMmesh.hpp>
using Mesh = CartesianSEMmesh<float, int, int, 2>;
#endif

#endif // MODEL_H_
