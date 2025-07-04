#ifndef MODEL_H_
#define MODEL_H_

#ifdef USE_SIMPLE_MESH
#include <SEMmesh.hpp>
using Mesh = SEMmesh<float, int, int>;
#endif

#endif // MODEL_H_
