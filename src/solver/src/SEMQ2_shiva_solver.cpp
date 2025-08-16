#include "SEMsolver_impl.hpp"

constexpr int ORDER = 2;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::SHIVA>::type, CartesianSEMmesh<float, int, ORDER> >;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::SHIVA>::type, CartesianUnstructMesh<float, int, ORDER> >;
