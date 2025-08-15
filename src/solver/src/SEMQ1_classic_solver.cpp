#include "SEMsolver_impl.hpp"

constexpr int ORDER = 1;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::CLASSIC >::type, CartesianSEMmesh<float, int, ORDER> >;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::CLASSIC >::type, CartesianUnstructMesh<float, int, ORDER> >;
