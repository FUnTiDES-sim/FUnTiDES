#include "SEMsolver_impl.hpp"

constexpr int ORDER = 3;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::OPTIM >::type, CartesianSEMmesh<float, int, ORDER> >;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::OPTIM>::type, CartesianUnstructMesh<float, int, ORDER> >;
