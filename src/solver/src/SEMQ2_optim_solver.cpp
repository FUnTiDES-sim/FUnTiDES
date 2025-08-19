#include "SEMsolver_impl.hpp"

constexpr int ORDER = 2;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::OPTIM >::type, CartesianSEMmesh<float, int, ORDER> >;
