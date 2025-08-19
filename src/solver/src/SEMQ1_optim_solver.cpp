#include "SEMsolver_impl.hpp"

constexpr int ORDER = 1;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::OPTIM >::type, CartesianSEMmesh<float, int, ORDER> >;
