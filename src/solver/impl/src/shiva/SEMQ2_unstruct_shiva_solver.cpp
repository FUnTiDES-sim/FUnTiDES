#include <model_unstruct.h>
#include "SEMsolver_impl.hpp"

constexpr int ORDER = 2;
template class SEMsolver< ORDER, IntegralTypeSelector< ORDER, IntegralType::SHIVA>::type, model::ModelUnstruct<float, int> >;
