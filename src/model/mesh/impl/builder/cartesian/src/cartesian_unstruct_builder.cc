/**
 * @file cartesian_unstruct_builder.cc
 * @brief Explicit instantiation of CartesianUnstructBuilder for
 *        ScalarType = int and FloatType = float.
 */
#include "cartesian_unstruct_builder.h"

template class model::CartesianUnstructBuilder<float, int>;
