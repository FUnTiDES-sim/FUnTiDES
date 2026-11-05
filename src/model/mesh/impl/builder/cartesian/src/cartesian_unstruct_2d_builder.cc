/**
 * @file cartesian_unstruct_2d_builder.cc
 * @brief Explicit instantiation of CartesianUnstruct2DBuilder for
 *        ScalarType = int and FloatType = float (built with COMPILE_SEM_2D).
 */
#include "cartesian_unstruct_2d_builder.h"

template class model::CartesianUnstruct2DBuilder<float, int>;
