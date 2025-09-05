#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "model_unstruct.h"

namespace py = pybind11;

// Define enums for template selection
enum class FloatType { Float, Double };
enum class IntType { Int, Long };

PYBIND11_MODULE(model, m) {

    // Create submodule 'model'
    m.attr("__name__") = "pyproxys.model";

    // Bind enums
    py::enum_<FloatType>(m, "FloatType")
        .value("FLOAT", FloatType::Float)
        .value("DOUBLE", FloatType::Double);

    py::enum_<IntType>(m, "IntType")
        .value("INT", IntType::Int)
        .value("LONG", IntType::Long);

}