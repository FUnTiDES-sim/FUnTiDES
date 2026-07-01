#ifndef FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_DIFFERENTIATOR_H_
#define FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_DIFFERENTIATOR_H_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <Kokkos_Core.hpp>
#include <memory>

#include "data_type.h"
#include "differentiator.h"
#include "differentiator_data_acoustic.h"
#include "differentiator_data_elastic.h"
#include "differentiator_factory.h"
#include "gradient_acoustic.h"
#include "gradient_elastic.h"
#include "wavefield_view_backward_acoustic.h"
#include "wavefield_view_backward_elastic.h"
#include "wavefield_view_forward_acoustic.h"
#include "wavefield_view_forward_elastic.h"

namespace py = pybind11;

namespace gradient {

void bind_data_struct(py::module_& m) {
  py::class_<Differentiator::DataStruct, std::shared_ptr<Differentiator::DataStruct>>(m, "DataStruct")
      .def("print", &Differentiator::DataStruct::print);
}

void bind_gradient_data_acoustic(py::module_& m) {
  py::class_<GradientDataAcoustic, Differentiator::DataStruct, std::shared_ptr<GradientDataAcoustic>>(
      m, "GradientDataAcoustic")
      .def(py::init<const WavefieldViewForwardAcoustic&, const WavefieldViewBackwardAcoustic&,
                    const GradientAcoustic&>(),
           py::arg("fwd"), py::arg("bwd"), py::arg("gradient"))
      .def("print", &GradientDataAcoustic::print);
}

void bind_gradient_data_elastic(py::module_& m) {
  py::class_<GradientDataElastic, Differentiator::DataStruct, std::shared_ptr<GradientDataElastic>>(
      m, "GradientDataElastic")
      .def(py::init<const WavefieldViewForwardElastic&, const WavefieldViewBackwardElastic&, const GradientElastic&>(),
           py::arg("fwd"), py::arg("bwd"), py::arg("gradient"))
      .def("print", &GradientDataElastic::print);
}

void bind_differentiator_base(py::module_& m) {
  py::class_<Differentiator, std::shared_ptr<Differentiator>>(m, "Differentiator")
      .def("compute", &Differentiator::compute, py::arg("mesh"), py::arg("data"), py::arg("dt"))
      .def("init_geometric_mass_matrix", &Differentiator::initGeometricMassMatrix, py::arg("mesh"))
      .def(
          "get_geometric_mass_matrix",
          [](Differentiator& self) -> Kokkos::Experimental::python_view_type_t<vectorReal> {
            return self.getGeometricMassMatrix();
          },
          py::return_value_policy::reference_internal)
      .def("print", &Differentiator::print);
}

void bind_differentiator_factory(py::module_& m) {
  m.def(
      "create_differentiator",
      [](utils::enums::implemType implem, utils::enums::meshType mesh, utils::enums::modelLocationType modelLocation,
         utils::enums::physicType physic, int order) {
        auto diff = createDifferentiator(implem, mesh, modelLocation, physic, order);
        return std::shared_ptr<Differentiator>(std::move(diff));
      },
      py::arg("implem_type"), py::arg("mesh_type"), py::arg("model_location"), py::arg("physic_type"),
      py::arg("order"));
}

}  // namespace gradient

#endif  // FUNTIDES_GRADIENT_PYWRAP_INCLUDE_BINDINGS_DIFFERENTIATOR_H_