#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cartesian_struct_builder.h"
#include "cartesian_unstruct_builder.h"
#include "model_struct.h"
#include "model_unstruct.h"

namespace py = pybind11;

// Define enums for template selection
enum class FloatType
{
  Float,
  Double
};
enum class IntType
{
  Int,
  Long
};

PYBIND11_MODULE(model, m)
{
  // Create submodule 'model'
  m.attr("__name__") = "pyproxys.model";

  // Bind enums
  py::enum_<FloatType>(m, "FloatType")
      .value("FLOAT", FloatType::Float)
      .value("DOUBLE", FloatType::Double);

  py::enum_<IntType>(m, "IntType")
      .value("INT", IntType::Int)
      .value("LONG", IntType::Long);

  // TODO (see https://github.com/proxySem/proxys/issues/52)
  // until we get a factory that returns a shared pointer we have to bind
  // specific template instances so we only map what is used in pyfwi and python
  // examples

  // Bind ModelApi<float, int>
  using ModelAPIFI = model::ModelApi<float, int>;
  py::class_<ModelAPIFI, std::shared_ptr<ModelAPIFI>>(m, "ModelApi");

  // Bind ModelUnstruct<float, int> (only one used by pyfwi so far)
  using ModelUnstructFI = model::ModelUnstruct<float, int>;
  py::class_<ModelUnstructFI, ModelAPIFI, std::shared_ptr<ModelUnstructFI>>(
      m, "ModelUnstruct")
      .def(py::init<>())
      .def("get_number_of_elements", &ModelUnstructFI::getNumberOfElements)
      .def("get_number_of_nodes", &ModelUnstructFI::getNumberOfNodes)
      .def("get_number_of_points_per_element",
           &ModelUnstructFI::getNumberOfPointsPerElement)
      .def("get_order", &ModelUnstructFI::getOrder)
      .def("domain_size", &ModelUnstructFI::domainSize)
      .def("node_coord", &ModelUnstructFI::nodeCoord)
      .def("global_node_index", &ModelUnstructFI::globalNodeIndex)
      .def("get_model_vp_on_nodes", &ModelUnstructFI::getModelVpOnNodes)
      .def("get_model_vp_on_element", &ModelUnstructFI::getModelVpOnElement)
      .def("get_model_rho_on_nodes", &ModelUnstructFI::getModelRhoOnNodes)
      .def("get_model_rho_on_element", &ModelUnstructFI::getModelRhoOnElement)
      .def("boundary_type", &ModelUnstructFI::boundaryType);

  // Bind ModelStruct<float, int, order> for orders 1, 2, and 3
  using ModelStructFI1 = model::ModelStruct<float, int, 1>;
  py::class_<ModelStructFI1, ModelAPIFI, std::shared_ptr<ModelStructFI1>>(
      m, "ModelStructFI1")
      .def(py::init<>())
      .def(py::init<const model::ModelStructData<float, int>&>())
      .def("get_number_of_elements", &ModelStructFI1::getNumberOfElements)
      .def("get_number_of_nodes", &ModelStructFI1::getNumberOfNodes)
      .def("get_number_of_points_per_element",
           &ModelStructFI1::getNumberOfPointsPerElement)
      .def("get_order", &ModelStructFI1::getOrder)
      .def("domain_size", &ModelStructFI1::domainSize)
      .def("node_coord", &ModelStructFI1::nodeCoord)
      .def("global_node_index", &ModelStructFI1::globalNodeIndex)
      .def("get_model_vp_on_nodes", &ModelStructFI1::getModelVpOnNodes)
      .def("get_model_vp_on_element", &ModelStructFI1::getModelVpOnElement)
      .def("get_model_rho_on_nodes", &ModelStructFI1::getModelRhoOnNodes)
      .def("get_model_rho_on_element", &ModelStructFI1::getModelRhoOnElement)
      .def("boundary_type", &ModelStructFI1::boundaryType);

  using ModelStructFI2 = model::ModelStruct<float, int, 2>;
  py::class_<ModelStructFI2, ModelAPIFI, std::shared_ptr<ModelStructFI2>>(
      m, "ModelStructFI2")
      .def(py::init<>())
      .def(py::init<const model::ModelStructData<float, int>&>())
      .def("get_number_of_elements", &ModelStructFI2::getNumberOfElements)
      .def("get_number_of_nodes", &ModelStructFI2::getNumberOfNodes)
      .def("get_number_of_points_per_element",
           &ModelStructFI2::getNumberOfPointsPerElement)
      .def("get_order", &ModelStructFI2::getOrder)
      .def("domain_size", &ModelStructFI2::domainSize)
      .def("node_coord", &ModelStructFI2::nodeCoord)
      .def("global_node_index", &ModelStructFI2::globalNodeIndex)
      .def("get_model_vp_on_nodes", &ModelStructFI2::getModelVpOnNodes)
      .def("get_model_vp_on_element", &ModelStructFI2::getModelVpOnElement)
      .def("get_model_rho_on_nodes", &ModelStructFI2::getModelRhoOnNodes)
      .def("get_model_rho_on_element", &ModelStructFI2::getModelRhoOnElement)
      .def("boundary_type", &ModelStructFI2::boundaryType);

  using ModelStructFI3 = model::ModelStruct<float, int, 3>;
  py::class_<ModelStructFI3, ModelAPIFI, std::shared_ptr<ModelStructFI3>>(
      m, "ModelStructFI3")
      .def(py::init<>())
      .def(py::init<const model::ModelStructData<float, int>&>())
      .def("get_number_of_elements", &ModelStructFI3::getNumberOfElements)
      .def("get_number_of_nodes", &ModelStructFI3::getNumberOfNodes)
      .def("get_number_of_points_per_element",
           &ModelStructFI3::getNumberOfPointsPerElement)
      .def("get_order", &ModelStructFI3::getOrder)
      .def("domain_size", &ModelStructFI3::domainSize)
      .def("node_coord", &ModelStructFI3::nodeCoord)
      .def("global_node_index", &ModelStructFI3::globalNodeIndex)
      .def("get_model_vp_on_nodes", &ModelStructFI3::getModelVpOnNodes)
      .def("get_model_vp_on_element", &ModelStructFI3::getModelVpOnElement)
      .def("get_model_rho_on_nodes", &ModelStructFI3::getModelRhoOnNodes)
      .def("get_model_rho_on_element", &ModelStructFI3::getModelRhoOnElement)
      .def("boundary_type", &ModelStructFI3::boundaryType);

  // Bind ModelStructData<float, int>
  py::class_<model::ModelStructData<float, int>>(m, "ModelStructData")
      .def(py::init<>())
      .def_readwrite("ex", &model::ModelStructData<float, int>::ex_)
      .def_readwrite("ey", &model::ModelStructData<float, int>::ey_)
      .def_readwrite("ez", &model::ModelStructData<float, int>::ez_)
      .def_readwrite("dx", &model::ModelStructData<float, int>::dx_)
      .def_readwrite("dy", &model::ModelStructData<float, int>::dy_)
      .def_readwrite("dz", &model::ModelStructData<float, int>::dz_)
      .def_readwrite("order", &model::ModelStructData<float, int>::order_);

  // Bind CartesianStructBuilder<float, int, order> for orders 1, 2, and 3
  using CartesianStructBuilderFI1 =
      model::CartesianStructBuilder<float, int, 1>;
  py::class_<CartesianStructBuilderFI1,
             std::shared_ptr<CartesianStructBuilderFI1>>(
      m, "CartesianStructBuilderFI1")
      .def(py::init<int, float, int, float, int, float>(), py::arg("ex"),
           py::arg("hx"), py::arg("ey"), py::arg("hy"), py::arg("ez"),
           py::arg("hz"), py::arg("isModelOnNodes"))
      .def("get_model", &CartesianStructBuilderFI1::getModel);

  using CartesianStructBuilderFI2 =
      model::CartesianStructBuilder<float, int, 2>;
  py::class_<CartesianStructBuilderFI2,
             std::shared_ptr<CartesianStructBuilderFI2>>(
      m, "CartesianStructBuilderFI2")
      .def(py::init<int, float, int, float, int, float>(), py::arg("ex"),
           py::arg("hx"), py::arg("ey"), py::arg("hy"), py::arg("ez"),
           py::arg("hz"), py::arg("isModelOnNodes"))
      .def("get_model", &CartesianStructBuilderFI2::getModel);

  using CartesianStructBuilderFI3 =
      model::CartesianStructBuilder<float, int, 3>;
  py::class_<CartesianStructBuilderFI3,
             std::shared_ptr<CartesianStructBuilderFI3>>(
      m, "CartesianStructBuilderFI3")
      .def(py::init<int, float, int, float, int, float>(), py::arg("ex"),
           py::arg("hx"), py::arg("ey"), py::arg("hy"), py::arg("ez"),
           py::arg("hz"), py::arg("isModelOnNodes"))
      .def("get_model", &CartesianStructBuilderFI3::getModel);

  // Bind CartesianParams<float, int>
  using CartesianParamsFI = model::CartesianParams<float, int>;
  py::class_<CartesianParamsFI, std::shared_ptr<CartesianParamsFI>>(
      m, "CartesianParams")
      .def(py::init<>())
      .def(py::init<int, int, int, int, float, float, float>(),
           py::arg("order"), py::arg("ex"), py::arg("ey"), py::arg("ez"),
           py::arg("lx"), py::arg("ly"), py::arg("lz"))
      .def_readwrite("order", &CartesianParamsFI::order)
      .def_readwrite("ex", &CartesianParamsFI::ex)
      .def_readwrite("ey", &CartesianParamsFI::ey)
      .def_readwrite("ez", &CartesianParamsFI::ez)
      .def_readwrite("lx", &CartesianParamsFI::lx)
      .def_readwrite("ly", &CartesianParamsFI::ly)
      .def_readwrite("lz", &CartesianParamsFI::lz)
      .def_readwrite("isModelOnNodes", &CartesianParamsFI::isModelOnNodes);

  // Bind CartesianUnstructBuilder<float, int>
  using CartesianUnstructBuilderFI =
      model::CartesianUnstructBuilder<float, int>;
  py::class_<CartesianUnstructBuilderFI,
             std::shared_ptr<CartesianUnstructBuilderFI>>(
      m, "CartesianUnstructBuilder")
      .def(py::init<>())
      .def(py::init<const model::CartesianParams<float, int>&>())
      .def("get_model", &CartesianUnstructBuilderFI::getModel);
}
