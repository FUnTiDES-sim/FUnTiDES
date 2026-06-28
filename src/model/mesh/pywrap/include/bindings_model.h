#ifndef FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_MODEL_H_
#define FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_MODEL_H_

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>

#include "bindings_utils.h"
#include "common_macros.h"
#include "face_connectivity_unstruct.h"
#include "model.h"
#include "model_struct.h"
#include "model_unstruct.h"

namespace py = pybind11;

namespace model {

// template binder for ModelAPI
template <typename FloatType, typename ScalarType>
void bind_modelapi(py::module_ &m) {
  using T = model::ModelApi<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("ModelApi");

  py::class_<T, std::shared_ptr<T>>(m, name.c_str())
      .def("node_coord", &T::nodeCoord)
      .def("global_node_index", &T::globalNodeIndex)
      .def("get_model_vp_on_node", &T::getModelVpOnNodes)
      .def("get_model_vp_on_element", &T::getModelVpOnElement)
      .def("get_model_rho_on_node", &T::getModelRhoOnNodes)
      .def("get_model_rho_on_element", &T::getModelRhoOnElement)
      .def("get_model_vs_on_node", &T::getModelVsOnNodes)
      .def("get_model_vs_on_element", &T::getModelVsOnElement)
      .def("get_model_qp_on_node", &T::getModelQpOnNodes)
      .def("get_model_qp_on_element", &T::getModelQpOnElement)
      .def("get_model_qs_on_node", &T::getModelQsOnNodes)
      .def("get_model_qs_on_element", &T::getModelQsOnElement)
      .def("get_model_delta_on_node", &T::getModelDeltaOnNodes)
      .def("get_model_delta_on_element", &T::getModelDeltaOnElement)
      .def("get_model_epsilon_on_node", &T::getModelEpsilonOnNodes)
      .def("get_model_epsilon_on_element", &T::getModelEpsilonOnElement)
      .def("get_model_gamma_on_node", &T::getModelGammaOnNodes)
      .def("get_model_gamma_on_element", &T::getModelGammaOnElement)
      .def("get_model_phi_on_node", &T::getModelPhiOnNodes)
      .def("get_model_phi_on_element", &T::getModelPhiOnElement)
      .def("get_model_theta_on_node", &T::getModelThetaOnNodes)
      .def("get_model_theta_on_element", &T::getModelThetaOnElement)
      .def("get_number_of_elements", &T::getNumberOfElements)
      .def("get_number_of_nodes", &T::getNumberOfNodes)
      .def("get_number_of_points_per_element", &T::getNumberOfPointsPerElement)
      .def("get_order", &T::getOrder)
      .def("boundary_type", &T::boundaryType)
      .def("face_normal", &T::faceNormal)
      .def("domain_size", &T::domainSize)
      .def("get_min_spacing", &T::getMinSpacing)
      .def("get_max_speed", &T::getMaxSpeed)
      // Build the face connectivity and register a Kokkos finalize hook that
      // releases its internally-owned device Views before Kokkos::finalize().
      // The client (pykokkos) may call kokkos.finalize() before this model is
      // destroyed; without this hook, the managed Views would be deallocated
      // after finalize and trigger a Kokkos host_abort. The hook captures a
      // weak_ptr so it is a no-op if the model is already destroyed.
      .def("build_face_connectivity",
           [](std::shared_ptr<T> self) {
             self->buildFaceConnectivity();
             std::weak_ptr<T> weak = self;
             Kokkos::push_finalize_hook([weak]() {
               if (auto locked = weak.lock()) locked->releaseFaceConnectivity();
             });
           })
      .def("get_number_of_faces", &T::getNumberOfFaces)
      .def("is_boundary_face", &T::isBoundaryFace)
      .def("get_global_node_from_face", &T::getGlobalNodeFromFace)
      .def("get_global_face", &T::getGlobalFace)
      .def("set_quality_factors", &T::setQualityFactors, py::arg("qp"), py::arg("qs"))
      .def("is_free_surface", &T::isFreeSurface);
}

// templated binder for one ModelStruct instantiation
template <typename FloatType, typename ScalarType, int Order>
void bind_modelstruct(py::module_ &m) {
  using Base = model::ModelApi<FloatType, ScalarType>;
  using T = model::ModelStruct<FloatType, ScalarType, Order>;
  using Data = model::ModelStructData<FloatType, ScalarType>;

  std::string name = model_class_name<FloatType, ScalarType, Order>("ModelStruct");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str()).def(py::init<const Data &>());
}

// templated binder for ModelStructData
template <typename FloatType, typename ScalarType>
void bind_modelstructdata(py::module_ &m) {
  using Data = model::ModelStructData<FloatType, ScalarType>;
  std::string name = model_class_name<FloatType, ScalarType>("ModelStructData");

  py::class_<Data>(m, name.c_str())
      .def(py::init<>())
      .def_readwrite("ex", &Data::ex_)
      .def_readwrite("ey", &Data::ey_)
      .def_readwrite("ez", &Data::ez_)
      .def_readwrite("dx", &Data::dx_)
      .def_readwrite("dy", &Data::dy_)
      .def_readwrite("dz", &Data::dz_);
}

// templated binder for ModelUnstruct
template <typename FloatType, typename ScalarType>
void bind_modelunstruct(py::module_ &m) {
  using Base = model::ModelApi<FloatType, ScalarType>;
  using T = model::ModelUnstruct<FloatType, ScalarType>;
  using Data = model::ModelUnstructData<FloatType, ScalarType>;

  std::string name = model_class_name<FloatType, ScalarType>("ModelUnstruct");

  py::class_<T, Base, std::shared_ptr<T>>(m, name.c_str()).def(py::init<const Data &>());
}

// template binder for ModelUnstructData
template <typename FloatType, typename ScalarType>
void bind_modelunstructdata(py::module_ &m) {
  using Data = model::ModelUnstructData<FloatType, ScalarType>;

  std::string name = model_class_name<FloatType, ScalarType>("ModelUnstructData");

  using PArrInt = Kokkos::Experimental::python_view_type_t<arrayInt>;
  using PVecReal = Kokkos::Experimental::python_view_type_t<vectorReal>;
  using PArr3Real = Kokkos::Experimental::python_view_type_t<array3DReal>;
  using PVecInt = Kokkos::Experimental::python_view_type_t<vectorInt>;

  py::class_<Data, std::shared_ptr<Data>>(m, name.c_str())
      // Constructeur existant INCHANGÉ (sans face_connectivity). Construit via
      // une fabrique afin d'enregistrer un finalize hook : le client détruit ce
      // model_data APRÈS kokkos.finalize(), donc on relâche ici (via weak_ptr)
      // les Views Kokkos gérées avant finalize pour éviter un host_abort.
      .def(py::init([](ScalarType order, ScalarType n_element, ScalarType n_node, FloatType lx, FloatType ly,
                       FloatType lz, bool is_model_on_nodes, bool is_elastic, PArrInt global_node_index,
                       PVecReal nodes_coords_x, PVecReal nodes_coords_y, PVecReal nodes_coords_z, PVecReal model_vp_node,
                       PVecReal model_vp_element, PVecReal model_rho_node, PVecReal model_rho_element,
                       PVecReal model_vs_node, PVecReal model_vs_element, PVecReal model_delta_node,
                       PVecReal model_delta_element, PVecReal model_epsilon_node, PVecReal model_epsilon_element,
                       PVecReal model_gamma_node, PVecReal model_gamma_element, PVecReal model_theta_node,
                       PVecReal model_theta_element, PVecReal model_phi_node, PVecReal model_phi_element,
                       PArr3Real model_C_tensor_element, PVecInt boundaries_t) {
             auto data = std::make_shared<Data>(
                 order, n_element, n_node, lx, ly, lz, is_model_on_nodes, is_elastic, global_node_index, nodes_coords_x,
                 nodes_coords_y, nodes_coords_z, model_vp_node, model_vp_element, model_rho_node, model_rho_element,
                 model_vs_node, model_vs_element, model_delta_node, model_delta_element, model_epsilon_node,
                 model_epsilon_element, model_gamma_node, model_gamma_element, model_theta_node, model_theta_element,
                 model_phi_node, model_phi_element, model_C_tensor_element, boundaries_t);
             std::weak_ptr<Data> weak = data;
             Kokkos::push_finalize_hook([weak]() {
               if (auto locked = weak.lock()) locked->release();
             });
             return data;
           }),
           py::arg("order"), py::arg("n_element"), py::arg("n_node"), py::arg("lx"), py::arg("ly"), py::arg("lz"),
           py::arg("is_model_on_nodes"), py::arg("is_elastic"), py::arg("global_node_index"), py::arg("nodes_coords_x"),
           py::arg("nodes_coords_y"), py::arg("nodes_coords_z"), py::arg("model_vp_node"), py::arg("model_vp_element"),
           py::arg("model_rho_node"), py::arg("model_rho_element"), py::arg("model_vs_node"),
           py::arg("model_vs_element"), py::arg("model_delta_node"), py::arg("model_delta_element"),
           py::arg("model_epsilon_node"), py::arg("model_epsilon_element"), py::arg("model_gamma_node"),
           py::arg("model_gamma_element"), py::arg("model_theta_node"), py::arg("model_theta_element"),
           py::arg("model_phi_node"), py::arg("model_phi_element"), py::arg("model_C_tensor_element"),
           py::arg("boundaries_t"))

      .def_readwrite("face_connectivity", &Data::face_connectivity_);
}

}  // namespace model

#endif  // FUNTIDES_MODEL_MESH_PYWRAP_INCLUDE_BINDINGS_MODEL_H_
