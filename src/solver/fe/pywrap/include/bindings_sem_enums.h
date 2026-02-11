#ifndef SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_ENUMS_H_
#define SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_ENUMS_H_
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sem_enums.h"

namespace py = pybind11;

namespace solver
{
namespace fe
{

void bind_method_type(py::module_ &m)
{
  py::enum_<enums::methodType>(m, "MethodType")
      .value("SEM", enums::methodType::kSem)
      .value("DG", enums::methodType::kDg);
}

void bind_implem_type(py::module_ &m)
{
  py::enum_<enums::implemType>(m, "ImplemType")
      .value("MAKUTU", enums::implemType::kMakutu);
}

void bind_mesh_type(py::module_ &m)
{
  py::enum_<enums::meshType>(m, "MeshType")
      .value("STRUCT", enums::meshType::kStruct)
      .value("UNSTRUCT", enums::meshType::kUnstruct);
}

void bind_model_location_type(py::module_ &m)
{
  py::enum_<enums::modelLocationType>(m, "ModelLocationType")
      .value("ONNODES", enums::modelLocationType::kOnNodes)
      .value("ONELEMENTS", enums::modelLocationType::kOnElements)
      .export_values();
}

void bind_physic_type(py::module_ &m)
{
  py::enum_<enums::physicType>(m, "PhysicType")
      .value("ACOUSTIC", enums::physicType::kAcoustic)
      .value("ELASTIC", enums::physicType::kElastic)
      .export_values();
}

void bind_all_sem_enums(py::module_ &m)
{
  bind_method_type(m);
  bind_implem_type(m);
  bind_mesh_type(m);
  bind_model_location_type(m);
  bind_physic_type(m);
}

}  // namespace fe
}  // namespace solver
#endif  // SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_ENUMS_H_
