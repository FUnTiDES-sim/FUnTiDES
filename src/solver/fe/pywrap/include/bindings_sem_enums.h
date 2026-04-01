#ifndef FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_ENUMS_H_
#define FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_ENUMS_H_
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
  py::enum_<utils::enums::methodType>(m, "MethodType")
      .value("SEM", utils::enums::methodType::kSem)
      .value("DG", utils::enums::methodType::kDg);
}

void bind_implem_type(py::module_ &m)
{
  py::enum_<utils::enums::implemType>(m, "ImplemType")
      .value("MAKUTU", utils::enums::implemType::kMakutu);
}

void bind_mesh_type(py::module_ &m)
{
  py::enum_<utils::enums::meshType>(m, "MeshType")
      .value("STRUCT", utils::enums::meshType::kStruct)
      .value("UNSTRUCT", utils::enums::meshType::kUnstruct);
}

void bind_model_location_type(py::module_ &m)
{
  py::enum_<utils::enums::modelLocationType>(m, "ModelLocationType")
      .value("ONNODES", utils::enums::modelLocationType::kOnNodes)
      .value("ONELEMENTS", utils::enums::modelLocationType::kOnElements)
      .export_values();
}

void bind_physic_type(py::module_ &m)
{
  py::enum_<utils::enums::physicType>(m, "PhysicType")
      .value("ACOUSTIC", utils::enums::physicType::kAcoustic)
      .value("ELASTIC", utils::enums::physicType::kElastic)
      .value("ACOUSTOELASTIC", utils::enums::physicType::kAcoustoElastic)
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
#endif  // FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_SEM_ENUMS_H_
