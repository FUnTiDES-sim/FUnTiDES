#ifndef FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_DEVICE_UTILS_H_
#define FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_DEVICE_UTILS_H_
#include <pybind11/pybind11.h>

#include <KokkosExp_InterOp.hpp>

#include "data_type.h"

namespace py = pybind11;

namespace solver {
namespace fe {

/**
 * @brief Scatter/gather copy between two per-element DG fields, restricted to listed dofs.
 *
 * Device-side equivalent of the numpy fancy-indexing ghost-layer exchange
 * dst[np.ix_(dst_elems, dst_dofs)] = src[np.ix_(src_elems, src_dofs)] used by the
 * Python-driven multi-solver pipelines: keeps the per-time-step face exchange on the
 * device instead of pulling UVM pages to the host every step.
 *
 * @param dst       Destination field (nElem x nDofPerElem).
 * @param dst_elems Destination element indices (one per exchanged element).
 * @param src       Source field (nElem x nDofPerElem).
 * @param src_elems Source element indices (same length as @p dst_elems).
 * @param dst_dofs  Destination local-dof indices (one per exchanged face dof).
 * @param src_dofs  Source local-dof indices (same length as @p dst_dofs).
 */
inline void copyFaceDofsImpl(arrayReal dst, vectorInt dst_elems, arrayReal src, vectorInt src_elems,
                             vectorInt dst_dofs, vectorInt src_dofs) {
  int const n_d = static_cast<int>(dst_dofs.extent(0));
  int const n = static_cast<int>(dst_elems.extent(0)) * n_d;
  Kokkos::parallel_for(
      "PyCopyFaceDofs", n, KOKKOS_LAMBDA(const int idx) {
        int const e = idx / n_d;
        int const d = idx % n_d;
        dst(dst_elems(e), dst_dofs(d)) = src(src_elems(e), src_dofs(d));
      });
}

/**
 * @brief Gather field values at listed nodes into one column of a 2D trace buffer.
 *
 * out(i, col) = field(nodes(i)) — device-side receiver-trace extraction: the trace
 * accumulates on the device over the whole run and is read back to the host once at
 * the end, instead of one numpy/UVM read per time step.
 *
 * @param out   Trace buffer (nReceivers x nSamples).
 * @param col   Time-sample column to fill.
 * @param field Nodal field to sample (e.g. a SEM pressure field).
 * @param nodes Global node index per receiver.
 */
inline void gatherColumnImpl(arrayReal out, int col, vectorReal field, vectorInt nodes) {
  int const n = static_cast<int>(nodes.extent(0));
  Kokkos::parallel_for(
      "PyGatherColumn", n, KOKKOS_LAMBDA(const int i) { out(i, col) = field(nodes(i)); });
}

/**
 * @brief Max of |field| over a per-element 2D field, computed on the device.
 * @param field Field to reduce (nElem x nDofPerElem).
 * @return max(|field|).
 */
inline float maxAbsArrayImpl(arrayReal field) {
  int const n1 = static_cast<int>(field.extent(1));
  int const n = static_cast<int>(field.extent(0)) * n1;
  float result = 0.0f;
  Kokkos::parallel_reduce(
      "PyMaxAbs2D", n,
      KOKKOS_LAMBDA(const int idx, float &lmax) {
        float const v = Kokkos::fabs(field(idx / n1, idx % n1));
        if (v > lmax) lmax = v;
      },
      Kokkos::Max<float>(result));
  return result;
}

/**
 * @brief Max of |field| over a nodal 1D field, computed on the device.
 * @param field Field to reduce (nNodes).
 * @return max(|field|).
 */
inline float maxAbsVectorImpl(vectorReal field) {
  int const n = static_cast<int>(field.extent(0));
  float result = 0.0f;
  Kokkos::parallel_reduce(
      "PyMaxAbs1D", n,
      KOKKOS_LAMBDA(const int i, float &lmax) {
        float const v = Kokkos::fabs(field(i));
        if (v > lmax) lmax = v;
      },
      Kokkos::Max<float>(result));
  return result;
}

/**
 * @brief Register the device-side utility functions for Python-driven pipelines.
 *
 * These helpers let driver scripts keep every per-time-step data motion on the device
 * (ghost-layer exchange, receiver-trace extraction, |p|_max diagnostics), so wavefield
 * buffers are never touched through numpy/UVM inside the time loop.
 */
inline void bind_device_utils(py::module_ &m) {
  m.def(
      "copy_face_dofs",
      [](Kokkos::Experimental::python_view_type_t<arrayReal> dst,
         Kokkos::Experimental::python_view_type_t<vectorInt> dst_elems,
         Kokkos::Experimental::python_view_type_t<arrayReal> src,
         Kokkos::Experimental::python_view_type_t<vectorInt> src_elems,
         Kokkos::Experimental::python_view_type_t<vectorInt> dst_dofs,
         Kokkos::Experimental::python_view_type_t<vectorInt> src_dofs) {
        copyFaceDofsImpl(dst, dst_elems, src, src_elems, dst_dofs, src_dofs);
      },
      py::arg("dst"), py::arg("dst_elems"), py::arg("src"), py::arg("src_elems"), py::arg("dst_dofs"),
      py::arg("src_dofs"));

  m.def(
      "gather_column",
      [](Kokkos::Experimental::python_view_type_t<arrayReal> out, int col,
         Kokkos::Experimental::python_view_type_t<vectorReal> field,
         Kokkos::Experimental::python_view_type_t<vectorInt> nodes) { gatherColumnImpl(out, col, field, nodes); },
      py::arg("out"), py::arg("col"), py::arg("field"), py::arg("nodes"));

  m.def(
      "max_abs",
      [](Kokkos::Experimental::python_view_type_t<arrayReal> field) { return maxAbsArrayImpl(field); },
      py::arg("field"));
  m.def(
      "max_abs",
      [](Kokkos::Experimental::python_view_type_t<vectorReal> field) { return maxAbsVectorImpl(field); },
      py::arg("field"));

  m.def("device_fence", []() { Kokkos::fence(); });
}

}  // namespace fe
}  // namespace solver
#endif  // FUNTIDES_SOLVER_FE_PYWRAP_INCLUDE_BINDINGS_DEVICE_UTILS_H_
