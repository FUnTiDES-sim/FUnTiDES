#ifndef FUNTIDES_IO_PYWRAP_INCLUDE_BINDINGS_IO_H_
#define FUNTIDES_IO_PYWRAP_INCLUDE_BINDINGS_IO_H_

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "io_controller_base.h"

namespace py = pybind11;

namespace io_bindings {

/// Scalar type of a snapshot, as stored by every backend.
using Real = funtides::io::HostVectorReal::non_const_value_type;

/// 1D C-contiguous NumPy array of Real. Without forcecast: a wrong dtype is
/// rejected instead of being silently converted into a temporary copy.
using RealArray = py::array_t<Real, py::array::c_style>;

/**
 * @brief Wraps a NumPy buffer in an unmanaged host view, without copy.
 *
 * The view does not own the memory: the caller keeps `arr` alive for as long
 * as the backend may read the view.
 */
inline funtides::io::HostVectorReal asHostView(const RealArray& arr) {
  if (arr.ndim() != 1) {
    throw py::value_error("funtides.io: a snapshot must be a 1D array, got " + std::to_string(arr.ndim()) +
                          " dimensions");
  }
  // The view type has a non-const value type; write paths never modify it.
  return funtides::io::HostVectorReal(const_cast<Real*>(arr.data()), static_cast<std::size_t>(arr.shape(0)));
}

/**
 * @brief Python-side owner of an IOControllerBase.
 *
 * Two duties the C++ interface leaves to its caller, and that a Python caller
 * cannot fulfil by itself:
 * - BUFFER CONTRACT. writeSnapshot() may read `field` after it returns, until
 *   the next writeSnapshot(), flush() or close(). In Python, nothing keeps the
 *   array alive in the meantime, so this class holds a reference to the last
 *   written array and drops it only once the contract releases it.
 * - GIL. Every call that may block on I/O (write back-pressure, read, flush,
 *   close) releases the GIL, so other Python threads keep running. The
 *   controller never calls back into Python, which makes this safe.
 */
class PyIOController {
 public:
  explicit PyIOController(std::unique_ptr<funtides::io::IOControllerBase> ctrl) : ctrl_(std::move(ctrl)) {}

  void writeSnapshot(const RealArray& field) {
    auto& ctrl = get();
    const funtides::io::HostVectorReal view = asHostView(field);
    {
      py::gil_scoped_release release;
      ctrl.writeSnapshot(view);
    }
    // The previous array is released by this call, the new one is held until the next.
    pinned_ = field;
  }

  void readSnapshot(RealArray& field, std::size_t index) {
    auto& ctrl = get();
    if (!field.writeable()) throw py::value_error("funtides.io: read_snapshot needs a writeable array");
    const funtides::io::HostVectorReal view = asHostView(field);
    py::gil_scoped_release release;
    ctrl.readSnapshot(view, index);
  }

  void flush() {
    auto& ctrl = get();
    {
      py::gil_scoped_release release;
      ctrl.flush();
    }
    pinned_ = py::none();
  }

  /// Idempotent, like IOControllerBase::close(). The C++ object is kept: its
  /// own methods report "after close" errors with the usual messages.
  void close() {
    if (closed_) return;
    closed_ = true;
    try {
      py::gil_scoped_release release;
      ctrl_->close();
    } catch (...) {
      pinned_ = py::none();
      throw;
    }
    pinned_ = py::none();
  }

  bool closed() const { return closed_; }

 private:
  funtides::io::IOControllerBase& get() {
    if (!ctrl_) throw std::runtime_error("funtides.io: controller is not initialized");
    return *ctrl_;
  }

  std::unique_ptr<funtides::io::IOControllerBase> ctrl_;
  py::object pinned_ = py::none();  ///< Last array handed to writeSnapshot().
  bool closed_{false};
};

/**
 * @brief Binds BackendKind and OpenMode, with their values exported to the
 *        module scope.
 * @param[in,out] m Python module receiving the enums.
 */
inline void bind_io_enums(py::module_& m) {
  py::enum_<funtides::io::BackendKind>(m, "BackendKind")
      .value("POSIX", funtides::io::BackendKind::kPosix)
      .value("ADIOS2", funtides::io::BackendKind::kAdios2)
      .export_values();

  py::enum_<funtides::io::OpenMode>(m, "OpenMode")
      .value("WRITE", funtides::io::OpenMode::kWrite)
      .value("READ", funtides::io::OpenMode::kRead)
      .export_values();
}

/**
 * @brief Binds IOConfig with snake_case read/write attributes.
 *
 * The std::vector fields are converted by value (pybind11/stl.h): assign a
 * whole list, `cfg.local_dims = [4, 5, 6]`; `cfg.local_dims.append(7)`
 * modifies a temporary copy and is lost.
 * Of DistributedContext, only the rank is read by the backends; it is exposed
 * directly as `rank`.
 * @param[in,out] m Python module receiving the class.
 */
inline void bind_io_config(py::module_& m) {
  using funtides::io::IOConfig;

  py::class_<IOConfig>(m, "IOConfig")
      .def(py::init<>())
      .def_property(
          "rank", [](const IOConfig& c) { return c.ctx.rank; },
          [](IOConfig& c, decltype(std::declval<IOConfig&>().ctx.rank) rank) { c.ctx.rank = rank; },
          "Rank of this process; part of every file name.")
      .def_readwrite("output_dir", &IOConfig::output_dir)
      .def_readwrite("prefix", &IOConfig::prefix)
      .def_readwrite("shot_id", &IOConfig::shot_id)
      .def_readwrite("global_dims", &IOConfig::global_dims)
      .def_readwrite("start_offsets", &IOConfig::start_offsets)
      .def_readwrite("local_dims", &IOConfig::local_dims)
      .def_readwrite("nt", &IOConfig::nt)
      .def_readwrite("nb_receiver", &IOConfig::nb_receiver)
      .def_readwrite("async_snapshots", &IOConfig::async_snapshots)
      .def("__repr__", [](const IOConfig& c) {
        return "IOConfig(rank=" + std::to_string(c.ctx.rank) + ", output_dir='" + c.output_dir + "', prefix='" +
               c.prefix + "', shot_id='" + c.shot_id + "', async_snapshots=" + (c.async_snapshots ? "True" : "False") +
               ")";
      });
}

/**
 * @brief Binds the I/O controller and its factory.
 *
 * Both backends are reached through the same class and the factory, as in
 * C++; the concrete backend classes are not exposed. `IOController` is also a
 * context manager that calls close() on exit, because the destructor runs at
 * a time Python does not guarantee and can only print close errors.
 * @param[in,out] m Python module receiving the class and functions.
 */
inline void bind_io_controller(py::module_& m) {
  py::class_<PyIOController>(m, "IOController")
      .def("write_snapshot", &PyIOController::writeSnapshot, py::arg("field"),
           "Writes the next snapshot from a 1D C-contiguous array of the solver's real type.")
      .def("read_snapshot", &PyIOController::readSnapshot, py::arg("field").noconvert(), py::arg("index"),
           "Reads snapshot `index` in place into `field` (1D, C-contiguous, writeable, exact dtype).")
      .def("flush", &PyIOController::flush)
      .def("close", &PyIOController::close)
      .def_property_readonly("closed", &PyIOController::closed)
      .def(
          "__enter__", [](PyIOController& self) -> PyIOController& { return self; },
          py::return_value_policy::reference_internal)
      .def("__exit__", [](PyIOController& self, const py::object&, const py::object&, const py::object&) {
        self.close();
        return false;  // never swallow the exception of the with block
      });

  m.def(
      "make_io_controller",
      [](funtides::io::BackendKind kind, funtides::io::OpenMode mode, const funtides::io::IOConfig& config) {
        return PyIOController(funtides::io::makeIOController(kind, mode, config));
      },
      py::arg("kind"), py::arg("mode"), py::arg("config"),
      "Creates an I/O controller. Raises RuntimeError if the backend is not built in.");

  m.def("is_backend_available", &funtides::io::isBackendAvailable, py::arg("kind"));

  // Exposed so Python code can allocate snapshot buffers with the right dtype.
  m.attr("real_dtype") = py::dtype::of<Real>();
}

}  // namespace io_bindings

#endif  // FUNTIDES_IO_PYWRAP_INCLUDE_BINDINGS_IO_H_
