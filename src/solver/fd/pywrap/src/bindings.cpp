#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <KokkosExp_InterOp.hpp>

#include <memory>
#include <string>

#include "data_type.h"
#include "fd_abckernels.h"
#include "fd_grids.h"
#include "fd_kernels.h"
#include "fd_options.h"
#include "fd_solver.h"
#include "fd_source_receivers.h"
#include "fd_stencils.h"

namespace py = pybind11;

PYBIND11_MODULE(fd_solver, m)
{
  // Create submodule 'fd_solver'
  m.attr("__name__") = "pyfuntides.fd_solver";
  m.doc() = "Python bindings for FUnTiDES Finite Difference (FD) solver";

  // ============================================================================
  // FDTD Options and Configuration
  // ============================================================================

  // Bind nested structs for FdtdOptions
  py::class_<fdtd::options::FdtdOptions::GridParams>(m, "GridParams")
      .def(py::init<>())
      .def_readwrite("nx", &fdtd::options::FdtdOptions::GridParams::nx,
                     "Grid size in x direction")
      .def_readwrite("ny", &fdtd::options::FdtdOptions::GridParams::ny,
                     "Grid size in y direction")
      .def_readwrite("nz", &fdtd::options::FdtdOptions::GridParams::nz,
                     "Grid size in z direction")
      .def_readwrite("dx", &fdtd::options::FdtdOptions::GridParams::dx,
                     "Grid spacing in x direction (m)")
      .def_readwrite("dy", &fdtd::options::FdtdOptions::GridParams::dy,
                     "Grid spacing in y direction (m)")
      .def_readwrite("dz", &fdtd::options::FdtdOptions::GridParams::dz,
                     "Grid spacing in z direction (m)")
      .def_readwrite("mesh", &fdtd::options::FdtdOptions::GridParams::mesh,
                     "Mesh type (cartesian)");

  py::class_<fdtd::options::FdtdOptions::StencilParams>(m, "StencilParams")
      .def(py::init<>())
      .def_readwrite("lx", &fdtd::options::FdtdOptions::StencilParams::lx,
                     "Half stencil length in x direction")
      .def_readwrite("ly", &fdtd::options::FdtdOptions::StencilParams::ly,
                     "Half stencil length in y direction")
      .def_readwrite("lz", &fdtd::options::FdtdOptions::StencilParams::lz,
                     "Half stencil length in z direction")
      .def_readwrite("implem", &fdtd::options::FdtdOptions::StencilParams::implem,
                     "Stencil implementation (remez|taylor)");

  py::class_<fdtd::options::FdtdOptions::SourceParams>(m, "SourceParams")
      .def(py::init<>())
      .def_readwrite("xs", &fdtd::options::FdtdOptions::SourceParams::xs,
                     "Source x position (grid index)")
      .def_readwrite("ys", &fdtd::options::FdtdOptions::SourceParams::ys,
                     "Source y position (grid index)")
      .def_readwrite("zs", &fdtd::options::FdtdOptions::SourceParams::zs,
                     "Source z position (grid index)")
      .def_readwrite("f0", &fdtd::options::FdtdOptions::SourceParams::f0,
                     "Source frequency (Hz)")
      .def_readwrite("source_order",
                     &fdtd::options::FdtdOptions::SourceParams::source_order,
                     "Source derivative order");

  py::class_<fdtd::options::FdtdOptions::VelocityParams>(m, "VelocityParams")
      .def(py::init<>())
      .def_readwrite("vmin", &fdtd::options::FdtdOptions::VelocityParams::vmin,
                     "Minimum velocity (m/s)")
      .def_readwrite("vmax", &fdtd::options::FdtdOptions::VelocityParams::vmax,
                     "Maximum velocity (m/s)")
      .def_readwrite("file_model",
                     &fdtd::options::FdtdOptions::VelocityParams::file_model,
                     "Path to velocity model file")
      .def_readwrite("use_file_model",
                     &fdtd::options::FdtdOptions::VelocityParams::use_file_model,
                     "Use file-based model");

  py::class_<fdtd::options::FdtdOptions::TimeParams>(m, "TimeParams")
      .def(py::init<>())
      .def_readwrite("time_step",
                     &fdtd::options::FdtdOptions::TimeParams::time_step,
                     "Time step size (s)")
      .def_readwrite("time_max", &fdtd::options::FdtdOptions::TimeParams::time_max,
                     "Maximum simulation time (s)")
      .def_readwrite("method", &fdtd::options::FdtdOptions::TimeParams::method,
                     "Time stepping method");

  py::class_<fdtd::options::FdtdOptions::BoundaryParams>(m, "BoundaryParams")
      .def(py::init<>())
      .def_readwrite("use_pml",
                     &fdtd::options::FdtdOptions::BoundaryParams::use_pml,
                     "Enable PML boundary conditions")
      .def_readwrite("use_sponge",
                     &fdtd::options::FdtdOptions::BoundaryParams::use_sponge,
                     "Enable sponge boundary conditions")
      .def_readwrite("pml_size",
                     &fdtd::options::FdtdOptions::BoundaryParams::pml_size,
                     "Number of PML layers")
      .def_readwrite("sponge_size",
                     &fdtd::options::FdtdOptions::BoundaryParams::sponge_size,
                     "Number of sponge layers")
      .def_readwrite("sponge_alpha",
                     &fdtd::options::FdtdOptions::BoundaryParams::sponge_alpha,
                     "Sponge damping coefficient");

  py::class_<fdtd::options::FdtdOptions::OutputParams>(m, "OutputParams")
      .def(py::init<>())
      .def_readwrite("save_snapshots",
                     &fdtd::options::FdtdOptions::OutputParams::save_snapshots,
                     "Enable wavefield snapshot saving")
      .def_readwrite("snapshot_interval",
                     &fdtd::options::FdtdOptions::OutputParams::snapshot_interval,
                     "Time steps between snapshots");

  // Bind FdtdOptions (main configuration struct)
  py::class_<fdtd::options::FdtdOptions>(m, "FdtdOptions")
      .def(py::init<>())
      .def_readwrite("grid", &fdtd::options::FdtdOptions::grid,
                     "Grid configuration")
      .def_readwrite("stencil", &fdtd::options::FdtdOptions::stencil,
                     "Stencil configuration")
      .def_readwrite("source", &fdtd::options::FdtdOptions::source,
                     "Source configuration")
      .def_readwrite("velocity", &fdtd::options::FdtdOptions::velocity,
                     "Velocity model configuration")
      .def_readwrite("time", &fdtd::options::FdtdOptions::time,
                     "Time stepping configuration")
      .def_readwrite("boundary", &fdtd::options::FdtdOptions::boundary,
                     "Boundary condition configuration")
      .def_readwrite("output", &fdtd::options::FdtdOptions::output,
                     "Output configuration")
      .def("Validate", &fdtd::options::FdtdOptions::Validate,
           "Validate all configuration parameters");

  // ============================================================================
  // FDTD Components
  // ============================================================================

  // Bind FdtdSourceReceivers
  py::class_<FdtdSourceReceivers>(m, "FdtdSourceReceivers")
      .def(py::init<>())
      .def_readwrite("xsrc", &FdtdSourceReceivers::xsrc,
                     "Source x coordinate (grid index)")
      .def_readwrite("ysrc", &FdtdSourceReceivers::ysrc,
                     "Source y coordinate (grid index)")
      .def_readwrite("zsrc", &FdtdSourceReceivers::zsrc,
                     "Source z coordinate (grid index)");

  // Bind FdtdStencils
  py::class_<fdtd::stencils::FdtdStencils>(m, "FdtdStencils")
      .def(py::init<>())
      .def_readwrite("lx", &fdtd::stencils::FdtdStencils::lx,
                     "Half stencil length in x")
      .def_readwrite("ly", &fdtd::stencils::FdtdStencils::ly,
                     "Half stencil length in y")
      .def_readwrite("lz", &fdtd::stencils::FdtdStencils::lz,
                     "Half stencil length in z")
      .def_readwrite("coef0", &fdtd::stencils::FdtdStencils::coef0,
                     "Central stencil coefficient")
      .def_readonly("coefx", &fdtd::stencils::FdtdStencils::coefx,
                    "Stencil coefficients in x direction")
      .def_readonly("coefy", &fdtd::stencils::FdtdStencils::coefy,
                    "Stencil coefficients in y direction")
      .def_readonly("coefz", &fdtd::stencils::FdtdStencils::coefz,
                    "Stencil coefficients in z direction")
      .def("initStencilsCoefficients",
           &fdtd::stencils::FdtdStencils::initStencilsCoefficients,
           py::arg("options"), py::arg("dx"), py::arg("dy"), py::arg("dz"),
           "Initialize stencil coefficients");

  // Bind FdtdKernels
  py::class_<fdtd::kernel::FdtdKernels>(m, "FdtdKernels")
      .def(py::init<>())
      .def("initFieldsArrays",
           &fdtd::kernel::FdtdKernels::initFieldsArrays, py::arg("nx"),
           py::arg("ny"), py::arg("nz"), py::arg("lx"), py::arg("ly"),
           py::arg("lz"), "Initialize wavefield arrays")
      .def_readonly("pnGlobal", &fdtd::kernel::FdtdKernels::pnGlobal,
                    "Pressure field array (2D: [n_points, 2])")
      .def_readonly("RHSTerm", &fdtd::kernel::FdtdKernels::RHSTerm,
                    "Source RHS term (1D: [n_time_steps])")
      .def_readonly("phi", &fdtd::kernel::FdtdKernels::phi,
                    "PML auxiliary field (1D: [n_points])");

  // Bind FdtdAbcKernels
  py::class_<fdtd::abckernel::FdtdAbcKernels>(m, "FdtdAbcKernels")
      .def(py::init<>())
      .def("defineSpongeBoundary",
           &fdtd::abckernel::FdtdAbcKernels::defineSpongeBoundary,
           py::arg("nx"), py::arg("ny"), py::arg("nz"),
           "Define sponge boundary conditions")
      .def_readonly("eta", &fdtd::abckernel::FdtdAbcKernels::eta,
                    "PML damping coefficient array")
      .def_readonly("spongeArray", &fdtd::abckernel::FdtdAbcKernels::spongeArray,
                    "Sponge damping array");

  // Bind FdtdGrids
  py::class_<model::fdgrid::FdtdGrids>(m, "FdtdGrids")
      .def(py::init<>())
      .def("InitGrid", &model::fdgrid::FdtdGrids::InitGrid, py::arg("options"),
           "Initialize grid geometry")
      .def("InitModelArrays", &model::fdgrid::FdtdGrids::InitModelArrays,
           py::arg("options"), "Initialize model arrays")
      .def("nx", &model::fdgrid::FdtdGrids::nx, "Get grid size in x")
      .def("ny", &model::fdgrid::FdtdGrids::ny, "Get grid size in y")
      .def("nz", &model::fdgrid::FdtdGrids::nz, "Get grid size in z")
      .def("dx", &model::fdgrid::FdtdGrids::dx, "Get grid spacing in x")
      .def("dy", &model::fdgrid::FdtdGrids::dy, "Get grid spacing in y")
      .def("dz", &model::fdgrid::FdtdGrids::dz, "Get grid spacing in z");

  // ============================================================================
  // FdtdSolver - Main Solver Class
  // ============================================================================

  py::class_<FdtdSolver>(m, "FdtdSolver")
      .def(py::init<model::fdgrid::FdtdGrids&, fdtd::kernel::FdtdKernels&,
                    fdtd::abckernel::FdtdAbcKernels&,
                    fdtd::stencils::FdtdStencils&, FdtdSourceReceivers&>(),
           py::arg("grids"), py::arg("kernels"), py::arg("abckernels"),
           py::arg("stencils"), py::arg("source_receivers"),
           "Construct FDTD solver with all required components")
      .def("compute_one_stepSB", &FdtdSolver::compute_one_stepSB,
           py::arg("itime"), py::arg("i1"), py::arg("i2"),
           "Compute one time step with sponge boundary conditions")
      .def("compute_one_stepPML", &FdtdSolver::compute_one_stepPML,
           py::arg("itime"), py::arg("i1"), py::arg("i2"),
           "Compute one time step with PML boundary conditions");

  // ============================================================================
  // Helper Functions
  // ============================================================================

  m.def(
      "create_fd_solver",
      [](model::fdgrid::FdtdGrids& grids, fdtd::kernel::FdtdKernels& kernels,
         fdtd::abckernel::FdtdAbcKernels& abckernels,
         fdtd::stencils::FdtdStencils& stencils,
         FdtdSourceReceivers& source_receivers) {
        return FdtdSolver(grids, kernels, abckernels, stencils,
                          source_receivers);
      },
      py::arg("grids"), py::arg("kernels"), py::arg("abckernels"),
      py::arg("stencils"), py::arg("source_receivers"),
      "Factory function to create an FdtdSolver instance");

  m.def(
      "initialize_rhs_term",
      [](fdtd::kernel::FdtdKernels& kernels, int n_time_steps, float dt,
         float f0, int source_order) {
        // Allocate RHS term
        kernels.RHSTerm = allocateVector<vectorReal>(n_time_steps, "RHSTerm");

        // Compute Ricker wavelet
        float t0 = 1.2f / f0;  // Time delay
        for (int i = 0; i < n_time_steps; i++) {
          float t = i * dt;
          float tau = t - t0;
          float omega = 2.0f * M_PI * f0;
          float omega2 = omega * omega;
          float tau2 = tau * tau;

          if (source_order == 0) {
            // Ricker wavelet (2nd derivative of Gaussian)
            kernels.RHSTerm(i) =
                (1.0f - 2.0f * omega2 * tau2) * std::exp(-omega2 * tau2);
          } else if (source_order == 1) {
            // 1st derivative
            kernels.RHSTerm(i) =
                -4.0f * omega2 * tau * (1.5f - omega2 * tau2) *
                std::exp(-omega2 * tau2);
          } else if (source_order == 2) {
            // 2nd derivative
            kernels.RHSTerm(i) =
                4.0f * omega2 * (omega2 * tau2 * (omega2 * tau2 - 3.0f) + 1.5f) *
                std::exp(-omega2 * tau2);
          } else {
            throw std::runtime_error("Unsupported source order");
          }
        }
      },
      py::arg("kernels"), py::arg("n_time_steps"), py::arg("dt"),
      py::arg("f0"), py::arg("source_order") = 2,
      "Initialize RHS term (source wavelet) using Ricker wavelet\n\n"
      "Parameters:\n"
      "  kernels: FdtdKernels object\n"
      "  n_time_steps: Number of time steps\n"
      "  dt: Time step size (s)\n"
      "  f0: Source frequency (Hz)\n"
      "  source_order: Derivative order (0, 1, or 2)\n");
}
