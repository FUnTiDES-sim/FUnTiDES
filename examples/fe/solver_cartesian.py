"""
This module runs the solver using a cartesian model with:
  - Kokkos GPU or CPU memory space
  - Structured or unstructured cartesian mesh
  - Polynomial order 1, 2 or 3
  - Implementation type: CLASSIC, MAKUTU, OPTIM or SHIVA
  - Only acoustic physics is implemented for now
It demonstrates usage of pybind11 wrapped C++ classes and functions from the proxys library.
For help run with the --help option.
"""

import argparse
import os
import time
from datetime import datetime
from enum import Enum
import sys

import matplotlib.pyplot as plt
import numpy as np
import kokkos  # must be imported after matplotlib to avoid conflicts
import pykokkos as pk
import pyfuntides.model as Model
import pyfuntides.solver as Solver

# Create alias to use float32 and int32 types
CartesianStructBuilderFI1 = Model.CartesianStructBuilder_f32_i32_O1
CartesianStructBuilderFI2 = Model.CartesianStructBuilder_f32_i32_O2
CartesianStructBuilderFI3 = Model.CartesianStructBuilder_f32_i32_O3
CartesianUnstructBuilder = Model.CartesianUnstructBuilder_f32_i32
CartesianParams = Model.CartesianParams_f32_i32

# Avoid taking the enitre dev node for this example
os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OMP_THREAD_LIMIT", "6")
os.environ.setdefault("KOKKOS_NUM_THREADS", "6")


class MemSpace(Enum):
    """
    Memory space options for Kokkos.

    Attributes
    ----------
    CPU : str
        Host memory space ("HostSpace").
    GPU : str
        CUDA memory space ("CudaUVMSpace").
    """

    CPU = "HostSpace"
    GPU = "CudaUVMSpace"


class ModelType(Enum):
    """
    Cartesian model type options.

    Attributes
    ----------
    STRUCTURED : str
        On-the-fly generated mesh, no stored arrays.
    UNSTRUCTURED : str
        Stores mesh in arrays (coordinates, etc.)
    """

    STRUCTURED = "Structured"
    UNSTRUCTURED = "Unstructured"


class ImplemType(Enum):
    """
    Implementation type options.

    Attributes
    ----------
    CLASSIC : str
        Classic implementation.
    GEOS : str
        Geos implementation.
    OPTIM : str
        Optim implementation.
    SHIVA : str
        Shiva implementation.
    """

    CLASSIC = "Classic"
    MAKUTU = "MAKUTU"
    OPTIM = "Optim"
    SHIVA = "Shiva"


def detect_default_memspace():
    """
    Attempts to detect the default memory space expected by the C++ solver bindings
    by inspecting the docstring of the data constructor.
    Returns MemSpace.GPU.name if CUDA types are found, otherwise MemSpace.CPU.name.
    """
    try:
        # Check the docstring of the __init__ method for C++ signature
        doc = Solver.WavefieldAcoustic.__init__.__doc__
        if doc and ("CudaUVMSpace" in doc or "CudaUVMSpace" in doc):
            return MemSpace.GPU.name
        if doc and "HostSpace" in doc:
            return MemSpace.CPU.name
    except Exception:
        pass
    return MemSpace.CPU.name


def parse_args():
    """
    Parses command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """

    # Auto-detect default memory space based on compiled C++ bindings
    default_mem = detect_default_memspace()

    parser = argparse.ArgumentParser(
        description="Run FE Cartesian solver with Kokkos memspace selection."
    )
    parser.add_argument(
        "--mem",
        choices=[e.name for e in MemSpace],
        default=default_mem,
        help=f"Choose Kokkos memspace: {', '.join(e.name for e in MemSpace)} (default: {default_mem} [auto-detected])",
    )
    parser.add_argument(
        "--model",
        choices=[e.name for e in ModelType],
        default=ModelType.STRUCTURED.name,
        help=f"Choose model type: {', '.join(e.name for e in ModelType)} (default: {ModelType.STRUCTURED.name})",
    )
    parser.add_argument(
        "--impl",
        choices=[e.name for e in ImplemType],
        default=ImplemType.MAKUTU.name,
        help=f"Choose implementation type: {', '.join(e.name for e in ImplemType)} (default: {ImplemType.MAKUTU.name})",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=2,
        choices=range(1, 4),
        help="Polynomial order of the elements (default: 2, max 3)",
    )
    parser.add_argument(
        "--domain_size",
        type=float,
        default=1500.0,
        help="Size of the cubic domain (default: 1500.0)",
    )
    parser.add_argument(
        "--ex",
        type=int,
        default=100,
        help="Number of elements in x-direction (default: 50)",
    )
    parser.add_argument(
        "--ey",
        type=int,
        default=100,
        help="Number of elements in y-direction (default: 50)",
    )
    parser.add_argument(
        "--ez",
        type=int,
        default=100,
        help="Number of elements in z-direction (default: 50)",
    )
    parser.add_argument(
        "--f0",
        type=float,
        default=5.0,
        help="Peak frequency for the Ricker source term (default: 5.0)",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=0.001,
        help="Time step size (default: 0.001)",
    )
    parser.add_argument(
        "--n_time_steps",
        type=int,
        default=1500,
        help="Number of time steps to run (default: 1500)",
    )
    parser.add_argument(
        "--n_rhs",
        type=int,
        default=2,
        help="Number of right-hand side sources (default: 2)",
    )
    parser.add_argument(
        "--on_nodes",
        action="store_true",
        default=False,
        help="Whether to apply model on nodes (default: False)",
    )
    # Sponge boundary arguments
    parser.add_argument(
        "--boundaries_size",
        type=float,
        default=0.0,
        help="Size of absorbing boundaries in meters (default: 0)",
    )
    parser.add_argument(
        "--surface_sponge",
        action="store_true",
        default=False,
        help="Enable sponge at the free surface (default: False)",
    )
    parser.add_argument(
        "--taper_delta",
        type=float,
        default=0.015,
        help="Taper delta for sponge boundaries (default: 0.015)",
    )
    return parser.parse_args()


def select_kokkos_memspace(memspace_arg):
    """
    Select the Kokkos memory space and layout.

    Parameters
    ----------
    memspace_arg : str
        The memory space argument, either 'CPU' or 'GPU'.

    Returns
    -------
    memspace : kokkos.Space
        The selected Kokkos memory space.
    layout : kokkos.Layout
        The selected Kokkos layout.
    """

    try:
        enum_value = MemSpace[memspace_arg]
    except KeyError:
        raise ValueError(f"Unknown python memory space: {memspace_arg}")
    if enum_value == MemSpace.CPU:
        memspace = kokkos.HostSpace
        layout = kokkos.LayoutRight
    else:
        # Use standard VRAM space
        memspace = kokkos.CudaUVMSpace
        layout = kokkos.LayoutLeft
    return memspace, layout


def get_solver_model_type(model_type):
    """
    Map a ModelType value (or name) to the corresponding Solver.MeshType.
    """
    try:
        enum_value = ModelType[model_type]
    except KeyError:
        raise ValueError(f"Unknown python model type: {model_type}")
    if enum_value == ModelType.STRUCTURED:
        return Solver.MeshType.STRUCT
    elif enum_value == ModelType.UNSTRUCTURED:
        return Solver.MeshType.UNSTRUCT
    else:
        raise ValueError(f"Unknown solver model type for: {enum_value.name}")


def get_solver_implem_type(implem_type):
    """
    Map an implementation identifier (name or enum) to the corresponding Solver.ImplemType.
    """
    try:
        enum_value = ImplemType[implem_type]
    except KeyError:
        raise ValueError(f"Unknown python implementation type: {implem_type}")
    if enum_value == ImplemType.CLASSIC:
        return Solver.ImplemType.CLASSIC
    elif enum_value == ImplemType.MAKUTU:
        return Solver.ImplemType.MAKUTU
    elif enum_value == ImplemType.OPTIM:
        return Solver.ImplemType.OPTIM
    elif enum_value == ImplemType.SHIVA:
        return Solver.ImplemType.SHIVA
    else:
        raise ValueError(
            f"Unknown solver implementation type for: {enum_value.name}"
        )


def create_model(model_type, e, h, l, order, on_nodes):
    """
    Create a Cartesian model based on the specified type.
    """
    try:
        enum_value = ModelType[model_type]
    except KeyError:
        raise ValueError(f"Unknown python model type: {model_type}")
    if enum_value == ModelType.STRUCTURED:
        return create_structured_model(e, l, order, on_nodes)
    elif enum_value == ModelType.UNSTRUCTURED:
        return create_unstructured_model(e, l, order, on_nodes)
    else:
        raise ValueError(f"Unknown model type: {enum_value.name}")


def create_structured_model(e, l, order, on_nodes):
    """
    Create a structured Cartesian model based on the specified order.
    """
    if order == 1:
        builder = CartesianStructBuilderFI1(
            e[0], l[0], e[1], l[1], e[2], l[2], on_nodes, False
        )
    elif order == 2:
        builder = CartesianStructBuilderFI2(
            e[0], l[0], e[1], l[1], e[2], l[2], on_nodes, False
        )
    elif order == 3:
        builder = CartesianStructBuilderFI3(
            e[0], l[0], e[1], l[1], e[2], l[2], on_nodes, False
        )
    else:
        raise ValueError(
            f"Order {order} is not wrapped by pybind11 (only 1, 2, 3 supported)"
        )
    return builder.get_model()


def create_unstructured_model(e, l, order, on_nodes):
    """
    Create an unstructured Cartesian model.
    """
    if order not in (1, 2, 3):
        raise ValueError(
            f"Order {order} is not wrapped by pybind11 (only 1, 2, 3 supported)"
        )
    params = CartesianParams()
    params.ex, params.ey, params.ez = e
    params.lx, params.ly, params.lz = l
    params.order = order
    params.is_model_on_nodes = on_nodes
    params.is_elastic = False
    builder = CartesianUnstructBuilder(params)
    return builder.get_model()


def create_solver(implem_type, model_type, order, on_nodes):
    """
    Create a solver based on the specified implementation type.
    """
    impl = get_solver_implem_type(implem_type)
    model = get_solver_model_type(model_type)
    model_location = (
        Solver.ModelLocationType.ONNODES
        if on_nodes
        else Solver.ModelLocationType.ONELEMENTS
    )
    physic_type = Solver.PhysicType.ACOUSTIC

    return Solver.create_solver(
        Solver.MethodType.SEM, impl, model, model_location, Solver.PhysicType.ACOUSTIC, order
    )


def source_term(time_n, f0):
    """
    Computes the source term value at a given time for a Ricker wavelet.
    """
    o_tpeak = 1.0 / f0
    pulse = 0.0
    if time_n <= -0.9 * o_tpeak or time_n >= 2.9 * o_tpeak:
        return pulse

    pi = 3.14157
    lam = (f0 * pi) * (f0 * pi)
    pulse = (
        2.0
        * lam
        * (2.0 * lam * (time_n - o_tpeak) * (time_n - o_tpeak) - 1.0)
        * np.exp(-lam * (time_n - o_tpeak) * (time_n - o_tpeak))
    )
    return pulse


def get_snapshot(nx, ny, nz, pnGlobal, normalize=False):
    """
    Extracts a 2D snapshot from a 3D global array at a specified index, with optional normalization.
    """
    offset = nx * nz * (int(ny / 2) - 1)
    grid = np.zeros((nx, nz))
    for I in range(offset, offset + nx * nz):
        i = (I - offset) % nx
        j = int((I - offset - i) / nx)
        grid[i, j] = pnGlobal[I]

    if normalize:
        maxvalue = np.abs(grid).max()
        if maxvalue != 0:
            grid = grid / maxvalue

    return grid


def setup_plot(nx, nz, cmpvalue=0.15):
    """
    Set up a matplotlib plot for a 2D slice of a Float32 array.
    """
    grid = np.zeros((nx, nz))
    fig, ax = plt.subplots()
    im = ax.imshow(grid, cmap="viridis", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Intensity")
    plt.title("2D Slice of a Float32 Array")
    plt.xlabel("X-axis")
    plt.ylabel("Z-axis")
    plt.ioff()  # Prevent showing the plot interactively
    return fig, ax, im


def plot_snapshot(nx, ny, nz, kk_pnGlobal, im, t):
    """
    Plot and save a snapshot of the simulation at the given time step.
    Passes the Kokkos array so it can be safely mirrored to the host.
    """
    # 1. Create a host mirror view
    kk_pnGlobal_host = kokkos.create_mirror_view(kk_pnGlobal)

    # 2. Deep copy from device to host
    kokkos.deep_copy(kk_pnGlobal_host, kk_pnGlobal)

    # 3. Create numpy wrapper over the safe host memory
    pnGlobal_host = np.array(kk_pnGlobal_host, copy=False)

    grid = get_snapshot(nx, ny, nz, pnGlobal_host, False)
    im.set_array(grid)  # Update plot with new values
    plt.draw()  # Redraw the figure with updated data
    plt.ioff()
    plt.savefig(f"snap0{t:0{5}d}.png")


def allocate_pressure(n_dof, memspace, layout):
    """
    Allocate and initialize to 0 the pressure arrays.
    """
    kk_pnGlobalPrev = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    kk_pnGlobalCurr = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pk.fence()
    # Note: These numpy arrays point to device memory if using CudaUVMSpace!
    # Do not loop over them directly in Python.
    pnGlobalPrev = np.array(kk_pnGlobalPrev, copy=False)
    pnGlobalCurr = np.array(kk_pnGlobalCurr, copy=False)
    return kk_pnGlobalPrev, pnGlobalPrev, kk_pnGlobalCurr, pnGlobalCurr


def allocate_rhs_term(n_rhs, n_time_steps, dt, f0, memspace, layout):
    """
    Allocate and fill the RHSTerm array for the source term.
    """
    kk_RHSTerm = kokkos.array(
        [n_rhs, n_time_steps], dtype=kokkos.float32, space=memspace, layout=layout
    )

    # Create Host Mirror for initialization
    kk_RHSTerm_host = kokkos.create_mirror_view(kk_RHSTerm)
    RHSTerm_host = np.array(kk_RHSTerm_host, copy=False)

    for i in range(n_time_steps):
        RHSTerm_host[0, i] = source_term(i * dt, f0)
        RHSTerm_host[1, i] = source_term(i * dt, f0)

    # Copy initialized data to device
    kokkos.deep_copy(kk_RHSTerm, kk_RHSTerm_host)

    return kk_RHSTerm, RHSTerm_host


def allocate_rhs_weight(n_rhs, model, memspace, layout):
    """
    Allocate and fill the RHSWeights array.
    """
    nb_points = model.get_number_of_points_per_element()
    kk_RHSWeights = kokkos.array(
        [n_rhs, nb_points],
        dtype=kokkos.float32,
        space=memspace,
        layout=layout,
    )

    # Create Host Mirror for initialization
    kk_RHSWeights_host = kokkos.create_mirror_view(kk_RHSWeights)
    RHSWeights_host = np.array(kk_RHSWeights_host, copy=False)

    weight_val = 1.0 / nb_points
    for i in range(n_rhs):
        for j in range(nb_points):
            RHSWeights_host[i, j] = weight_val

    # Copy initialized data to device
    kokkos.deep_copy(kk_RHSWeights, kk_RHSWeights_host)

    return kk_RHSWeights, RHSWeights_host


def allocate_rhs_element(n_rhs, ex, ey, ez, memspace, layout):
    """
    Allocate and fill the RHSElement array.
    """
    kk_RHSElement = kokkos.array(
        [n_rhs], dtype=kokkos.int32, space=memspace, layout=layout
    )

    # Create Host Mirror for safe potential manipulation
    kk_RHSElement_host = kokkos.create_mirror_view(kk_RHSElement)
    RHSElement_host = np.array(kk_RHSElement_host, copy=False)

    # If you need to manipulate RHSElement in python, do it here on the host array
    # RHSElement_host[0] = int(ex / 2 + ey / 2 * ex + ez / 2 * ey * ex)
    # RHSElement_host[1] = int(ex / 3 + ey / 2 * ex + ez / 2 * ey * ex)

    # Copy to device
    kokkos.deep_copy(kk_RHSElement, kk_RHSElement_host)

    return kk_RHSElement, RHSElement_host


def create_solver_data(kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RHSElement, kk_RHSWeights):
    """
    Create SEMsolverData instance and associated wavefield and rhs.
    """
    wavefield = Solver.WavefieldAcoustic(kk_pnGlobalPrev, kk_pnGlobalCurr)
    rhs = Solver.RhsAcoustic(kk_RHSTerm, kk_RHSElement, kk_RHSWeights)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)
    data.print()

    return wavefield, rhs, data


def compute_step(
    time_sample,
    dt,
    solver,
    data,
    iteration_times,
    n_time_steps,
    nx,
    ny,
    nz,
    kk_pnGlobal,  # Pass the Kokkos array!
    im,
):
    """
    Perform a single time step of the simulation, update timing, plot, and swap indices.
    """
    iter_start = time.time()

    # 1. Compute forces (RHS of the equation)
    solver.compute_forces(dt, time_sample, data)

    # 2. Synchronize boundaries (Placeholder for distributed logic)
    # m_syncer->synchronize(m_solver->getForceVector(c), par_topology_);

    # 3. Update solution using mass matrix and accumulated forces
    solver.update_solution(dt, data)

    iter_time = time.time() - iter_start
    iteration_times.append(iter_time)

    if time_sample % 1000 == 0:
        print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
        print()
    if time_sample % 100 == 0:
        print(f"Time {time_sample} / {n_time_steps}")
    if time_sample % 10 == 0:
        # Pass the Kokkos array so the plotter can mirror it back to the host
        plot_snapshot(nx, ny, nz, kk_pnGlobal, im, time_sample)


def main():
    # Parse command line arguments
    args = parse_args()

    # Initialize global parameters from command-line arguments
    on_nodes = args.on_nodes
    f0 = args.f0
    dt = args.dt
    n_time_steps = args.n_time_steps
    n_rhs = args.n_rhs
    order = args.order
    domain_size = args.domain_size
    lx = ly = lz = domain_size
    ex = args.ex
    ey = args.ey
    ez = args.ez
    hx = lx / ex
    hy = ly / ey
    hz = lz / ez
    nx = ex * order + 1
    ny = ey * order + 1
    nz = ez * order + 1
    n_dof = nx * ny * nz
    n_elements = ex * ey * ez
    n_points_per_elements = (order + 1) * (order + 1) * (order + 1)

    print("==========SIMULATION PARAMETERS==========")
    print(f"order                        : {order}")
    print(f"on_nodes                     : {on_nodes}")
    print(f"memspace                     : {args.mem}")
    print(f"impl                         : {args.impl}")
    print(f"model                        : {args.model}")
    print(f"number of elements           : {n_elements}")
    print(f"number of points per element : {n_points_per_elements}")
    print(f"f0                           : {f0}")
    print(f"dt                           : {dt}")
    print(f"n_time_steps                 : {n_time_steps}")
    print(f"n_rhs                        : {n_rhs}")
    print(f"boundaries size              : {args.boundaries_size}")
    print("=========================================")

    # Setup graphic display
    print("Setting up plot...")
    _, _, im = setup_plot(nx, nz)
    print("Plot set up")

    # Initialize Kokkos
    kokkos.initialize()
    print("Kokkos initialized")
    memspace, layout = select_kokkos_memspace(args.mem)

    # Add timing variables
    start_time = time.time()
    simulation_start = datetime.now()
    iteration_times = []
    print(f"Simulation started at: {simulation_start}")

    # Create model
    print("Creating model...")
    model = create_model(
        args.model, (ex, ey, ez), (hx, hy, hz), (lx, ly, lz), order, on_nodes
    )
    print("Model created")

    # Create solver
    print("Creating solver...")
    solver = create_solver(args.impl, args.model, order, on_nodes)
    print("Solver created")

    # Initialize model
    print("Initializing model...")
    sponge_size = [args.boundaries_size, args.boundaries_size, args.boundaries_size]
    solver.compute_fe_init(model, sponge_size, args.surface_sponge, args.taper_delta)
    print("Model initialized")

    # allocate pressure
    print("Allocating Pressure...")
    kk_pnGlobalPrev, pnGlobalPrev, kk_pnGlobalCurr, pnGlobalCurr = allocate_pressure(n_dof, memspace, layout)
    print("Pressure allocated")

    # allocate RHS arrays
    print("Allocating RHS element...")
    kk_RHSElement, RHSElement = allocate_rhs_element(
        n_rhs, ex, ey, ez, memspace, layout
    )
    print("RHS element allocated")

    print("Allocating RHS weights...")
    kk_RHSWeights, rhsWeights = allocate_rhs_weight(n_rhs, model, memspace, layout)
    print("RHS weights allocated")

    print("Allocating RHS term...")
    kk_RHSTerm, rhsTerm = allocate_rhs_term(
        n_rhs, n_time_steps, dt, f0, memspace, layout
    )
    print("RHS term allocated")

    # Create solver data instance
    print("Creating solver data...")
    wavefield, rhs, data = create_solver_data(kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RHSElement, kk_RHSWeights)
    print("Solver data created")

    # Loop over time steps
    for time_sample in range(n_time_steps):
        compute_step(
            time_sample,
            dt,
            solver,
            data,
            iteration_times,
            n_time_steps,
            nx,
            ny,
            nz,
            kk_pnGlobalPrev,  # Pass Kokkos array instead of numpy view!
            im,
        )

        # Swap pressure arrays for next iteration
        data.swap_wavefields()

    # Print final timing statistics
    end_time = time.time()
    simulation_end = datetime.now()
    total_time = end_time - start_time

    print("==========SIMULATION STATISTICS==========")
    print(f"{'Start time:':<25} {simulation_start}")
    print(f"{'End time:':<25} {simulation_end}")
    print(f"{'Total runtime:':<25} {total_time:.2f} seconds")
    print(f"{'Average iteration time:':<25} {np.mean(iteration_times):.4f} seconds")
    print(f"{'Min iteration time:':<25} {np.min(iteration_times):.4f} seconds")
    print(f"{'Max iteration time:':<25} {np.max(iteration_times):.4f} seconds")
    print("=========================================")

    # release kokkos arrays and vectors
    del data
    del solver
    del model

    # release explicit views
    del kk_pnGlobalPrev
    del kk_pnGlobalCurr
    del kk_RHSTerm
    del kk_RHSElement
    del kk_RHSWeights

    # release numpy wrappers which might hold references to views
    del pnGlobalPrev
    del pnGlobalCurr
    del rhsTerm
    del RHSElement
    del rhsWeights

    kokkos.finalize()
    print("End of  computation")


if __name__ == "__main__":
    main()
