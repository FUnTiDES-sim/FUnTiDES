#!/usr/bin/env python3
"""
This module runs the solver using a cartesian model with:
  - Kokkos GPU or CPU memory space
  - Structured or unstructured cartesian mesh
  - Polynomial order 1, 2 or 3
  - Implementation type: CLASSIC, GEOS, OPTIM or SHIVA
It demonstrates usage of pybind11 wrapped C++ classes and functions from the proxys library.
For help run with the --help option.
"""

import argparse
import time
from datetime import datetime
from enum import Enum

import kokkos
import matplotlib.pyplot as plt
import numpy as np
import pyproxys.model as Model
import pyproxys.solver as Solver


class MemSpace(Enum):
    """
    Memory space options for Kokkos.

    Attributes
    ----------
    CPU : str
        Host memory space ("HostSpace").
    GPU : str
        CUDA Unified Virtual Memory space ("CudaUVMSpace").
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
    GEOS = "Geos"
    OPTIM = "Optim"
    SHIVA = "Shiva"


def parse_args():
    """
    Parses command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Run FE Cartesian solver with Kokkos memspace selection."
    )
    parser.add_argument(
        "--mem",
        choices=[e.value for e in MemSpace],
        default=MemSpace.CPU.value,
        help="Choose Kokkos memspace: 'HostSpace' (CPU, default) or 'CudaUVMSpace' (GPU)",
    )
    parser.add_argument(
        "--model",
        choices=[e.value for e in ModelType],
        default=ModelType.STRUCTURED.value,
        help="Choose model type: 'Structured' (default) or 'Unstructured'",
    )
    parser.add_argument(
        "--impl",
        choices=[e.name for e in ImplemType],
        default=ImplemType.SHIVA.value,
        help=f"Choose implementation type: {', '.join(e.name for e in ImplemType)} (default: SHIVA)",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=2,
        choices=range(1, 4),
        help="Polynomial order of the elements (default: 3, max 3)",
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
    in_memspace = MemSpace(memspace_arg)
    if in_memspace == MemSpace.CPU:
        memspace = kokkos.HostSpace
        layout = kokkos.LayoutRight
    else:
        memspace = kokkos.CudaUVMSpace
        layout = kokkos.LayoutLeft
    return memspace, layout


def get_solver_model_type(model_type):
    """
    Map a ModelType value (or name) to the corresponding Solver.MeshType.

    Parameters
    ----------
    model_type : str or ModelType
        Model type name or ModelType enum value. Accepted values are
        'Structured' / ModelType.STRUCTURED and 'Unstructured' / ModelType.UNSTRUCTED.

    Returns
    -------
    Solver.MeshType
        Corresponding Solver.MeshType enum (Solver.MeshType.STRUCT or Solver.MeshType.UNSTRUCT).

    Raises
    ------
    ValueError
        If the provided model_type is unknown or unsupported.
    """
    match ModelType(model_type):
        case ModelType.STRUCTURED:
            return Solver.MeshType.STRUCT
        case ModelType.UNSTRUCTURED:
            return Solver.MeshType.UNSTRUCT
        case _:
            raise ValueError(f"Unknown model type: {model_type}")


def get_solver_implem_type(implem_type):
    """
    Map an implementation identifier (name or enum) to the corresponding Solver.ImplemType.

    Parameters
    ----------
    implem_type : str or ImplemType
        Implementation name or ImplemType enum. Accepted names are
        'CLASSIC', 'GEOS', 'OPTIM', 'SHIVA' (case-insensitive when passed as enum names).

    Returns
    -------
    Solver.ImplemType
        Corresponding Solver.ImplemType enum value.

    Raises
    ------
    ValueError
        If the provided implem_type is unknown or unsupported.
    """
    match ImplemType(implem_type):
        case ImplemType.CLASSIC:
            return Solver.ImplemType.CLASSIC
        case ImplemType.GEOS:
            return Solver.ImplemType.GEOS
        case ImplemType.OPTIM:
            return Solver.ImplemType.OPTIM
        case ImplemType.SHIVA:
            return Solver.ImplemType.SHIVA
        case _:
            raise ValueError(f"Unknown implementation type: {implem_type}")


def create_model(model_type, e, h, order):
    """
    Create a Cartesian model based on the specified type.

    Parameters
    ----------
    model_type : str
        The type of model to create, either 'Structured' or 'Unstructured'.
    e : int or tuple of int
        Number of elements in each dimension (ex, ey, ez).
    h : int or tuple of float
        Element sizes in each dimension (hx, hy, hz).
    order : int
        The polynomial order of the elements.

    Returns
    -------
    model : Model.ModelStruct or Model.ModelUnstruct
        The created Cartesian model.

    Raises
    ------
    ValueError
        If the model type is unknown.
    """
    match ModelType(model_type):
        case ModelType.STRUCTURED:
            return create_structured_model(e[0], h[0], order)
        case ModelType.UNSTRUCTURED:
            return create_unstructured_model(e, h, order)
        case _:
            raise ValueError(f"Unknown model type: {model_type}")


def create_structured_model(e, h, order):
    """
    Create a structured Cartesian model based on the specified order.

    Parameters
    ----------
    e : int
        Number of elements in each dimension (ex, ey, ez).
    h : int
        Element sizes in each dimension (hx, hy, hz).
    order : int
        The polynomial order of the elements.

    Returns
    -------
    model : Model.ModelStruct
        The created structured Cartesian model.

    Raises
    ------
    ValueError
        If the order is not 1, 2, or 3.
    """
    match order:
        case 1:
            builder = Model.CartesianStructBuilderFI1()
        case 2:
            builder = Model.CartesianStructBuilderFI2()
        case 3:
            builder = Model.CartesianStructBuilderFI3()
        case _:
            raise ValueError(
                f"Order {order} is not wrapped by pybind11 (only 1, 2, 3 supported)"
            )
    return builder.get_model(e, h)


def create_unstructured_model(e, h, order):
    """
    Create an unstructured Cartesian model.

    Parameters
    ----------
    e : tuple of int
        Number of elements in each dimension (ex, ey, ez).
    h : tuple of float
        Element sizes in each dimension (hx, hy, hz).
    order : int
        The polynomial order of the elements.

    Returns
    -------
    model : Model.ModelUnstruct
        The created unstructured Cartesian model.

    Raises
    ------
    ValueError
        If the order is not 1, 2, or 3.
    """
    if order not in (1, 2, 3):
        raise ValueError(
            f"Order {order} is not wrapped by pybind11 (only 1, 2, 3 supported)"
        )
    params = Model.CartesianParams()
    params.ex, params.ey, params.ez = e
    params.hx, params.hy, params.hz = h
    params.order = order
    builder = Model.CartesianUnstructBuilder(params)
    return builder.getModel()


def create_solver(implem_type, model_type, order):
    """
    Create a solver based on the specified implementation type.

    Parameters
    ----------
    implem_type : str
        The implementation type, one of 'CLASSIC', 'GEOS', 'OPTIM', or 'SHIVA'.
    model_type : str
        The model type, either 'Structured' or 'Unstructured'.
    order : int
        The polynomial order of the elements.

    Returns
    -------
    solver : Solver.Solver
        The created solver.

    Raises
    ------
    ValueError
        If the implementation type is unknown.
    """
    impl = get_solver_implem_type(implem_type)
    model = get_solver_model_type(model_type)

    return Solver.create_solver(Solver.MethodType.SEM, impl, model, order)


def source_term(time_n, f0):
    """
    Computes the source term value at a given time for a Ricker wavelet.

    Parameters
    ----------
    time_n : float
        The current time at which to evaluate the source term.
    f0 : float
        The peak frequency of the Ricker wavelet.

    Returns
    -------
    float
        The value of the source term at the specified time.

    Notes
    -----
    The function returns zero outside the interval [-0.9 * t_peak, 2.9 * t_peak], where t_peak = 1.0 / f0.
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


def get_snapshot(i1, nx, ny, nz, pnGlobal, normalize=False):
    """
    Extracts a 2D snapshot from a 3D global array at a specified index, with optional normalization.

    Parameters
    ----------
    i1 : int
        Index for the pressure field.
    nx : int
        Number of grid points in the x-direction.
    ny : int
        Number of grid points in the y-direction.
    nz : int
        Number of grid points in the z-direction.
    pnGlobal : np.ndarray
        The global 3D array of shape (nx * ny * nz, N), where N >= i1 + 1.
    normalize : bool
        If True, normalize the resulting 2D grid by its maximum absolute value (default: False).

    Returns
    -------
    grid : np.ndarray
        A 2D array of shape (nx, nz) representing the extracted snapshot.
    """
    offset = nx * nz * (int(ny / 2) - 1)
    grid = np.zeros((nx, nz))
    for I in range(offset, offset + nx * nz):
        i = (I - offset) % nx
        j = int((I - offset - i) / nx)
        grid[i, j] = pnGlobal[I, i1]

    if normalize:
        maxvalue = np.abs(grid).max()
        if maxvalue != 0:
            grid = grid / maxvalue

    return grid


def setup_plot(nx, nz, cmpvalue=0.15):
    """
    Set up a matplotlib plot for a 2D slice of a Float32 array.

    Parameters
    ----------
    nx : int
        Number of grid points in the x-direction.
    nz : int
        Number of grid points in the z-direction.
    cmpvalue : float, optional
        Color map value range (default: 0.15).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    ax : matplotlib.axes.Axes
        The created axes.
    im : matplotlib.image.AxesImage
        The image object for updating the plot.
    """
    grid = np.zeros((nx, nz))
    fig, ax = plt.subplots()
    im = ax.imshow(
        grid, cmap="viridis", interpolation="nearest"
    )
    plt.colorbar(im, ax=ax, label="Intensity")
    plt.title("2D Slice of a Float32 Array")
    plt.xlabel("X-axis")
    plt.ylabel("Z-axis")
    plt.ioff()  # Prevent showing the plot interactively
    return fig, ax, im


def plot_snapshot(i1, nx, ny, nz, pnGlobal, im, t):
    """
    Plot and save a snapshot of the simulation at the given time step.

    Parameters
    ----------
    i1 : int
        Index for the pressure field.
    nx, ny, nz : int
        Grid dimensions.
    pnGlobal : np.ndarray
        Pressure field array.
    im : matplotlib.image.AxesImage
        The image object for updating the plot.
    t : int
        Current time step.
    """
    grid = get_snapshot(i1, nx, ny, nz, pnGlobal, False)
    im.set_array(grid)  # Update plot with new values
    plt.draw()  # Redraw the figure with updated data
    plt.ioff()
    plt.savefig(f"snap0{t:0{5}d}.png")


def allocate_pressure(n_dof, memspace, layout):
    """
    Allocate and initialize to 0 the pressure arrays.

    Parameters
    ----------
    n_dof : int
        Number of degrees of freedom.
    memspace : kokkos.Space
        The selected Kokkos memory space.
    layout : kokkos.Layout
        The selected Kokkos layout.

    Returns
    -------
    kk_pnGlobal : kokkos array
        The Kokkos array for pressure.
    pnGlobal : np.ndarray
        The numpy array view of the pressure.
    """
    kk_pnGlobal = kokkos.array(
        [n_dof, 2], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pnGlobal = np.array(kk_pnGlobal, copy=False)
    pnGlobal[:] = 0.0

    return kk_pnGlobal, pnGlobal


def allocate_rhs_term(n_rhs, n_time_steps, dt, f0, memspace, layout):
    """
    Allocate and fill the RHSTerm array for the source term.

    Parameters
    ----------
    n_rhs : int
        Number of right-hand side sources.
    n_time_steps : int
        Number of time steps.
    dt : float
        Time step size.
    f0 : float
        Source frequency.
    memspace : kokkos.Space
        The selected Kokkos memory space.
    layout : kokkos.Layout
        The selected Kokkos layout.

    Returns
    -------
    kk_RHSTerm : kokkos array
        The Kokkos array for the source term.
    RHSTerm : np.ndarray
        The numpy array view of the source term.
    """
    kk_RHSTerm = kokkos.array(
        [n_rhs, n_time_steps], dtype=kokkos.float32, space=memspace, layout=layout
    )
    RHSTerm = np.array(kk_RHSTerm, copy=False)
    for i in range(n_time_steps):
        RHSTerm[0, i] = source_term(i * dt, f0)
        RHSTerm[1, i] = source_term(i * dt, f0)
    return kk_RHSTerm, RHSTerm


def allocate_rhs_weight(n_rhs, model, memspace, layout):
    """
    Allocate and fill the RHSWeights array.

    Parameters
    ----------
    n_rhs : int
        Number of right-hand side sources.
    model : Model.ModelStruct or Model.ModelUnstruct
        The model object with get_number_of_points_per_element methods.
    memspace : kokkos.Space
        The selected Kokkos memory space.
    layout : kokkos.Layout
        The selected Kokkos layout.

    Returns
    -------
    kk_RHSWeights : kokkos array
        The Kokkos array for the weights.
    RHSWeights : np.ndarray
        The numpy array view of the weights.
    """
    nb_points = model.get_number_of_points_per_element()
    kk_RHSWeights = kokkos.array(
        [n_rhs, nb_points],
        dtype=kokkos.float32,
        space=memspace,
        layout=layout,
    )
    RHSWeights = np.array(kk_RHSWeights, copy=False)
    for i in range(n_rhs):
        for j in range(model.get_number_of_points_per_element()):
            RHSWeights[i, j] = 1 / model.get_number_of_points_per_element()
    return kk_RHSWeights, RHSWeights


def allocate_rhs_element(n_rhs, ex, ey, ez, memspace, layout):
    """
    Allocate and fill the RHSElement array.

    Parameters
    ----------
    n_rhs : int
        Number of right-hand side sources.
    n_points_per_elements : int
        Number of points per element.
    memspace : kokkos.Space
        The selected Kokkos memory space.
    layout : kokkos.Layout
        The selected Kokkos layout.

    Returns
    -------
    kk_RHSElement : kokkos array
        The Kokkos array for the element indices.
    RHSElement : np.ndarray
        The numpy array view of the element indices.
    """
    kk_RHSElement = kokkos.array(
        [n_rhs], dtype=kokkos.int32, space=memspace, layout=layout
    )
    RHSElement = np.array(kk_RHSElement, copy=False)
    RHSElement[0] = ex/2 + ey/2 * ex + ez/2 * ey * ex  # one half of slice
    RHSElement[1] = ex/3 + ey/2 * ex + ez/2 * ey * ex  # one third of slice
    return kk_RHSElement, RHSElement


def create_solver_data(kk_RHSTerm, kk_pnGlobal, kk_RHSElement, kk_RHSWeights):
    """
    Create SEMsolverData instance and return it along with i1 and i2.

    Parameters
    ----------
    kk_RHSTerm : kokkos array
        The Kokkos array for the source term.
    kk_pnGlobal : kokkos array
        The Kokkos array for pressure.
    kk_RHSElement : kokkos array
        The Kokkos array for the element indices.
    kk_RHSWeights : kokkos array
        The Kokkos array for the weights.

    Returns
    -------
    data : Solver.SEMsolverData
        The SEMsolverData instance.
    """
    data = Solver.SEMsolverData(
        0,
        1,
        kk_RHSTerm,
        kk_pnGlobal,
        kk_RHSElement,
        kk_RHSWeights,
    )
    data.print()

    return data


def compute_step(
    time_sample,
    dt,
    solver,
    data,
    iteration_times,
    n_time_steps,
    i1, i2,
    nx, ny, nz,
    pnGlobal,
    im
):
    """
    Perform a single time step of the simulation, update timing, plot, and swap indices.

    Parameters
    ----------
    time_sample : int
        Current time step.
    dt : float
        Time step size.
    solver : Solver.Solver
        The solver instance.
    data : Solver.SEMsolverData
        The SEMsolverData instance.
    iteration_times : list
        List to append iteration time.
    n_time_steps : int
        Total number of time steps.
    i1, i2 : int
        Indices for pressure fields.
    nx, ny, nz : int
        Grid dimensions for the plot.
    pnGlobal : np.ndarray
        Pressure field array.
    im : matplotlib.image.AxesImage
        The image object for updating the plot.

    Returns
    -------
    i1, i2 : int
        Updated indices for pressure fields.
    """
    iter_start = time.time()
    print(f"Computing time step {time_sample + 1} / {n_time_steps}")
    solver.compute_one_step(dt, time_sample, data)
    print(f"Time step {time_sample + 1} computed")
    iter_time = time.time() - iter_start
    iteration_times.append(iter_time)
    if time_sample % 1000 == 0:
        print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
        print()
    if time_sample % 100 == 0:
        print(f"Time {time_sample} / {n_time_steps}")
    if time_sample % 10 == 0:
        plot_snapshot(i1, nx, ny, nz, pnGlobal, im, time_sample)
    # Swap pn and pn+1
    tmp = i1
    i1 = i2
    i2 = tmp
    data.i1 = i1
    data.i2 = i2
    return i1, i2


def main():
    # Parse command line arguments
    args = parse_args()

    # Initialize global parameters from command-line arguments
    f0 = args.f0
    dt = args.dt
    n_time_steps = args.n_time_steps
    n_rhs = args.n_rhs
    order = args.order
    domain_size = args.domain_size
    ex = args.ex
    ey = args.ey
    ez = args.ez
    hx = domain_size / ex
    hy = domain_size / ey
    hz = domain_size / ez
    nx = ex * order + 1
    ny = ey * order + 1
    nz = ez * order + 1
    n_dof = nx * ny * nz
    n_elements = ex * ey * ez
    n_points_per_elements = (order + 1) * (order + 1) * (order + 1)

    print("==========SIMULATION PARAMETERS==========")
    print(f"order                        : {order}")
    print(f"memspace                     : {args.mem}")
    print(f"impl                         : {args.impl}")
    print(f"model                        : {args.model}")
    print(f"number of elements           : {n_elements}")
    print(f"number of points per element : {n_points_per_elements}")
    print(f"f0                           : {f0}")
    print(f"dt                           : {dt}")
    print(f"n_time_steps                 : {n_time_steps}")
    print(f"n_rhs                        : {n_rhs}")
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
    model = create_model(args.model, (ex, ey, ez), (hx, hy, hz), order)
    print("Model created")

    # Create solver
    print("Creating solver...")
    solver = create_solver(args.impl, args.model, order)
    print("Solver created")

    # Initialize model
    print("Initializing model...")
    solver.compute_fe_init(model)
    print("Model initialized")

    # allocate pressure
    print("Allocating Pressure...")
    kk_pnGlobal, pnGlobal = allocate_pressure(n_dof, memspace, layout)
    print("Pressure allocated")

    # allocate RHS arrays
    print("Allocating RHS element...")
    kk_RHSElement, RHSElement = allocate_rhs_element(n_rhs, ex, ey, ez, memspace, layout)
    print("  - RHS element number 0", RHSElement[0])
    print("  - RHS element number 1", RHSElement[1])
    print("RHS element allocated")

    print("Allocating RHS weights...")
    kk_RHSWeights, rhsWeights = allocate_rhs_weight(n_rhs, model, memspace, layout)
    print("RHS weights allocated")

    print("Allocating RHS term...")
    kk_RHSTerm, rhsTerm = allocate_rhs_term(n_rhs, n_time_steps, dt, f0, memspace, layout)
    print("RHS term allocated")

    # Create solver data instance
    print("Creating solver data...")
    data = create_solver_data(
        kk_RHSTerm, kk_pnGlobal, kk_RHSElement, kk_RHSWeights
    )
    print("Solver data created")

    i1 = data.i1
    i2 = data.i2

    # Loop over time steps
    for time_sample in range(n_time_steps):

        i1, i2 = compute_step(
            time_sample,
            dt,
            solver,
            data,
            iteration_times,
            n_time_steps,
            i1, i2,
            nx, ny, nz,
            pnGlobal,
            im
        )

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
    del kk_pnGlobal
    del kk_RHSTerm
    del kk_RHSElement
    del kk_RHSWeights
    del solver
    del model

    kokkos.finalize()
    print("End of  computation")


if __name__ == "__main__":
    main()
