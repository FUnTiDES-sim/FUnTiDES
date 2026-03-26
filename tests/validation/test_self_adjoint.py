#!/usr/bin/env python3
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

import pytest
import os
import time
from datetime import datetime
from enum import Enum
import sys

import matplotlib.pyplot as plt
import numpy as np
import kokkos  # must be imported after matplotlib to avoid conflicts
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
        if doc and ("CudaUVMSpace" in doc or "CudaSpace" in doc):
            return MemSpace.GPU.name
        if doc and "HostSpace" in doc:
            return MemSpace.CPU.name
    except Exception:
        pass
    return MemSpace.CPU.name


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

    try:
        enum_value = ModelType[model_type]
    except KeyError:
        raise ValueError(f"Unknown python model type: {model_type}")
    print(enum_value)
    match enum_value:
        case ModelType.STRUCTURED:
            return Solver.MeshType.STRUCT
        case ModelType.UNSTRUCTURED:
            return Solver.MeshType.UNSTRUCT
        case _:
            raise ValueError(f"Unknown solver model type for: {enum_value.name}")


def get_solver_implem_type(implem_type):
    """
    Map an implementation identifier (name or enum) to the corresponding Solver.ImplemType.

    Parameters
    ----------
    implem_type : str or ImplemType
        Implementation name or ImplemType enum. Accepted names are
        'CLASSIC', 'MAKUTU', 'OPTIM', 'SHIVA' (case-insensitive when passed as enum names).

    Returns
    -------
    Solver.ImplemType
        Corresponding Solver.ImplemType enum value.

    Raises
    ------
    ValueError
        If the provided implem_type is unknown or unsupported.
    """

    try:
        enum_value = ImplemType[implem_type]
    except KeyError:
        raise ValueError(f"Unknown python implementation type: {implem_type}")
    match enum_value:
        case ImplemType.CLASSIC:
            return Solver.ImplemType.CLASSIC
        case ImplemType.MAKUTU:
            return Solver.ImplemType.MAKUTU
        case ImplemType.OPTIM:
            return Solver.ImplemType.OPTIM
        case ImplemType.SHIVA:
            return Solver.ImplemType.SHIVA
        case _:
            raise ValueError(
                f"Unknown solver implementation type for: {enum_value.name}"
            )


def create_model(model_type, e, h, l, order, on_nodes, elastic):
    """
    Create a Cartesian model based on the specified type.

    Parameters
    ----------
    model_type : str
        The type of model to create, either 'Structured' or 'Unstructured'.
    e : int or tuple of int
        Number of elements in each dimension (ex, ey, ez).
    h : int or tuple of float
        Element sizes in each dimension (hx, hy, hz). Required for structured models.
    l : tuple of float
        Domain sizes in each dimension (lx, ly, lz). Required for unstructured models.
    order : int
        The polynomial order of the elements.
    on_nodes : bool
        Whether to apply the model on nodes (True) or elements (False).
    elastic : bool
        Solving Elastic wave equation (True) or Acoustic wave equation (False).

    Returns
    -------
    model : Model.ModelStruct or Model.ModelUnstruct
        The created Cartesian model.

    Raises
    ------
    ValueError
        If the model type is unknown.
    """

    try:
        enum_value = ModelType[model_type]
    except KeyError:
        raise ValueError(f"Unknown python model type: {model_type}")
    match enum_value:
        case ModelType.STRUCTURED:
            return create_structured_model(e, l, order, on_nodes, elastic)
        case ModelType.UNSTRUCTURED:
            return create_unstructured_model(e, l, order, on_nodes, elastic)
        case _:
            raise ValueError(f"Unknown model type: {enum_value.name}")


def create_structured_model(e, l, order, on_nodes, elastic):
    """
    Create a structured Cartesian model based on the specified order.

    Parameters
    ----------
    e : int
        Number of elements in each dimension (ex, ey, ez).
    l : tuple of float
        Domain sizes in each dimension (lx, ly, lz).
    order : int
        The polynomial order of the elements.
    on_nodes: bool
        Whether to apply the model on nodes (True) or elements (False).
    elastic: bool
        Solving Elastic wave equation (True) or Acoustic wave equation (False).

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
            builder = CartesianStructBuilderFI1(
                e[0], l[0], e[1], l[1], e[2], l[2], on_nodes, elastic
            )
        case 2:
            builder = CartesianStructBuilderFI2(
                e[0], l[0], e[1], l[1], e[2], l[2], on_nodes, elastic
            )
        case 3:
            builder = CartesianStructBuilderFI3(
                e[0], l[0], e[1], l[1], e[2], l[2], on_nodes, elastic
            )
        case _:
            raise ValueError(
                f"Order {order} is not wrapped by pybind11 (only 1, 2, 3 supported)"
            )
    return builder.get_model()


def create_unstructured_model(e, l, order, on_nodes, elastic):
    """
    Create an unstructured Cartesian model.

    Parameters
    ----------
    e : tuple of int
        Number of elements in each dimension (ex, ey, ez).
    l : tuple of float
        Domain sizes in each dimension (lx, ly, lz).
    order : int
        The polynomial order of the elements.
    on_nodes: bool
        Whether to apply the model on nodes (True) or elements (False).
    elastic: bool
        Solving Elastic wave equation (True) or Acoustic wave equation (False).

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
    params = CartesianParams()
    params.ex, params.ey, params.ez = e
    params.lx, params.ly, params.lz = l
    params.order = order
    params.is_model_on_nodes = on_nodes
    params.is_elastic = elastic
    builder = CartesianUnstructBuilder(params)
    return builder.get_model()


def create_solver(implem_type, model_type, order, on_nodes, elastic):
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
    on_nodes : bool
        Whether the model is applied on nodes (True) or elements (False).
    elastic: bool
        Solving Elastic wave equation (True) or Acoustic wave equation (False).

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
    model_location = (
        Solver.ModelLocationType.ONNODES
        if on_nodes
        else Solver.ModelLocationType.ONELEMENTS
    )
    if elastic:
        physic_type = Solver.PhysicType.ELASTIC
    else:
        physic_type = Solver.PhysicType.ACOUSTIC

    return Solver.create_solver(
        Solver.MethodType.SEM, impl, model, model_location, physic_type, order
    )


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


def get_snapshot(nx, ny, nz, pnGlobal, normalize=False):
    """
    Extracts a 2D snapshot from a 3D global array at a specified index, with optional normalization.

    Parameters
    ----------
    nx : int
        Number of grid points in the x-direction.
    ny : int
        Number of grid points in the y-direction.
    nz : int
        Number of grid points in the z-direction.
    pnGlobal : np.ndarray
        The pressure field array.
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
        grid[i, j] = pnGlobal[I]

    if normalize:
        maxvalue = np.abs(grid).max()
        if maxvalue != 0:
            grid = grid / maxvalue

    return grid


def get_snapshot_param_model(nx, ny, nz, is_elastic, normalize=False):
    """
    Extracts a 2D snapshot from a 3D global array at a specified index, with optional normalization.

    Parameters
    ----------
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
        grid[i, j] = pnGlobal[I]

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
    im = ax.imshow(grid, cmap="viridis", interpolation="nearest", vmin=-cmpvalue, vmax=cmpvalue)
    plt.colorbar(im, ax=ax, label="Intensity")
    plt.title("2D Slice of a Float32 Array")
    plt.xlabel("X-axis")
    plt.ylabel("Z-axis")
    plt.ioff()  # Prevent showing the plot interactively
    return fig, ax, im


def plot_snapshot(nx, ny, nz, pnGlobal, im, t, prefix=""):
    """
    Plot and save a snapshot of the simulation at the given time step.

    Parameters
    ----------
    nx, ny, nz : int
        Grid dimensions.
    pnGlobal : np.ndarray
        Pressure field array.
    im : matplotlib.image.AxesImage
        The image object for updating the plot.
    t : int
        Current time step.
    prefix : string
        Prefix for snapshot name

    """

    grid = get_snapshot(nx, ny, nz, pnGlobal, False)
    im.set_array(grid)  # Update plot with new values
    plt.draw()  # Redraw the figure with updated data
    plt.ioff()
    plt.savefig(f"{prefix}snap0{t:0{5}d}.png")


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
    kk_pnGlobalPrev : kokkos array
        The Kokkos array for previous timestep pressure.
    pnGlobalPrev : np.ndarray
        The numpy array view of the previous timestep pressure.
    kk_pnGlobalCurr : kokkos array
        The Kokkos array for current timestep pressure.
    pnGlobalCurr : np.ndarray
        The numpy array view of the current timestep pressure.
    """

    kk_pnGlobalPrev = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pnGlobalPrev = np.array(kk_pnGlobalPrev, copy=False)
    pnGlobalPrev[:] = 0.0
    kk_pnGlobalCurr = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pnGlobalCurr = np.array(kk_pnGlobalCurr, copy=False)
    pnGlobalCurr[:] = 0.0

    return kk_pnGlobalPrev, pnGlobalPrev, kk_pnGlobalCurr, pnGlobalCurr


def allocate_displacement(n_dof, memspace, layout):
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

    kk_unGlobalPrev = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    unGlobalPrev = np.array(kk_unGlobalPrev, copy=False)
    unGlobalPrev[:] = 0.0
    kk_unGlobalCurr = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    unGlobalCurr = np.array(kk_unGlobalCurr, copy=False)
    unGlobalCurr[:] = 0.0

    return kk_unGlobalPrev, unGlobalPrev, kk_unGlobalCurr, unGlobalCurr


def allocate_rhs_term(n_rhs, n_time_steps, dt, f0, memspace, layout, src_file=None, backward=False):
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
    src_file : string
        path to ascii file.
    backward : bool
        Define source backward (starting from the end).
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
    RHSTerm[:,:] = 0.0
    if src_file is not None:
        if len(src_file)==n_rhs:
            for i in range(n_rhs):
                if os.path.exists(src_file[i]):
                    read_src=np.genfromtxt(src_file[i], delimiter="")
                    RHSTerm[i,:]=read_src[:,-1]
                else:
                    raise ValueError(
                        f"file {src_file[i]} does not exist"
                    )
        else:
            raise ValueError(
                f"src_file should be a list of size n_rhs:{n_rhs}"
            )
    else:
        src_file=[]
        for j in range(n_rhs):
            if os.path.exists(f"Src_{j}"):
                src_file.append(open(f"Src2_{j}",'w+'))
            else:
                src_file.append(open(f"Src_{j}",'w+'))
        if backward:
            for j in range(n_rhs):
                for i in range(n_time_steps-1):
                    ##To use same source as forward
                    #RHSTerm[j, i] = source_term((n_time_steps-i-2) * dt, f0)
                    ##To use different (amplitude and shifted) source as forward
                    RHSTerm[j, i] = 4.0 * source_term((n_time_steps-i-10-2) * dt, f0)
                    src_file[j].write(f"{i} {RHSTerm[j, i]}\n")
                src_file[j].write(f"{n_time_steps-1} 0.0\n") 
                RHSTerm[j, n_time_steps-1] = 0.0
        else:
            for j in range(n_rhs):
                src_file[j].write(f"0 0.0\n") 
                RHSTerm[j, 0] = 0.0
                for i in range(n_time_steps-1):
                    RHSTerm[j, i+1] = source_term(i * dt, f0)
                    src_file[j].write(f"{i+1} {RHSTerm[j, i+1]}\n")
        for j in range(n_rhs):
            src_file[j].close()
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
    RHSWeights[:, :] = 0.
    for i in range(n_rhs):
        ##Source taken on first node of element
        RHSWeights[i, :] = 0.
        RHSWeights[i, 0] = 1.
    return kk_RHSWeights, RHSWeights


def allocate_rcv_weight(n_rcv, model, memspace, layout):
    """
    Allocate and fill the RHSWeights array.

    Parameters
    ----------
    n_rcv : int
        Number of receivers.
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
    kk_RCVWeights = kokkos.array(
        [n_rcv, nb_points],
        dtype=kokkos.float32,
        space=memspace,
        layout=layout,
    )
    RCVWeights = np.array(kk_RCVWeights, copy=False)
    RCVWeights[:, :] = 0.
    for i in range(n_rcv):
        ##Receiver taken on first node of element
        RCVWeights[i, :] = 0.
        RCVWeights[i, 0] = 1.
            
    return kk_RCVWeights, RCVWeights


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
    RHSElement[:] =  0.
    ##For sources at 1/3 of first axis and middle of the 2 others
    for i in range(n_rhs) :
        RHSElement[i] =  ez / 3 + (i+1) * ey / (2*n_rhs) * ex +  ex / 2 * ey * ez

    return kk_RHSElement, RHSElement


def allocate_rcv_element(n_rcv, ex, ey, ez, memspace, layout):
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

    kk_RCVElement = kokkos.array(
        [n_rcv], dtype=kokkos.int32, space=memspace, layout=layout
    )
    RCVElement = np.array(kk_RCVElement, copy=False)
    RCVElement[:] =  0.
    ##For rcv at 2/3 of first axis and middle of the 2 others
    for i in range(n_rcv) :
        RCVElement[i] = 2 * ex / 3 + (i+1) * ey / (2*n_rcv) * ex + ez / 2 * ey * ex
    return kk_RCVElement, RCVElement


def create_solver_data_acoustic(kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RHSElement, kk_RHSWeights):
    """
    Create SEMsolverData instance and associated wavefield and rhs.

    Parameters
    ----------
    kk_RHSTerm : kokkos array
        The Kokkos array for the source term.
    kk_pnGlobalPrev : kokkos array
        The Kokkos array for previous pressure.
    kk_pnGlobalCurr : kokkos array
        The Kokkos array for current pressure.
    kk_RHSElement : kokkos array
        The Kokkos array for the element indices.
    kk_RHSWeights : kokkos array
        The Kokkos array for the weights.

    Returns
    -------
    wavefield : Solver.WavefieldAcoustic
        The Wavefield instance for acoustic propagation.
    rhs : Solver.RhsAcoustic
        The Rhs instance for acoustic propagation.
    data : Solver.SEMsolverData
        The SEMsolverData instance.
    """

    wavefield = Solver.WavefieldAcoustic(kk_pnGlobalPrev, kk_pnGlobalCurr)
    rhs = Solver.RhsAcoustic(kk_RHSTerm, kk_RHSElement, kk_RHSWeights)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)
    data.print()

    return wavefield, rhs, data


def create_solver_data_elastic(kk_RHSTermx, kk_RHSTermy, kk_RHSTermz, kk_uxnGlobalPrev, kk_uynGlobalPrev, kk_uznGlobalPrev, kk_uxnGlobalCurr, kk_uynGlobalCurr, kk_uznGlobalCurr, kk_RHSElement, kk_RHSWeights):
    """
    Create SEMsolverData instance and return it along with i1 and i2.

    Parameters
    ----------
    kk_RHSTermx : kokkos array
        The Kokkos array for the source term X.
    kk_RHSTermy : kokkos array
        The Kokkos array for the source term Y.
    kk_RHSTermz : kokkos array
        The Kokkos array for the source term Z.
    kk_uxnGlobal : kokkos array
        The Kokkos array for displacement X.
    kk_uynGlobal : kokkos array
        The Kokkos array for displacement Y.
    kk_uznGlobal : kokkos array
        The Kokkos array for displacement Z.
    kk_RHSElement : kokkos array
        The Kokkos array for the element indices.
    kk_RHSWeights : kokkos array
        The Kokkos array for the weights.

    Returns
    -------
    data : Solver.SEMsolverData
        The SEMsolverData instance.
    """

    wavefield = Solver.WavefieldElastic(kk_uxnGlobalPrev, kk_uxnGlobalCurr, kk_uynGlobalPrev, kk_uynGlobalCurr, kk_uznGlobalPrev, kk_uznGlobalCurr)
    rhs = Solver.RhsElastic(kk_RHSTermx, kk_RHSTermy, kk_RHSTermz, kk_RHSElement, kk_RHSWeights)
    data = Solver.SEMsolverDataElastic(wavefield, rhs)
    data.print()

    return wavefield, rhs, data


def compute_step_acoustic(
    time_sample,
    dt,
    solver,
    data,
    iteration_times,
    n_time_steps,
    nx,
    ny,
    nz,
    pnGlobalPrev,
    pnGlobalCurr,
    im,
    rcv_nodes,
    rcv_files,
    it,
    prefix="",
):
    """
    Perform a single time step of the simulation, update timing, plot, and swap indices.
    This now uses the split-phase approach (compute_forces + update_solution)
    which mimics the Domain Decomposition logic in C++.

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
    rcv : array of elem recv number

    Returns
    -------
    i1, i2 : int
        Updated indices for pressure fields.
    """

    iter_start = time.time()
    truedt=abs(dt)
    forward=(dt > 0)
    compute_one_step(solver, truedt, time_sample, data, forward)
    iter_time = time.time() - iter_start
    iteration_times.append(iter_time)
    if time_sample % 1000 == 0:
        print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
        print()
    if time_sample % 100 == 0:
        print(f"Time {time_sample} / {n_time_steps}")
        backtime_sample=time_sample if dt > 0 else (n_time_steps-time_sample)
        print(f"{backtime_sample} {pnGlobalPrev[rcv_nodes[0]]} {pnGlobalCurr[rcv_nodes[0]]}")
        for f in rcv_files: f.flush()
    ##plot output results in debug
    if time_sample % 100 == 0 and debug:
        plot_snapshot(nx, ny, nz, pnGlobalCurr, im, time_sample, prefix=prefix)
    if it % 2 == 1: # and time_sample > 0 and time_sample < n_time_steps:
        first=-1 if dt > 0 else 1
        ##write rcv files
        for i in range(len(rcv_nodes)):
            rcv_files[i].write(f"{time_sample+first} {pnGlobalPrev[rcv_nodes[i]]}\n")
            rcv_files[i].write(f"{time_sample} {pnGlobalCurr[rcv_nodes[i]]}\n")
            
    return


def compute_step_elastic(
    time_sample,
    dt,
    solver,
    data,
    iteration_times,
    n_time_steps,
    nx,
    ny,
    nz,
    uxnGlobalPrev,
    uynGlobalPrev,
    uznGlobalPrev,
    uxnGlobalCurr,
    uynGlobalCurr,
    uznGlobalCurr,
    im,
    rcv_nodes,
    rcv_files,
    it,
    prefix="",
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
    nx, ny, nz : int
        Grid dimensions for the plot.
    uxnGlobal : np.ndarray
        Displacement X field array.
    uynGlobal : np.ndarray
        Displacement Y field array.
    uznGlobal : np.ndarray
        Displacement Z field array.
    im : matplotlib.image.AxesImage
        The image object for updating the plot.

    """

    iter_start = time.time()
    truedt=abs(dt)
    forward=(dt > 0)
    compute_one_step(solver, truedt, time_sample, data, forward)
    iter_time = time.time() - iter_start
    iteration_times.append(iter_time)
    if time_sample % 1000 == 0:
        print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
        print()
    if time_sample % 100 == 0:
        print(f"Time {time_sample} / {n_time_steps}")
        backtime_sample=time_sample if dt > 0 else (n_time_steps-time_sample)
        print(f"{backtime_sample} {uxnGlobalPrev[rcv_nodes[0]]} {uxnGlobalCurr[rcv_nodes[0]]}")
        for f in rcv_files: f.flush()
    ##plot output results in debug
    if time_sample % 10 == 0 and debug:
        plot_snapshot(nx, ny, nz, uxnGlobalCurr, im, time_sample, f"{prefix}Ux_")
        plot_snapshot(nx, ny, nz, uynGlobalCurr, im, time_sample, f"{prefix}Uy_")
        plot_snapshot(nx, ny, nz, uznGlobalCurr, im, time_sample, f"{prefix}Uz_")
    if time_sample % 2 == 1:
        first=-1 if dt > 0 else 1
        ##write rcv files
        for i in range(len(rcv_nodes)):
            rcv_files[i*3  ].write(f"{time_sample+first} {uxnGlobalPrev[rcv_nodes[i]]}\n")
            rcv_files[i*3+1].write(f"{time_sample+first} {uynGlobalPrev[rcv_nodes[i]]}\n")
            rcv_files[i*3+2].write(f"{time_sample+first} {uznGlobalPrev[rcv_nodes[i]]}\n")
            rcv_files[i*3  ].write(f"{time_sample} {uxnGlobalCurr[rcv_nodes[i]]}\n")
            rcv_files[i*3+1].write(f"{time_sample} {uynGlobalCurr[rcv_nodes[i]]}\n")
            rcv_files[i*3+2].write(f"{time_sample} {uznGlobalCurr[rcv_nodes[i]]}\n")

    return


def compute_one_step(solver, dt, time_sample, data, forward):
    
    # 1. Compute forces (RHS of the equation)
    solver.compute_forces(dt, time_sample, data)

    # 2. Synchronize boundaries (Placeholder for distributed logic)
    # In C++, the BoundarySynchronizer is called here.
    # m_syncer->synchronize(m_solver->getForceVector(c), par_topology_);

    # 3. Update solution using mass matrix and accumulated forces
    solver.update_solution_forward(dt, data)

    return

    
def run_acoustic(model, solver, dt, n_time_steps, kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RHSElement, kk_RHSWeights, kk_RCVElement, kk_RCVWeights, is_backward, im, it_times, ndim, prefix):
    
    # Create solver data instance
    print("Creating solver data...")
    wavefield, rhs, data = create_solver_data_acoustic(kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RHSElement, kk_RHSWeights)
    print("Solver data created")

    if is_backward :
        dt=-dt
        time_list=list(range(n_time_steps-1,-1,-1))
    else :
        time_list=list(range(n_time_steps))

    print(f"min time : {min(time_list)} max time : {max(time_list)}")
        
    rcv_files=[]
    rcv_fname=[]
    rcv_nodes=[]

    # Loop over time steps
    print(f"nb receivers:{kk_RCVElement.shape[0]}")
    for i in range(kk_RCVElement.shape[0]):
        fname=f"{prefix}rcv_P{i}"
        filei=open(fname,'w+')
        rcv_files.append(filei)
        rcv_fname.append(fname)
        rcv_nodes.append(model.global_node_index(kk_RCVElement[i],0,0,0))

    it=0
    for time_sample in time_list:
        compute_step_acoustic(
            time_sample,
            dt,
            solver,
            data,
            it_times,
            n_time_steps,
            ndim[0],
            ndim[1],
            ndim[2],
            kk_pnGlobalPrev,
            kk_pnGlobalCurr,
            im,
            rcv_nodes,
            rcv_files,
            it,
            prefix=prefix
        )
        it+=1

        # Swap pressure arrays for next iteration
        data.swap_wavefields()


    for i in range(len(rcv_files)):
        rcv_files[i].close()

        
def run_elastic(model, solver, dt, n_time_steps, kk_RHSTermx, kk_RHSTermy, kk_RHSTermz, kk_uxnGlobalPrev, kk_uynGlobalPrev, kk_uznGlobalPrev, kk_uxnGlobalCurr, kk_uynGlobalCurr, kk_uznGlobalCurr, kk_RHSElement, kk_RHSWeights, kk_RCVElement, kk_RCVWeights, is_backward, im, it_times, ndim, prefix):
    # Create solver data instance
    print("Creating solver data...")
    wavefield, rhs, data = create_solver_data_elastic(kk_RHSTermx, kk_RHSTermy, kk_RHSTermz, kk_uxnGlobalPrev, kk_uynGlobalPrev, kk_uznGlobalPrev, kk_uxnGlobalCurr, kk_uynGlobalCurr, kk_uznGlobalCurr, kk_RHSElement, kk_RHSWeights)
    print("Solver data created")

    if is_backward :
        dt=-dt
        time_list=list(range(n_time_steps-1,-1,-1))
    else :
        time_list=list(range(n_time_steps))

    rcv_files=[]
    rcv_nodes=[]

    # Loop over time steps
    for i in range(kk_RCVElement.shape[0]):
        fileiX=open(f"{prefix}rcv_Ux{i}",'w+')
        fileiY=open(f"{prefix}rcv_Uy{i}",'w+')
        fileiZ=open(f"{prefix}rcv_Uz{i}",'w+')
        rcv_files.append(fileiX)
        rcv_files.append(fileiY)
        rcv_files.append(fileiZ)
        rcv_nodes.append(model.global_node_index(kk_RCVElement[i],0,0,0))

    it=0
    for time_sample in time_list:
        compute_step_elastic(
            time_sample,
            dt,
            solver,
            data,
            it_times,
            n_time_steps,
            ndim[0],
            ndim[1],
            ndim[2],
            kk_uxnGlobalPrev,
            kk_uynGlobalPrev,
            kk_uznGlobalPrev,
            kk_uxnGlobalCurr,
            kk_uynGlobalCurr,
            kk_uznGlobalCurr,
            im,
            rcv_nodes,
            rcv_files,
            it,
            prefix=prefix,
        )
        it+=1

        # Swap displacement arrays for next iteration
        data.swap_wavefields()


    for i in range(len(rcv_files)):
        rcv_files[i].close()
        

def test_self_adjoint(on_nodes,is_elastic,is_backward,f0,dt,n_time_steps,n_rhs,n_rcv,order,domain_size,ex,ey,ez,mem,impl,model):

    if mem=="default_mem":
        mem=detect_default_memspace()

    model_type=model

    lx = ly = lz = domain_size
    hx = lx / ex
    hy = ly / ey
    hz = lz / ez
    nx = ex * order + 1
    ny = ey * order + 1
    nz = ez * order + 1
    ndim=[ nx, ny, nz ]
    n_dof = nx * ny * nz
    n_elements = ex * ey * ez
    n_points_per_elements = (order + 1) * (order + 1) * (order + 1)
    
    print("==========SIMULATION PARAMETERS==========")
    print(f"order                        : {order}")
    print(f"on_nodes                     : {on_nodes}")
    print(f"memspace                     : {mem}")
    print(f"impl                         : {impl}")
    print(f"model                        : {model}")
    print(f"number of elements           : {n_elements}")
    print(f"number of points per element : {n_points_per_elements}")
    print(f"f0                           : {f0}")
    print(f"dt                           : {dt}")
    print(f"n_time_steps                 : {n_time_steps}")
    print(f"n_rhs                        : {n_rhs}")
    print(f"is_elastic                   : {is_elastic}")
    print("=========================================")

    # Setup graphic display
    print("Setting up plot...")
    if is_elastic:
        cmp_value=1.e-7
    else:
        cmp_value=0.1
    _, _, im = setup_plot(nx, nz, cmpvalue=cmp_value)
    print("Plot set up")
    
    # Initialize Kokkos
    kokkos.initialize()
    print("Kokkos initialized")
    memspace, layout = select_kokkos_memspace(mem)

    # Add timing variables
    start_time = time.time()
    simulation_start = datetime.now()
    iteration_times = []
    print(f"Simulation started at: {simulation_start}")

    # Create model
    print("Creating model...")
    model = create_model(
        model_type, (ex, ey, ez), (hx, hy, hz), (lx, ly, lz), order, on_nodes, is_elastic
    )
    print("Model created")

    # Create solver
    print("Creating solver...")
    solver = create_solver(impl, model_type, order, on_nodes, is_elastic)
    print("Solver created")

    # Initialize model
    print("Initializing model...")
    solver.compute_fe_init(model)
    print("Model initialized")

    # allocate state variable array
    if is_elastic:
        print("Allocating Displacement...")
        kk_uxnGlobalPrev, uxnGlobalPrev, kk_uxnGlobalCurr, uxnGlobalCurr = allocate_displacement(n_dof, memspace, layout)
        kk_uynGlobalPrev, uynGlobalPrev, kk_uynGlobalCurr, uynGlobalCurr = allocate_displacement(n_dof, memspace, layout)
        kk_uznGlobalPrev, uznGlobalPrev, kk_uznGlobalCurr, uznGlobalCurr = allocate_displacement(n_dof, memspace, layout)
        print("Displacement allocated")
    else:
        print("Allocating Pressure...")
        kk_pnGlobalPrev, pnGlobalPrev, kk_pnGlobalCurr, pnGlobalCurr = allocate_pressure(n_dof, memspace, layout)
        print("Pressure allocated")

    # allocate RHS arrays
    print("Allocating RHS element...")
    kk_RHSElement, RHSElement = allocate_rhs_element(
        n_rhs, ex, ey, ez, memspace, layout
    )
    for i in range(n_rhs):
        print(f"  - RHS element number i", RHSElement[i])
    print("RHS element allocated")

    print("Allocating RHS weights...")
    kk_RHSWeights, RHSWeights = allocate_rhs_weight(n_rhs, model, memspace, layout)
    print("RHS weights allocated")

    print("Allocating RHS term...")
    if is_elastic:
        kk_RHSTermx, rhsTerm = allocate_rhs_term(
            n_rhs, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
        kk_RHSTermy, rhsTerm = allocate_rhs_term(
            n_rhs, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
        kk_RHSTermz, rhsTerm = allocate_rhs_term(
            n_rhs, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
    else:
        kk_RHSTerm, rhsTerm = allocate_rhs_term(
            n_rhs, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
    print("RHS term allocated")

    # allocate RCV arrays
    print("Allocating RCV element...")
    kk_RCVElement, RCVElement = allocate_rcv_element(
        n_rcv, ex, ey, ez, memspace, layout
    )
    for i in range(n_rcv):
        print(f"  - RCV element number {i}", RCVElement[i])
    print("RCV element allocated")

    print("Allocating RCV weights...")
    kk_RCVWeights, RCVWeights = allocate_rcv_weight(n_rcv, model, memspace, layout)
    print("RCV weights allocated")

    #---------------------------------------------------
    # check physics parameters in debug
    if debug:
        if is_elastic:
            if(on_nodes) :
                for n in range(n_dof):
                    print(f"Vp  on node {n}: {model.get_model_vp_on_node(n)}")
                    print(f"Vs  on node {n}: {model.get_model_vs_on_node(n)}")
                    print(f"Rho on node {n}: {model.get_model_rho_on_node(n)}")
            else:
                for n in range(n_elements):
                    print(f"Vp  on node {n}: {model.get_model_vp_on_element(n)}")
                    print(f"Vs  on node {n}: {model.get_model_vs_on_element(n)}")
                    print(f"Rho on node {n}: {model.get_model_rho_on_element(n)}")
        else:
            if(on_nodes) :
                for n in range(n_dof):
                    print(f"Vp  on node {n}: {model.get_model_vp_on_node(n)}")
                    print(f"Rho on node {n}: {model.get_model_rho_on_node(n)}")
            else:
                for n in range(n_elements):
                    print(f"Vp  on node {n}: {model.get_model_vp_on_element(n)}")
                    print(f"Rho on node {n}: {model.get_model_rho_on_element(n)}")            
    #---------------------------------------------------

    prefix="shot1_"
    if is_elastic:
        run_elastic(model, solver, dt, n_time_steps, kk_RHSTermx, kk_RHSTermy, kk_RHSTermz, kk_uxnGlobalPrev, kk_uynGlobalPrev, kk_uznGlobalPrev, kk_uxnGlobalCurr, kk_uynGlobalCurr, kk_uznGlobalCurr, kk_RHSElement, kk_RHSWeights, kk_RCVElement, kk_RCVWeights, is_backward, im, iteration_times, ndim, prefix)
    else:
        run_acoustic(model, solver, dt, n_time_steps, kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RHSElement, kk_RHSWeights, kk_RCVElement, kk_RCVWeights, is_backward, im, iteration_times, ndim, prefix)
    
    is_backward=(not is_backward)
    if is_elastic:
        del kk_RHSTermx
        del kk_RHSTermy
        del kk_RHSTermz
        kk_RHSTermx, rhsTerm = allocate_rhs_term(
            n_rcv, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
        kk_RHSTermy, rhsTerm = allocate_rhs_term(
            n_rcv, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
        kk_RHSTermz, rhsTerm = allocate_rhs_term(
            n_rcv, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
       
    else:
        del kk_RHSTerm
        kk_RHSTerm, rhsTerm = allocate_rhs_term(
            n_rcv, n_time_steps, dt, f0, memspace, layout, backward=is_backward
        )
        
    prefix="shot2_"
    if is_elastic:
        uxnGlobalPrev[:] = 0.0
        uynGlobalPrev[:] = 0.0
        uznGlobalPrev[:] = 0.0
        uxnGlobalCurr[:] = 0.0
        uynGlobalCurr[:] = 0.0
        uznGlobalCurr[:] = 0.0
        run_elastic(model, solver, dt, n_time_steps, kk_RHSTermx, kk_RHSTermy, kk_RHSTermz, kk_uxnGlobalPrev, kk_uynGlobalPrev, kk_uznGlobalPrev, kk_uxnGlobalCurr, kk_uynGlobalCurr, kk_uznGlobalCurr, kk_RCVElement, kk_RCVWeights, kk_RHSElement, kk_RHSWeights, is_backward, im, iteration_times, ndim, prefix)
    else:
        pnGlobalPrev[:] = 0.0
        pnGlobalCurr[:] = 0.0
        run_acoustic(model, solver, dt, n_time_steps, kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RCVElement, kk_RCVWeights, kk_RHSElement, kk_RHSWeights, is_backward, im, iteration_times, ndim, prefix)

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
    
    prod_u1f2=0.
    prod_u2f1=0.
    f1_norm=0.
    u2_norm=0.
    if not is_elastic:
        u1=np.genfromtxt("shot1_rcv_P0", delimiter="")
        f2=np.genfromtxt("Src2_0", delimiter="")
        u2=np.genfromtxt("shot2_rcv_P0", delimiter="")
        f1=np.genfromtxt("Src_0", delimiter="")
        for i in range(n_time_steps):
            prod_u1f2+=u1[i,-1]*f2[i,-1]
            prod_u2f1+=u2[i,-1]*f1[n_time_steps-i-1,-1]
            if debug:
                print(f"step {i} {u1[i,-1]} {f2[i,-1]} {u2[i,-1]} {f1[n_time_steps-i-1,-1]}")
                print(f"step {i} {prod_u1f2} {prod_u2f1}")
            f1_norm+=f1[i,-1]*f1[i,-1]
            u2_norm+=u2[i,-1]*u2[i,-1]
    else:
        f1=np.genfromtxt("Src_0", delimiter="")
        f2=np.genfromtxt("Src2_0", delimiter="")
        u1=np.genfromtxt("shot1_rcv_Ux0", delimiter="")
        u2=np.genfromtxt("shot2_rcv_Ux0", delimiter="")
        for i in range(n_time_steps):
            prod_u1f2+=u1[i,-1]*f2[i,-1]
            prod_u2f1+=u2[i,-1]*f1[n_time_steps-i-1,-1]
            f1_norm+=f1[i,-1]*f1[i,-1]
            u2_norm+=u2[i,-1]*u2[i,-1]
        print(f"after Ux : {prod_u1f2} {prod_u2f1}")
        u1=np.genfromtxt("shot1_rcv_Uy0", delimiter="")
        u2=np.genfromtxt("shot2_rcv_Uy0", delimiter="")
        for i in range(n_time_steps):
            prod_u1f2+=u1[i,-1]*f2[i,-1]
            prod_u2f1+=u2[i,-1]*f1[n_time_steps-i-1,-1]
            f1_norm+=f1[i,-1]*f1[i,-1]
            u2_norm+=u2[i,-1]*u2[i,-1]
        print(f"after Uy : {prod_u1f2} {prod_u2f1}")
        u1=np.genfromtxt("shot1_rcv_Uz0", delimiter="")
        u2=np.genfromtxt("shot2_rcv_Uz0", delimiter="")
        for i in range(n_time_steps):
            prod_u1f2+=u1[i,-1]*f2[i,-1]
            prod_u2f1+=u2[i,-1]*f1[n_time_steps-i-1,-1]
            f1_norm+=f1[i,-1]*f1[i,-1]
            u2_norm+=u2[i,-1]*u2[i,-1]
        
            
    print(f"Final : {prod_u1f2} {prod_u2f1}")
    print(f"scalar products relative diff to ||f||.||ub|| ({np.sqrt(f1_norm)*np.sqrt(u2_norm)}): {abs(prod_u1f2 - prod_u2f1)/(np.sqrt(f1_norm)*np.sqrt(u2_norm))}")
    
    if not is_elastic:
        try:
            os.remove("shot1_rcv_P0")
        except OSError:
            pass
        try:
            os.remove("shot2_rcv_P0")
        except OSError:
            pass
    else:
        try:
            os.remove("shot1_rcv_Ux0")
        except OSError:
            pass
        try:
            os.remove("shot2_rcv_Ux0")
        except OSError:
            pass
        try:
            os.remove("shot1_rcv_Uy0")
        except OSError:
            pass
        try:
            os.remove("shot2_rcv_Uy0")
        except OSError:
            pass
        try:
            os.remove("shot1_rcv_Uz0")
        except OSError:
            pass
        try:
            os.remove("shot2_rcv_Uz0")
        except OSError:
            pass
        
    try:
        os.remove("Src2_0")
    except OSError:
        pass
    try:
        os.remove("Src_0")
    except OSError:
        pass
            
    if not is_elastic:
        tolerance=1e-6
    else:
        tolerance=1e-4
    assert abs(prod_u1f2 - prod_u2f1)/(np.sqrt(f1_norm)*np.sqrt(u2_norm)) < tolerance
    
    # release kokkos arrays and vectors
    if is_elastic:
        del uxnGlobalPrev
        del uynGlobalPrev
        del uznGlobalPrev
        del uxnGlobalCurr
        del uynGlobalCurr
        del uznGlobalCurr
        del kk_uxnGlobalPrev
        del kk_uynGlobalPrev
        del kk_uznGlobalPrev
        del kk_uxnGlobalCurr
        del kk_uynGlobalCurr
        del kk_uznGlobalCurr
        del kk_RHSTermx
        del kk_RHSTermy
        del kk_RHSTermz
    else:
        del pnGlobalPrev
        del pnGlobalCurr
        del kk_pnGlobalPrev
        del kk_pnGlobalCurr
        del kk_RHSTerm
    del rhsTerm
    del RHSElement
    del RHSWeights
    del RCVElement
    del RCVWeights
    del kk_RHSElement
    del kk_RHSWeights
    del kk_RCVElement
    del kk_RCVWeights
    del solver
    del model

    kokkos.finalize()
    print("End of  computation")

    if (abs(prod_u1f2 - prod_u2f1)/(np.sqrt(f1_norm)*np.sqrt(u2_norm)) > 1e-4):
        print(f"scalar products relative diff too big (> 1e-4): {abs(prod_u1f2 - prod_u2f1)/(np.sqrt(f1_norm)*np.sqrt(u2_norm))}")
        sys.exit(-1)

if __name__ == "__main__":
    main()
