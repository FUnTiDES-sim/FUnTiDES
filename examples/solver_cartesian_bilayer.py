#!/usr/bin/env python3
"""
This module runs the solver using a cartesian model with:
  - Kokkos GPU or CPU memory space
  - Structured or unstructured cartesian mesh
  - Polynomial order 1, 2 or 3
  - Implementation type: CLASSIC, MAKUTU
  - Only accoustic physics is implemented for now
It demonstrates usage of pybind11 wrapped C++ classes and functions from the proxys library.
For help run with par=<parameter_file> option.
"""

import os
import time
from datetime import datetime
from enum import Enum
import sys

import matplotlib.pyplot as plt
import numpy as np
import argparse

# FUnTiDES
import kokkos  # must be imported after matplotlib to avoid conflicts
import pyfuntides.model as Model
import pyfuntides.solver as Solver

# Create alias to use float32 and int32 types
CartesianUnstructBuilder = Model.CartesianUnstructBuilder_f32_i32
CartesianParams = Model.CartesianParams_f32_i32

ModelUnstruct = Model.ModelUnstruct_f32_i32
# Model data type alias
ModelUnstructData = Model.ModelUnstructData_f32_i32

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
    MAKUTU : str
        Makutu implementation.
    """

    CLASSIC = "Classic"
    MAKUTU = "MAKUTU"


def compute_gll_weights(src_coords_x, src_coords_y, src_coords_z, elem_to_node, coord_node, elem_indices, ngllx, nglly, ngllz):
    """
    Compute Lagrange interpolation weights at source coordinates for GLL nodes.
    
    Parameters
    ----------
    src_coords_x : array-like
        X coordinates of sources.
    src_coords_y : array-like
        Y coordinates of sources.
    src_coords_z : array-like
        Z coordinates of sources.
    elem_to_node : np.ndarray
        Element-to-node connectivity, shape (n_elem, n_gll_per_elem) for bilayer.
    coord_node : tuple of np.ndarray
        Node coordinates (x_coords, y_coords, z_coords).
    elem_indices : array-like
        Element indices (0-based) for each source.
    ngllx, nglly, ngllz : int
        Number of GLL points in each direction.
    
    Returns
    -------
    weights : np.ndarray
        Lagrange weights, shape (n_rhs, ngllx*nglly*ngllz).
    """
    n_rhs = len(src_coords_x)
    n_gll = ngllx * nglly * ngllz
    weights = np.zeros((n_rhs, n_gll))
    
    for rhs_idx in range(n_rhs):
        src_x = src_coords_x[rhs_idx]
        src_y = src_coords_y[rhs_idx]
        src_z = src_coords_z[rhs_idx]
        elem_idx = elem_indices[rhs_idx]
        
        # Get GLL node coordinates for this element
        gll_coords_x = np.zeros((ngllx, nglly, ngllz))
        gll_coords_y = np.zeros((ngllx, nglly, ngllz))
        gll_coords_z = np.zeros((ngllx, nglly, ngllz))
        
        # Get node indices for this element
        node_indices = elem_to_node[elem_idx, :]
        
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    # Local index within element (C-style: i varies fastest)
                    local_idx = i + j * ngllx + k * ngllx * nglly
                    node_idx = node_indices[local_idx]
                    gll_coords_x[i, j, k] = coord_node[0][node_idx]
                    gll_coords_y[i, j, k] = coord_node[1][node_idx]
                    gll_coords_z[i, j, k] = coord_node[2][node_idx]
        
        # Extract 1D GLL coordinates along each axis
        xi_nodes = gll_coords_x[:, 0, 0]   # x-coordinates at (i, 0, 0)
        eta_nodes = gll_coords_y[0, :, 0]  # y-coordinates at (0, j, 0)
        zeta_nodes = gll_coords_z[0, 0, :]  # z-coordinates at (0, 0, k)
        
        # Compute 1D Lagrange polynomials for each direction
        def lagrange_1d(x, x_nodes):
            """Evaluate all Lagrange basis functions at point x."""
            n = len(x_nodes)
            weights = np.zeros(n)
            for i in range(n):
                weight = 1.0
                for j in range(n):
                    if i != j:
                        weight *= (x - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
                weights[i] = weight
            return weights
        
        # Evaluate 1D weights
        w_xi = lagrange_1d(src_x, xi_nodes)
        w_eta = lagrange_1d(src_y, eta_nodes)
        w_zeta = lagrange_1d(src_z, zeta_nodes)
        
        # Compute 3D tensor product weights for this source
        idx = 0
        for k in range(ngllz):
            for j in range(nglly):
                for i in range(ngllx):
                    weights[rhs_idx, idx] = w_xi[i] * w_eta[j] * w_zeta[k]
                    idx += 1
    
    return weights


def find_source_elements(src_x, src_y, src_z, elem_to_node, coord_node, ngllx, nglly, ngllz):
    """
    Find the element indices containing the given source coordinates for bilayer mesh.
    
    Parameters
    ----------
    src_x : array-like
        X-coordinates of sources.
    src_y : array-like
        Y-coordinates of sources.
    src_z : array-like
        Z-coordinates of sources.
    elem_to_node : np.ndarray
        Element-to-node connectivity array, shape (n_elem, n_gll_per_elem).
    coord_node : tuple of np.ndarray
        Tuple of (x_coords, y_coords, z_coords) for all nodes.
    ngllx : int
        Number of GLL points in x direction.
    nglly : int
        Number of GLL points in y direction.
    ngllz : int
        Number of GLL points in z direction.
    
    Returns
    -------
    elem_indices : np.ndarray
        Array of element indices (0-based for solver).
    """
    n_sources = len(src_x)
    n_elem = elem_to_node.shape[0]
    elem_indices = np.zeros(n_sources, dtype=np.int32)
    
    # For structured bilayer mesh
    for src_idx in range(n_sources):
        x, y, z = src_x[src_idx], src_y[src_idx], src_z[src_idx]
        
        # Search through elements
        found = False
        for elem_idx in range(n_elem):
            # Get node indices for this element
            node_indices = elem_to_node[elem_idx, :]
            
            # Get bounding box from corner nodes (first and last in each direction)
            # For a hexahedral element with (order+1)^3 nodes
            corners = [
                node_indices[0],  # (0,0,0)
                node_indices[ngllx-1],  # (ngllx-1, 0, 0)
                node_indices[ngllx * (nglly-1)],  # (0, nglly-1, 0)
                node_indices[ngllx * nglly - 1],  # (ngllx-1, nglly-1, 0)
                node_indices[ngllx * nglly * (ngllz-1)],  # (0, 0, ngllz-1)
                node_indices[ngllx * nglly * (ngllz-1) + ngllx - 1],  # (ngllx-1, 0, ngllz-1)
                node_indices[ngllx * nglly * (ngllz-1) + ngllx * (nglly-1)],  # (0, nglly-1, ngllz-1)
                node_indices[-1],  # (ngllx-1, nglly-1, ngllz-1)
            ]
            
            # Get bounding box from corners
            x_corners = [coord_node[0][c] for c in corners]
            y_corners = [coord_node[1][c] for c in corners]
            z_corners = [coord_node[2][c] for c in corners]
            
            x_min, x_max = min(x_corners), max(x_corners)
            y_min, y_max = min(y_corners), max(y_corners)
            z_min, z_max = min(z_corners), max(z_corners)
            
            # Check if point is within bounding box (with small tolerance)
            tol = 1e-3
            if (x_min - tol <= x <= x_max + tol and
                y_min - tol <= y <= y_max + tol and
                z_min - tol <= z <= z_max + tol):
                elem_indices[src_idx] = elem_idx  # 0-based indexing for C++ solver
                found = True
                break
        
        if not found:
            print(f"Warning: No element found for source {src_idx} at ({x}, {y}, {z})")
            # Use a default element (center of domain)
            elem_indices[src_idx] = n_elem // 2
    
    return elem_indices


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
    if enum_value == ModelType.STRUCTURED:
        return Solver.MeshType.STRUCT
    elif enum_value == ModelType.UNSTRUCTURED:
        return Solver.MeshType.UNSTRUCT
    else:
        raise ValueError(f"Unknown solver model type for: {enum_value.name}")


def get_solver_implem_type(implem_type):
    """
    Map an implementation identifier (name or enum) to the corresponding Solver.ImplemType.

    Parameters
    ----------
    implem_type : str or ImplemType
        Implementation name or ImplemType enum. Accepted names are
        'CLASSIC', 'MAKUTU' (case-insensitive when passed as enum names).

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
    if enum_value == ImplemType.CLASSIC:
        return Solver.ImplemType.CLASSIC
    elif enum_value == ImplemType.MAKUTU:
        return Solver.ImplemType.MAKUTU
    else:
        raise ValueError(
            f"Unknown solver implementation type for: {enum_value.name}"
        )


def create_bilayer_model_data_from_params(order, n_elems, L, vels, z_elem_interface, memspace, layout):
    """
    Create model data from user parameters.

    !!! Only bilayer, model on elements, and acoustic are supported for now !!!

    Parameters
    ----------
    order : int
        The polynomial order of the elements.
    n_elems : tuple of int
        Number of elements in each dimension of the mesh.
    L : tuple of float
        Domain sizes in each dimension (lx, ly, lz).
    vels : tuple of float
        Velocities for the two layers (vel_layer1, vel_layer2).
    z_elem_interface : int
        The index of the element interface in the z-direction.
    memspace : kokkos.Space
        The Kokkos memory space.
    layout : kokkos.Layout
        The Kokkos layout.
    """

    lx = L[0]
    ly = L[1]
    lz = L[2]

    ex = n_elems[0]
    ey = n_elems[1]
    ez = n_elems[2]

    n_elem = ex * ey * ez

    # Create 3D squared elements mesh with bilayer velocity model
    mesh_model_vp = np.zeros((n_elem,), dtype=np.float32)
    mesh_model_rho = np.ones((n_elem,), dtype=np.float32) # always 1.0 for acoustic

    for k in range(ez):
        for j in range(ey):
            for i in range(ex):
                e = i + j * ex + k * ex * ey
                if k < z_elem_interface:
                    mesh_model_vp[e] = vels[0]
                else:
                    mesh_model_vp[e] = vels[1]

    # Create elem_to_nodes and nodes_coords arrays
    # Total number of nodes in each direction
    nx_nodes = order * ex + 1
    ny_nodes = order * ey + 1
    nz_nodes = order * ez + 1
    n_node = nx_nodes * ny_nodes * nz_nodes

    # Node spacing
    dx_node = lx / (order * ex)
    dy_node = ly / (order * ey)
    dz_node = lz / (order * ez)

    # Create coordinate arrays
    nodes_coords_x = np.zeros(n_node, dtype=np.float32)
    nodes_coords_y = np.zeros(n_node, dtype=np.float32)
    nodes_coords_z = np.zeros(n_node, dtype=np.float32)

    # Fill coordinates (C-style: x varies fastest)
    for k in range(nz_nodes):
        for j in range(ny_nodes):
            for i in range(nx_nodes):
                idx = i + j * nx_nodes + k * nx_nodes * ny_nodes
                nodes_coords_x[idx] = i * dx_node
                nodes_coords_y[idx] = j * dy_node
                nodes_coords_z[idx] = k * dz_node


    n_gll_per_elem = (order + 1) * (order + 1) * (order + 1)
    elem_to_nodes = np.zeros((n_elem, n_gll_per_elem), dtype=np.int32)

    # Build connectivity for each element
    for ez_i in range(ez):
        for ey_i in range(ey):
            for ex_i in range(ex):
                e = ex_i + ey_i * ex + ez_i * ex * ey  # Element index
            
                # For each GLL node within the element
                for k_loc in range(order + 1):
                    for j_loc in range(order + 1):
                        for i_loc in range(order + 1):
                            # Global node indices
                            i_glob = ex_i * order + i_loc
                            j_glob = ey_i * order + j_loc
                            k_glob = ez_i * order + k_loc
                        
                            # Global node index (C-style)
                            node_idx = i_glob + j_glob * nx_nodes + k_glob * nx_nodes * ny_nodes
                        
                            # Local index within element
                            local_idx = i_loc + j_loc * (order + 1) + k_loc * (order + 1) * (order + 1)
                        
                            elem_to_nodes[e, local_idx] = node_idx

    # Create Kokkos views
    kk_elem_to_nodes = kokkos.array(elem_to_nodes, dtype=kokkos.int32, space=memspace, layout=layout)

    # 1D Kokkos views for coordinates
    kk_nodes_coords_x = kokkos.array([n_node], dtype=kokkos.float32, space=memspace, layout=layout)
    kk_nodes_coords_y = kokkos.array([n_node], dtype=kokkos.float32, space=memspace, layout=layout)
    kk_nodes_coords_z = kokkos.array([n_node], dtype=kokkos.float32, space=memspace, layout=layout)
    np.array(kk_nodes_coords_x, copy=False)[:] = np.array(nodes_coords_x, dtype=np.float32)
    np.array(kk_nodes_coords_y, copy=False)[:] = np.array(nodes_coords_y, dtype=np.float32)
    np.array(kk_nodes_coords_z, copy=False)[:] = np.array(nodes_coords_z, dtype=np.float32)

    print(f"Diagnostics for model creation:")
    print(f" nodes_coords_x: shape={np.array(kk_nodes_coords_x, copy=False).shape}, dtype={np.array(kk_nodes_coords_x, copy=False).dtype}")
    print(f" nodes_coords_y: shape={np.array(kk_nodes_coords_y, copy=False).shape}, dtype={np.array(kk_nodes_coords_y, copy=False).dtype}")
    print(f" nodes_coords_z: shape={np.array(kk_nodes_coords_z, copy=False).shape}, dtype={np.array(kk_nodes_coords_z, copy=False).dtype}")
    print(f"  nodes_coords_x: min={np.array(kk_nodes_coords_x, copy=False).min():.2f}, max={np.array(kk_nodes_coords_x, copy=False).max():.2f}")
    print(f"  nodes_coords_y: min={np.array(kk_nodes_coords_y, copy=False).min():.2f}, max={np.array(kk_nodes_coords_y, copy=False).max():.2f}")
    print(f"  nodes_coords_z: min={np.array(kk_nodes_coords_z, copy=False).min():.2f}, max={np.array(kk_nodes_coords_z, copy=False).max():.2f}")
    print(f" elem_to_nodes: shape={elem_to_nodes.shape}, dtype={elem_to_nodes.dtype}")
    print(f"  elem_to_nodes: min={elem_to_nodes.min()}, max={elem_to_nodes.max()}")

    # 1D Kokkos views for model properties
    kk_model_vp_element = kokkos.array([n_elem], dtype=kokkos.float32, space=memspace, layout=layout)
    kk_model_rho_element = kokkos.array([n_elem], dtype=kokkos.float32, space=memspace, layout=layout)
    np.array(kk_model_vp_element, copy=False)[:] = np.array(mesh_model_vp, dtype=np.float32)
    np.array(kk_model_rho_element, copy=False)[:] = np.array(mesh_model_rho, dtype=np.float32)
    # kk_model_vs_element set to empty for acoustic tests
    kk_model_vs_element = kokkos.array([1], dtype=kokkos.float32, space=memspace, layout=layout)

    # Check for invalid model values
    vp_element_np = np.array(kk_model_vp_element, copy=False)
    rho_element_np = np.array(kk_model_rho_element, copy=False)
    print(f"\nModel diagnostics:")
    print(f" vp_element: shape={vp_element_np.shape}, dtype={vp_element_np.dtype}")
    print(f" rho_element: shape={rho_element_np.shape}, dtype={rho_element_np.dtype}")
    print(f"  vp_element: min={vp_element_np.min():.2f}, max={vp_element_np.max():.2f}, has_zeros={np.any(vp_element_np == 0)}, has_nan={np.isnan(vp_element_np).any()}")
    print(f"  rho_element: min={rho_element_np.min():.2f}, max={rho_element_np.max():.2f}, has_zeros={np.any(rho_element_np == 0)}, has_nan={np.isnan(rho_element_np).any()}")
    if np.any(vp_element_np <= 0) or np.any(rho_element_np <= 0):
        print("  WARNING: Model has zero or negative values!")

    # Create empty Kokkos views for node-based properties (set to None if not used)
    kk_model_vp_node = kokkos.array([1], dtype=kokkos.float32, space=memspace, layout=layout)
    kk_model_rho_node = kokkos.array([1], dtype=kokkos.float32, space=memspace, layout=layout)
    kk_model_vs_node = kokkos.array([1], dtype=kokkos.float32, space=memspace, layout=layout)

    # Empty views for optional parameters
    kk_empty_1d_float = kokkos.array([1], dtype=kokkos.float32, space=memspace, layout=layout)
    kk_empty_1d_int = kokkos.array([0], dtype=kokkos.int32, space=memspace, layout=layout)
    kk_empty_3d = kokkos.array([1, 1, 1], dtype=kokkos.float32, space=memspace, layout=layout)

    model_data = ModelUnstructData(
        order,
        n_elem,
        n_node,
        lx,
        ly,
        lz,
        True,
        False,
        kk_elem_to_nodes,
        kk_nodes_coords_x,
        kk_nodes_coords_y,
        kk_nodes_coords_z,
        kk_model_vp_node,
        kk_model_vp_element,
        kk_model_rho_node,
        kk_model_rho_element,
        kk_model_vs_node,
        kk_model_vs_element,
        kk_empty_1d_float,  # model_delta_node
        kk_empty_1d_float,  # model_delta_element
        kk_empty_1d_float,  # model_epsilon_node
        kk_empty_1d_float,  # model_epsilon_element
        kk_empty_1d_float,  # model_gamma_node
        kk_empty_1d_float,  # model_gamma_element
        kk_empty_1d_float,  # model_theta_node
        kk_empty_1d_float,  # model_theta_element
        kk_empty_1d_float,  # model_phi_node
        kk_empty_1d_float,  # model_phi_element
        kk_empty_3d,  # model_C_tensor_element
        kk_empty_1d_int,  # boundaries_t (MUST be int32, not float32, lenght 0 for now since we are not testing boundaries in this example)
    )
    # Create ModelUnstruct from ModelUnstructData
    # This mimics what CartesianUnstructBuilder::getModel() does in C++
    model_obj = ModelUnstruct(model_data)
    print("ModelUnstruct created")
    
    # Build face connectivity - CRITICAL: initializes internal mesh structures
    # This mimics what CartesianUnstructBuilder::getModel() calls after creating ModelUnstruct
    print("  Building face connectivity...")
    model_obj.build_face_connectivity()
    print("  Face connectivity built")
    
    return model_obj, elem_to_nodes, (nodes_coords_x, nodes_coords_y, nodes_coords_z)


def create_solver(implem_type, model_type, order, on_nodes):
    """
    Create a solver based on the specified implementation type.

    Parameters
    ----------
    implem_type : str
        The implementation type, one of 'CLASSIC' or 'MAKUTU'.
    model_type : str
        The model type, either 'Structured' or 'Unstructured'.
    order : int
        The polynomial order of the elements.
    on_nodes : bool
        Whether the model is applied on nodes (True) or elements (False).

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
    # We are running an acoustic simulation in this example
    physic_type = Solver.PhysicType.ACOUSTIC

    return Solver.create_solver(
        Solver.MethodType.SEM, impl, model, model_location, Solver.PhysicType.ACOUSTIC, order
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

    # Extract the middle y-slice: y_index = ny // 2
    # Data layout: index = ix + iy * nx + iz * nx * ny (C-style row-major)
    iy = ny // 2
    grid = np.zeros((nx, nz))
    for iz in range(nz):
        for ix in range(nx):
            index = ix + iy * nx + iz * nx * ny
            grid[ix, iz] = pnGlobal[index]

    if normalize:
        maxvalue = np.abs(grid).max()
        if maxvalue != 0:
            grid = grid / maxvalue

    return grid


def setup_plot(nx, nz, cmpvalue=0.15, lx=None, lz=None):
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
    lx : float, optional
        Domain size in the x-direction (default None).
    lz : float, optional
        Domain size in the z-direction (default None).

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
    
    if lx is not None and lz is not None:
        im = ax.imshow(grid.T, cmap="viridis", interpolation="nearest", extent=[0, lx, 0, lz])
        plt.xlabel("X-axis (m)")
        plt.ylabel("Z-axis (m)")
    else:
        im = ax.imshow(grid.T, cmap="viridis", interpolation="nearest")
        plt.xlabel("X-axis")
        plt.ylabel("Z-axis")
    
    plt.colorbar(im, ax=ax, label="Intensity")
    plt.title(f"2D Slice of Pressure Array. Max value: {np.max(grid):.2e}, Min value: {np.min(grid):.2e}")
    plt.ioff()  # Prevent showing the plot interactively
    return fig, ax, im


def plot_snapshot(nx, ny, nz, pnGlobal, im, t, lx=None, lz=None):
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
    lx : float, optional
        Domain size in the x-direction (default None).
    lz : float, optional
        Domain size in the z-direction (default None).
    """

    grid = get_snapshot(nx, ny, nz, pnGlobal, False)
    if lz is not None and lx is not None:
        im.set_extent([0, lx, 0, lz])
        ax = im.axes
        ax.set_xlim(0, lx)
        ax.set_ylim(0, lz)
    im.set_array(grid.T)  # Update plot with new values
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
        for j in range(n_rhs):
            RHSTerm[j, i] = source_term(i * dt, f0)
    return kk_RHSTerm, RHSTerm


def allocate_rhs_weight(n_rhs, model, memspace, layout, user_weights=None):
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
    user_weights : np.ndarray, optional
        User-provided projection weights. If None, uniform weights are used.

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

    if user_weights is not None:
        # Use user-input projection weights
        assert n_rhs == user_weights.shape[0], "Number of RHS does not match user weights"
        for i in range(n_rhs):
            RHSWeights[i, :] = user_weights[i, :]
    else:
        # Uniform weights: source in centre of element
        for i in range(n_rhs):
            for j in range(nb_points):
                RHSWeights[i, j] = 1 / nb_points

    return kk_RHSWeights, RHSWeights


def allocate_rhs_element(n_rhs, ex, ey, ez, memspace, layout, src_elem_indices=None):
    """
    Allocate and fill the RHSElement array.

    Parameters
    ----------
    n_rhs : int
        Number of right-hand side sources.
    ex : int
        Number of elements in the x-direction.
    ey : int
        Number of elements in the y-direction.
    ez : int
        Number of elements in the z-direction.
    memspace : kokkos.Space
        The selected Kokkos memory space.
    layout : kokkos.Layout
        The selected Kokkos layout.
    src_elem_indices : np.ndarray, optional
        Array of source element indices (0-based for solver). If None, default positions are used.

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

    if src_elem_indices is not None:
        # Use user-input source element indices
        assert len(src_elem_indices) == n_rhs, "Number of source elements does not match n_rhs"
        for i in range(n_rhs):
            RHSElement[i] = int(src_elem_indices[i])
    else:
        # Place sources near center: ix, iy, iz -> element_index = ix + iy*ex + iz*ey*ex
        # Using ey//2 ensures snapshot slice at ny/2 will show the sources
        RHSElement[0] = ex // 2 + (ey // 2) * ex + (ez // 2) * ey * ex
        RHSElement[1] = ex // 3 + (ey // 2) * ex + (ez // 2) * ey * ex
    
    return kk_RHSElement, RHSElement


def create_solver_data(kk_RHSTerm, kk_pnGlobalPrev, kk_pnGlobalCurr, kk_RHSElement, kk_RHSWeights):
    """
    Create SEMsolverData instance and associated wavefield and rhs.

    Parameters
    ----------
    kk_RHSTerm : kokkos array
        The Kokkos array for the source term.
    kk_pnGlobalPrev : kokkos array
        The Kokkos array for previous timestep pressure.
    kk_pnGlobalCurr : kokkos array
        The Kokkos array for current timestep pressure.
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
    data : Solver.SEMsolverDataAcoustic
        The SEMsolverData instance for acoustic propagation.
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
    pnGlobal,
    im,
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
    data : Solver.SEMsolverDataAcoustic
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
        Pressure field array (for python plots).
    im : matplotlib.image.AxesImage
        The image object for updating the plot.

    Returns
    -------
    i1, i2 : int
        Updated indices for pressure fields.
    """

    iter_start = time.time()

    # 1. Compute forces (RHS of the equation)
    solver.compute_forces(dt, time_sample, data)

    # 2. Synchronize boundaries (Placeholder for distributed logic)
    # In C++, the BoundarySynchronizer is called here.
    # m_syncer->synchronize(m_solver->getForceVector(c), par_topology_);

    # 3. Update solution using mass matrix and accumulated forces
    solver.update_solution_forward(dt, data)

    iter_time = time.time() - iter_start
    iteration_times.append(iter_time)

    if time_sample % 1000 == 0:
        print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
        print()
    if time_sample % 100 == 0:
        print(f"Time {time_sample} / {n_time_steps}")
    if time_sample % 100 == 0:
        plot_snapshot(nx, ny, nz, pnGlobal, im, time_sample)
        print(f"Max pressure at time {time_sample}: {np.max(np.abs(pnGlobal))}",flush=True)
        print(f"Min pressure at time {time_sample}: {np.min(np.abs(pnGlobal))}",flush=True)
        print(f"Index of max pressure at time {time_sample}: {np.unravel_index(np.argmax(np.abs(pnGlobal)), pnGlobal.shape)}",flush=True)


def parse_param_file(filename):
    """
    Parse a parameter file with key=value pairs.
    
    Parameters
    ----------
    filename : str
        Path to the parameter file.
    
    Returns
    -------
    dict
        Dictionary of parameter names and values.
    """
    params = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            # Parse key=value
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                # Remove quotes if present
                if value.startswith('"') and value.endswith('"'):
                    value = value[1:-1]
                elif value.startswith("'") and value.endswith("'"):
                    value = value[1:-1]
                params[key] = value
    return params


def main():

    # Auto-detect default memory space
    default_mem = detect_default_memspace()
    
    # Initialize argparse parser
    parser = argparse.ArgumentParser(
        description="FUnTiDES Cartesian Bilayer Solver - Run SEM solver with bilayer velocity model"
    )
    
    # Parameter file option
    parser.add_argument("--param_file", type=str, default=None,
                       help="Path to parameter file with key=value pairs")
    
    # Solver parameters
    parser.add_argument("--mem", type=str, default=default_mem,
                       help="Kokkos memspace: CPU or GPU")
    parser.add_argument("--model", type=str, default=ModelType.UNSTRUCTURED.name,
                       help="Model type: STRUCTURED or UNSTRUCTURED (default: UNSTRUCTURED for bilayer)")
    parser.add_argument("--impl", type=str, default=ImplemType.MAKUTU.name,
                       help="Implementation type: CLASSIC or MAKUTU")
    parser.add_argument("--order", type=int, default=2,
                       help="Polynomial order of the elements (1, 2, or 3)")
    parser.add_argument("--domain_sizeX", type=float, default=3500.0,
                       help="Size of the domain in X direction")
    parser.add_argument("--domain_sizeY", type=float, default=3500.0,
                       help="Size of the domain in Y direction")
    parser.add_argument("--domain_sizeZ", type=float, default=1500.0,
                       help="Size of the domain in Z direction")
    parser.add_argument("--ex", type=int, default=350,
                       help="Number of elements in x-direction")
    parser.add_argument("--ey", type=int, default=350,
                       help="Number of elements in y-direction")
    parser.add_argument("--ez", type=int, default=150,
                       help="Number of elements in z-direction")
    parser.add_argument("--f0", type=float, default=5.0,
                       help="Peak frequency for the Ricker source term")
    parser.add_argument("--dt", type=float, default=0.001,
                       help="Time step size")
    parser.add_argument("--n_time_steps", type=int, default=1500,
                       help="Number of time steps to run")
    parser.add_argument("--n_rhs", type=int, default=2,
                       help="Number of right-hand side sources")
    parser.add_argument("--on_nodes", action="store_true",
                       help="Whether to apply model on nodes")
    parser.add_argument("--boundaries_size", type=float, default=0.0,
                       help="Size of absorbing boundaries in meters")
    parser.add_argument("--surface_sponge", action="store_true",
                       help="Enable sponge at the free surface")
    parser.add_argument("--taper_delta", type=float, default=0.015,
                       help="Taper delta for sponge boundaries")
    parser.add_argument("--vel1", type=float, default=1500.0,
                       help="Velocity of the first layer in m/s")
    parser.add_argument("--vel2", type=float, default=2500.0,
                       help="Velocity of the second layer in m/s")
    parser.add_argument("--z_interface", type=int, default=None,
                       help="Depth index of the interface (default: ez/2)")
    parser.add_argument("--src_x", type=float, default=None,
                       help="X coordinate of the source (default: domain_sizeX/2)")
    parser.add_argument("--src_y", type=float, default=None,
                       help="Y coordinate of the source (default: domain_sizeY/2)")
    parser.add_argument("--src_z", type=float, default=None,
                       help="Z coordinate of the source (default: domain_sizeZ/2)")
    
    # Parse arguments
    args = parser.parse_args()
    
    # If param_file is provided, read parameters from it and update args
    if args.param_file:
        file_params = parse_param_file(args.param_file)
        # Convert file parameters to appropriate types and update args
        for key, value in file_params.items():
            if hasattr(args, key):
                # Get the current argument type
                current_val = getattr(args, key)
                if isinstance(current_val, bool):
                    # Handle boolean values
                    setattr(args, key, value.lower() in ('true', '1', 'yes', 'on'))
                elif isinstance(current_val, int):
                    setattr(args, key, int(value))
                elif isinstance(current_val, float):
                    setattr(args, key, float(value))
                elif current_val is None:
                    # For None defaults, try to infer type from parameter name patterns
                    if key in ['z_interface', 'ex', 'ey', 'ez', 'order', 'n_rhs', 'n_time_steps']:
                        setattr(args, key, int(value))
                    elif key in ['src_x', 'src_y', 'src_z', 'domain_sizeX', 'domain_sizeY', 'domain_sizeZ', 
                                'f0', 'dt', 'boundaries_size', 'taper_delta', 'vel1', 'vel2']:
                        setattr(args, key, float(value))
                    else:
                        setattr(args, key, value)
                else:
                    setattr(args, key, value)
    
    # Extract parameters
    mem = args.mem
    model = args.model
    impl = args.impl
    order = args.order
    domain_sizeX = args.domain_sizeX
    domain_sizeY = args.domain_sizeY
    domain_sizeZ = args.domain_sizeZ
    ex = args.ex
    ey = args.ey
    ez = args.ez
    f0 = args.f0
    dt = args.dt
    n_time_steps = args.n_time_steps
    n_rhs = args.n_rhs
    on_nodes = args.on_nodes
    boundaries_size = args.boundaries_size
    surface_sponge = args.surface_sponge
    taper_delta = args.taper_delta
    vel1 = args.vel1
    vel2 = args.vel2
    
    # Set defaults that depend on other parameters
    z_interface = args.z_interface if args.z_interface is not None else int(ez/2)
    src_info_x = args.src_x if args.src_x is not None else domain_sizeX/2
    src_info_y = args.src_y if args.src_y is not None else domain_sizeY/2
    src_info_z = args.src_z if args.src_z is not None else domain_sizeZ/2

    if z_interface < 0 or z_interface >= ez:
        raise ValueError("z_interface must be between 0 and ez-1")

    # Print parameters
    print("Parsing parameters...")
    print(f" memspace          : {mem}")
    print(f" model             : {model}")
    print(f" implementation    : {impl}")
    print(f" order             : {order}")
    print(f" domain_sizeX     : {domain_sizeX}")
    print(f" domain_sizeY     : {domain_sizeY}")
    print(f" domain_sizeZ     : {domain_sizeZ}")
    print(f" ex                : {ex}")
    print(f" ey                : {ey}")
    print(f" ez                : {ez}")
    print(f" f0                : {f0}")
    print(f" dt                : {dt}")
    print(f" n_time_steps      : {n_time_steps}")
    print(f" n_rhs             : {n_rhs}")
    print(f" on_nodes          : {on_nodes}")
    print(f" boundaries_size   : {boundaries_size}")
    print(f" surface_sponge    : {surface_sponge}")
    print(f" taper_delta       : {taper_delta}")
    print(f" vel1              : {vel1}")
    print(f" vel2              : {vel2}")
    print(f" z_interface       : {z_interface}")
    print(f" src_info_x        : {src_info_x}")
    print(f" src_info_y        : {src_info_y}")
    print(f" src_info_z        : {src_info_z}")
    print("Parameters parsed")

    # Initialize global parameters from command-line arguments
    # Variables already defined from param.frompar() above
    lx = domain_sizeX
    ly = domain_sizeY
    lz = domain_sizeZ
    nx = ex * order + 1
    ny = ey * order + 1
    nz = ez * order + 1
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
    print(f"boundaries size              : {boundaries_size}")
    print("=========================================")

    # Setup graphic display
    print("Setting up plot...")
    _, _, im = setup_plot(nx, nz, lx=lx, lz=lz)
    print("Plot set up")

    # Initialize Kokkos
    kokkos.initialize()
    print("Kokkos initialized")
    memspace, layout = select_kokkos_memspace(mem)

    # Create bilayer model with custom velocities
    print("Creating bilayer model with velocities...")
    model_obj, elem_to_nodes, coord_node = create_bilayer_model_data_from_params(
        order, (ex, ey, ez), (lx, ly, lz), (vel1, vel2), z_interface, memspace, layout
    )
    print("Model created (bilayer with custom velocities)")
    
    print("Model created")

    print("Finding source element indices from coordinates...")
    print(f"  Source coordinates: x={src_info_x}, y={src_info_y}, z={src_info_z}")
    src_elem_indices = find_source_elements(
        [src_info_x], [src_info_y], [src_info_z],
        elem_to_nodes, coord_node,
        order+1, order+1, order+1
    )


    # Add timing variables
    start_time = time.time()
    simulation_start = datetime.now()
    iteration_times = []
    print(f"Simulation started at: {simulation_start}")

    # Create solver
    print("Creating solver...")
    solver = create_solver(impl, model, order, on_nodes)
    print("Solver created")

    # Initialize model
    print("Initializing model...")
    # compute_fe_init with sponge parameters
    sponge_size = [boundaries_size, boundaries_size, boundaries_size]
    # In C++, surface_sponge is a boolean (true/false).
    # Here we invert it if the argument logic matches "sponge-surface" vs "surface_sponge"
    # Logic in C++: const double distToFrontierX = (surface_sponge_) ? ... : ...
    # Here we pass it directly.
    solver.compute_fe_init(model_obj, sponge_size, surface_sponge, taper_delta)
    # m_syncer->synchronize(m_solver->getMassMatrix(c), par_topology_);
    print("Model initialized")

    # allocate pressure
    print("Allocating Pressure...")
    kk_pnGlobalPrev, pnGlobalPrev, kk_pnGlobalCurr, pnGlobalCurr = allocate_pressure(n_dof, memspace, layout)
    print("Pressure allocated")

    # allocate RHS arrays
    print("Allocating RHS element...")
    kk_RHSElement, RHSElement = allocate_rhs_element(
        n_rhs, ex, ey, ez, memspace, layout, src_elem_indices
    )
    print(f"  - RHS element: shape={RHSElement.shape}, values={RHSElement}")
    print("RHS element allocated")

    print("Allocating RHS weights...")
    # Create source coordinate arrays for all n_rhs sources
    src_x_array = np.full(n_rhs, src_info_x)
    src_y_array = np.full(n_rhs, src_info_y)
    src_z_array = np.full(n_rhs, src_info_z)
    
    computed_weights = compute_gll_weights(
            src_x_array, src_y_array, src_z_array,
            elem_to_nodes,
            coord_node,
            src_elem_indices,
            order+1, order+1, order+1
        )
    print(f"  Using user-provided projection weights (shape: {computed_weights.shape})")
    kk_RHSWeights, rhsWeights = allocate_rhs_weight(n_rhs, model_obj, memspace, layout, computed_weights)
    print(f"RHS weights allocated: shape={rhsWeights.shape}, min={rhsWeights.min():.6f}, max={rhsWeights.max():.6f}")
    print("RHS weights allocated")

    print("Allocating RHS term...")
    kk_RHSTerm, rhsTerm = allocate_rhs_term(
        n_rhs, n_time_steps, dt, f0, memspace, layout
    )
    print(f"RHS term allocated: shape={rhsTerm.shape}")

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
            pnGlobalPrev,
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

    # release solver-related python objects
    del wavefield
    del rhs
    del data
    del solver
    del model_obj

    # release explicit kokkos arrays
    del kk_pnGlobalPrev
    del kk_pnGlobalCurr
    del kk_RHSTerm
    del kk_RHSElement
    del kk_RHSWeights

    # release numpy wrappers/views
    del pnGlobalPrev
    del pnGlobalCurr
    del rhsTerm
    del RHSElement
    del rhsWeights

    # force destructor calls for pybind11-wrapped C++ objects
    import gc
    gc.collect()

    # now finalize Kokkos
    kokkos.finalize()
    print("End of computation")


if __name__ == "__main__":
    main()
