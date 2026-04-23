#!/usr/bin/env python3
"""
Full Waveform Inversion (FWI) Gradient Testing - UNSTRUCTURED ACOUSTIC MODEL

This script performs a single FWI iteration on an acoustic model:
1. Forward propagation #1 (true model): generates observed seismograms
2. Forward propagation #2 (initial model): generates synthetic seismograms and wavefields
3. Compute residuals (observed - synthetic)
4. Adjoint propagation: backpropagates residuals
5. Compute gradient: accumulates sensitivity over adjoint-forward interaction

Domain: 200 x 200 x 200 m cubic
Model:
  - Mesh: UNSTRUCTURED (created via CartesianUnstructBuilder with proper initialization)
  - Density: constant 1000 kg/m³
  - Vp: Constant velocity (can be set per model)

Acquisition:
  - Source: (75, 75, 75)
  - Receiver: (125, 125, 75)

Configuration:
  - MODEL_DISCRETIZATION: Set to "ONNODES" or "ONELEMENTS" to control gradient computation
    * "ONNODES": Model properties defined at nodes (shared between elements, uses ATOMICADD)
    * "ONELEMENTS": Model properties defined on elements (unique per element, no atomic ops)

Output:
  - observed_seismo.npy: Observed seismograms (true model)
  - synthetic_seismo.npy: Synthetic seismograms (initial model)
  - wavefield_snapshots.npy: Forward wavefields at snapshot times
  - gradient_kappa.npy: Inverse property gradient (bulk modulus)
  - gradient_buoyancy.npy: Inverse property gradient (density)
"""

import os
import time
import numpy as np
import gc

import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pyfuntides.gradient as Gradient

# Aliases for unstructured model creation
CartesianUnstructBuilder = Model.CartesianUnstructBuilder_f32_i32
CartesianParams = Model.CartesianParams_f32_i32
ModelUnstruct = Model.ModelUnstruct_f32_i32
ModelUnstructData = Model.ModelUnstructData_f32_i32


# =============================================================================
# Configuration
# =============================================================================

# **MODEL DISCRETIZATION CHOICE** - Set this to control gradient computation
MODEL_DISCRETIZATION = "ONELEMENTS"  # Options: "ONNODES" or "ONELEMENTS"

# Domain parameters
DOMAIN_SIZE = 200.0  # m
LX = LY = LZ = DOMAIN_SIZE

# Element and mesh parameters
ORDER = 1 # Polynomial order  
EX = EY = EZ = 100  # Elements per dimension

# Acquisition
SRC_COORD = np.array([75.013, 75.013, 75.013]) # Slightly offset to avoid exact node for better gradient structure
RCV_COORD = np.array([125.031, 125.031, 75.031]) # Slightly offset to avoid exact node for better gradient structure

# Time stepping
DT = 0.001*0.25  # seconds
N_TIME_STEPS = 1000
SNAPSHOT_INTERVAL = 2  # Save wavefields every N timesteps

# Source parameters
F0 = 10.0  # Hz - Ricker wavelet peak frequency

# Physical parameters - TRUE MODEL (observed data)
RHO_TRUE_LAYER1 = 1000.0  # kg/m³ 
RHO_TRUE_LAYER2 = 2000.0  # kg/m³ 
VP_TRUE_LAYER1 = 1500.0  # m/s - upper layer velocity for true model
VP_TRUE_LAYER2 = 1500.0 #2000.0  # m/s - lower layer velocity for true model
LAYER_BOUNDARY = EZ // 2  # Interface at middle (z = 100 m)

# Physical parameters - INITIAL MODEL (synthetic data)
# Used for inversion testing: should differ from true model to create a non-zero gradient
RHO_INITIAL_LAYER1 = 1400.0  # kg/m³ 
RHO_INITIAL_LAYER2 = 1600.0  # kg/m³ 
VP_INITIAL_LAYER1 = 1500.0 #1620.0  # m/s - upper layer velocity for initial model (perturbed)
VP_INITIAL_LAYER2 = 1500.0 #1920.0  # m/s - lower layer velocity for initial model (perturbed)

# Kokkos settings
MEM_SPACE = kokkos.HostSpace
LAYOUT = kokkos.LayoutRight

# Output directory
OUTPUT_DIR = "./fwi_outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# Utility Functions
# =============================================================================

def reset_solver_state(solver, data, fields):
    """
    Reset solver state between forward propagations.
    
    Explicitly zeros all Kokkos arrays in the fields dictionary to ensure
    a completely clean state for the next forward pass.
    
    Parameters
    ----------
    solver : Solver.Solver
        The solver instance to reset
    data : Solver.SEMsolverDataAcoustic
        The solver data to reset
    fields : dict
        Dictionary containing wavefields and all Kokkos arrays
    """
    print("  Resetting solver state...")
    
    # Zero out all wavefield Kokkos arrays
    for key in ['p_curr_kk', 'p_prev_kk']:
        if key in fields:
            kk_array = fields[key]
            np_array = np.array(kk_array, copy=False)
            np_array[:] = 0.0
            print(f"    {key}: zeroed, range=[{np_array.min():.6e}, {np_array.max():.6e}]")
    
    # Zero out RHS Kokkos arrays
    for key in ['rhs_kk', 'rhs_element_kk', 'rhs_weight_kk']:
        if key in fields:
            kk_array = fields[key]
            np_array = np.array(kk_array, copy=False)
            np_array[:] = 0.0
            print(f"    {key}: zeroed, size={np_array.size}")
    
    # Zero out all Kokkos arrays in model_kokkos_refs dict
    if 'model_kokkos_refs' in fields:
        model_refs = fields['model_kokkos_refs']
        for key, kk_array in model_refs.items():
            np_array = np.array(kk_array, copy=False)
            np_array[:] = 0.0
            print(f"    model_kokkos_refs[{key}]: zeroed, size={np_array.size}")
    
    # Swap wavefields to ensure consistency in leapfrog scheme
    data.swap_wavefields()
    print("    Wavefields swapped for consistency")


def source_term_ricker(time_n, f0):
    """Ricker wavelet source term."""
    t_peak = 1.2 / f0
    if time_n <= -0.9 * t_peak or time_n >= 2.9 * t_peak:
        return 0.0
    
    pi = np.pi
    lam = (f0 * pi) ** 2
    pulse = (
        2.0
        * lam
        * (2.0 * lam * (time_n - t_peak) ** 2 - 1.0)
        * np.exp(-lam * (time_n - t_peak) ** 2)
    )
    return pulse


def coord_to_grid_index(coord, domain_size, n_elements, order):
    """
    Convert physical coordinate to grid index.
    
    Parameters
    ----------
    coord : float
        Physical coordinate
    domain_size : float
        Domain size in this dimension
    n_elements : int
        Number of elements in this dimension
    order : int
        Polynomial order
    
    Returns
    -------
    int
        Grid index (clamped to valid range)
    """
    # Handle both ONNODES and ONELEMENTS discretization
    if MODEL_DISCRETIZATION == "ONNODES":
        n_points = n_elements * order + 1
        h = domain_size / (n_elements * order)  # Grid spacing (2m for our case)
    else:  # ONELEMENTS
        n_points = n_elements
        h = domain_size / n_elements  # Grid spacing (4m for our case)
    
    idx = int(np.round(coord / h))
    return np.clip(idx, 0, n_points - 1)


def coord_to_element_index(coord, domain_size, n_elements, order):
    """
    Convert physical coordinate to element index.
    
    Parameters
    ----------
    coord : ndarray
        (x, y, z) physical coordinates in meters
    domain_size : float
        Domain size (assumed cubic)
    n_elements : int
        Elements per dimension
    order : int
        Polynomial order
    
    Returns
    -------
    int
        1D element index in the flattened array
    """
    h_elem = domain_size / n_elements  # Element size
    
    # Find element index for each coordinate
    i = int(np.clip(coord[0] / h_elem, 0, n_elements - 1))
    j = int(np.clip(coord[1] / h_elem, 0, n_elements - 1))
    k = int(np.clip(coord[2] / h_elem, 0, n_elements - 1))
    
    # Row-major layout: index = i + j*n_elements + k*n_elements*n_elements
    idx_1d = i + j * n_elements + k * n_elements * n_elements
    return idx_1d


def get_receiver_index_3d(coord, domain_size, n_elements, order):
    """
    Get the 1D DOF index for a 3D receiver coordinate.
    
    RECEIVER DATA IS ALWAYS EXTRACTED FROM NODE-BASED WAVEFIELDS,
    regardless of where MODEL parameters are defined.
    
    Parameters
    ----------
    coord : ndarray
        (x, y, z) physical coordinates
    domain_size : float
        Domain size (assumed cubic)
    n_elements : int
        Elements per dimension
    order : int
        Polynomial order
    
    Returns
    -------
    int
        1D node-based DOF index in the wavefield array
    """
    # Wavefields are ALWAYS node-based, even if model is on elements
    # So receiver reading ALWAYS uses node-based grid spacing
    n_points = n_elements * order + 1  # Number of nodes per dimension
    h = domain_size / (n_elements * order)  # Grid spacing for nodes
    
    # Convert physical coordinates to node indices
    i = int(np.round(coord[0] / h))
    i = np.clip(i, 0, n_points - 1)
    
    j = int(np.round(coord[1] / h))
    j = np.clip(j, 0, n_points - 1)
    
    k = int(np.round(coord[2] / h))
    k = np.clip(k, 0, n_points - 1)
    
    # Row-major layout: index = i + j*nx + k*nx*ny
    idx_1d = i + j * n_points + k * n_points * n_points
    return idx_1d


def allocate_field(n_dof, init_value=0.0):
    """Allocate and initialize a Kokkos field."""
    kk_field = kokkos.array([n_dof], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    np_field = np.array(kk_field, copy=False)
    np_field[:] = init_value
    return kk_field, np_field


def create_bilayer_model_data_from_params(order, n_elems, L, vels, rhos, z_elem_interface, on_nodes=True):
    """
    Create bilayer unstructured model data with custom velocities.
    
    Parameters
    ----------
    order : int
        Polynomial order
    n_elems : tuple of int
        (ex, ey, ez) number of elements
    L : tuple of float
        (lx, ly, lz) domain sizes
    vels : tuple of float
        (vel_layer1, vel_layer2) velocities for two layers
    rhos : tuple of float
        (rho_layer1, rho_layer2) densities for two layers
    z_elem_interface : int
        Element index where interface starts in z-direction
    on_nodes : bool
        Whether model properties are defined on nodes (True) or elements (False)
    
    Returns
    -------
    tuple
        (model_obj, elem_to_nodes, coord_node) where:
        - model_obj: ModelUnstruct object (fully initialized with face connectivity)
        - elem_to_nodes: numpy array of element-to-node connectivity
        - coord_node: tuple of (x_coords, y_coords, z_coords)
    """
    lx, ly, lz = L
    ex, ey, ez = n_elems
    n_elem = ex * ey * ez

    # Create mesh with bilayer velocity model
    mesh_model_vp = np.zeros((n_elem,), dtype=np.float32)
    mesh_model_rho = np.ones((n_elem,), dtype=np.float32)

    for k in range(ez):
        for j in range(ey):
            for i in range(ex):
                e = i + j * ex + k * ex * ey
                if k < z_elem_interface:
                    mesh_model_vp[e] = vels[0]
                    mesh_model_rho[e] = rhos[0]
                else:
                    mesh_model_vp[e] = vels[1]
                    mesh_model_rho[e] = rhos[1]

    # Element-to-node connectivity and coordinates
    nx_nodes = order * ex + 1
    ny_nodes = order * ey + 1
    nz_nodes = order * ez + 1
    n_node = nx_nodes * ny_nodes * nz_nodes

    dx_node = lx / (order * ex)
    dy_node = ly / (order * ey)
    dz_node = lz / (order * ez)

    nodes_coords_x = np.zeros(n_node, dtype=np.float32)
    nodes_coords_y = np.zeros(n_node, dtype=np.float32)
    nodes_coords_z = np.zeros(n_node, dtype=np.float32)

    for k in range(nz_nodes):
        for j in range(ny_nodes):
            for i in range(nx_nodes):
                idx = i + j * nx_nodes + k * nx_nodes * ny_nodes
                nodes_coords_x[idx] = i * dx_node
                nodes_coords_y[idx] = j * dy_node
                nodes_coords_z[idx] = k * dz_node

    n_gll_per_elem = (order + 1) ** 3
    elem_to_nodes = np.zeros((n_elem, n_gll_per_elem), dtype=np.int32)

    for ez_i in range(ez):
        for ey_i in range(ey):
            for ex_i in range(ex):
                e = ex_i + ey_i * ex + ez_i * ex * ey
                for k_loc in range(order + 1):
                    for j_loc in range(order + 1):
                        for i_loc in range(order + 1):
                            i_glob = ex_i * order + i_loc
                            j_glob = ey_i * order + j_loc
                            k_glob = ez_i * order + k_loc
                            node_idx = i_glob + j_glob * nx_nodes + k_glob * nx_nodes * ny_nodes
                            local_idx = i_loc + j_loc * (order + 1) + k_loc * (order + 1) ** 2
                            elem_to_nodes[e, local_idx] = node_idx

    # Create Kokkos views
    kk_elem_to_nodes = kokkos.array(elem_to_nodes, dtype=kokkos.int32, space=MEM_SPACE, layout=LAYOUT)
    kk_nodes_coords_x = kokkos.array([n_node], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    kk_nodes_coords_y = kokkos.array([n_node], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    kk_nodes_coords_z = kokkos.array([n_node], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    np.array(kk_nodes_coords_x, copy=False)[:] = nodes_coords_x
    np.array(kk_nodes_coords_y, copy=False)[:] = nodes_coords_y
    np.array(kk_nodes_coords_z, copy=False)[:] = nodes_coords_z

    kk_model_vp_element = kokkos.array([n_elem], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    kk_model_rho_element = kokkos.array([n_elem], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    np.array(kk_model_vp_element, copy=False)[:] = mesh_model_vp
    np.array(kk_model_rho_element, copy=False)[:] = mesh_model_rho
    
    # Create node-based model arrays with correct size
    mesh_model_vp_node = np.zeros((n_node,), dtype=np.float32)
    mesh_model_rho_node = np.zeros((n_node,), dtype=np.float32)
    
    # Fill node-based arrays by checking node z-coordinate
    z_boundary_phys = (z_elem_interface + 0.5) * (lz / ez)  # Physical z coordinate of interface
    for node_idx in range(n_node):
        node_z = nodes_coords_z[node_idx]
        if node_z < z_boundary_phys:
            mesh_model_vp_node[node_idx] = vels[0]
            mesh_model_rho_node[node_idx] = rhos[0]
        else:
            mesh_model_vp_node[node_idx] = vels[1]
            mesh_model_rho_node[node_idx] = rhos[1]
    
    kk_model_vs_element = kokkos.array([1], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    kk_model_vp_node = kokkos.array([n_node], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    kk_model_rho_node = kokkos.array([n_node], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    kk_model_vs_node = kokkos.array([1], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    
    # Fill node-based Kokkos arrays
    np.array(kk_model_vp_node, copy=False)[:] = mesh_model_vp_node
    np.array(kk_model_rho_node, copy=False)[:] = mesh_model_rho_node
    
    # Diagnostic: Print node-based model values
    print(f"\n  [DIAGNOSTIC] Node-based model (Vp):")
    print(f"    z_boundary_phys = {z_boundary_phys} m")
    print(f"    Number of nodes: {n_node}")
    print(f"    Vp values range: [{mesh_model_vp_node.min():.1f}, {mesh_model_vp_node.max():.1f}] m/s")
    print(f"    Layer 1 (vels[0]): {vels[0]} m/s, count: {np.sum(mesh_model_vp_node == vels[0])}")
    print(f"    Layer 2 (vels[1]): {vels[1]} m/s, count: {np.sum(mesh_model_vp_node == vels[1])}")
    print(f"    Rho values range: [{mesh_model_rho_node.min():.1f}, {mesh_model_rho_node.max():.1f}] m/s")
    print(f"    Layer 1 (rhos[0]): {rhos[0]} m/s, count: {np.sum(mesh_model_rho_node == rhos[0])}")
    print(f"    Layer 2 (rhos[1]): {rhos[1]} m/s, count: {np.sum(mesh_model_rho_node == rhos[1])}")
    
    # Sample some node values to inspect distribution
    sample_indices = [0, n_node//4, n_node//2, 3*n_node//4, n_node-1]
    print(f"    Sample nodes (indices {sample_indices}):")
    for idx in sample_indices:
        k_z = idx // (nx_nodes * ny_nodes)
        z_phys = nodes_coords_z[idx]
        print(f"      node[{idx}]: z_idx={k_z}, z_phys={z_phys:.2f}m, vp={mesh_model_vp_node[idx]:.1f} m/s")
        print(f"      node[{idx}]: z_idx={k_z}, z_phys={z_phys:.2f}m, rho={mesh_model_rho_node[idx]:.1f} kg/m³")

    kk_empty_1d_float = kokkos.array([1], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    kk_empty_3d = kokkos.array([1, 1, 1], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    # FIX: Create TRULY empty boundaries_t_ array (extent 0, not 1)
    # This ensures isBoundaryFace() uses the correct face_connectivity_ logic
    # Create truly empty boundaries array (extent(0) == 0)
    # This is critical: must be truly empty so isBoundaryFace() uses correct logic branch
    # Essential for multiple solver runs to avoid state persistence issues
    kk_empty_boundaries_t = kokkos.array([0], dtype=kokkos.int32, space=MEM_SPACE, layout=LAYOUT)

    model_data = ModelUnstructData(
        order, n_elem, n_node, lx, ly, lz, on_nodes, False,
        kk_elem_to_nodes,
        kk_nodes_coords_x, kk_nodes_coords_y, kk_nodes_coords_z,
        kk_model_vp_node, kk_model_vp_element,
        kk_model_rho_node, kk_model_rho_element,
        kk_model_vs_node, kk_model_vs_element,
        kk_empty_1d_float, kk_empty_1d_float, kk_empty_1d_float, kk_empty_1d_float,
        kk_empty_1d_float, kk_empty_1d_float, kk_empty_1d_float, kk_empty_1d_float,
        kk_empty_1d_float, kk_empty_1d_float, kk_empty_3d, kk_empty_boundaries_t
    )
    
    model_obj = ModelUnstruct(model_data)
    
    # Build face connectivity - CRITICAL
    model_obj.build_face_connectivity()
    
    # The model stores references to these arrays internally, so we must keep
    # Python references alive throughout the solver's lifetime
    kokkos_refs = {
        'elem_to_nodes': kk_elem_to_nodes,
        'coords_x': kk_nodes_coords_x,
        'coords_y': kk_nodes_coords_y,
        'coords_z': kk_nodes_coords_z,
        'vp_element': kk_model_vp_element,
        'rho_element': kk_model_rho_element,
        'vs_element': kk_model_vs_element,
        'vp_node': kk_model_vp_node,
        'rho_node': kk_model_rho_node,
        'vs_node': kk_model_vs_node,
    }
    
    return model_obj, elem_to_nodes, (nodes_coords_x, nodes_coords_y, nodes_coords_z), kokkos_refs


def create_layer_unstruct_model(n_elements, order, on_nodes=True,
                                   vp_layer1=None, vp_layer2=None, z_boundary=None,
                                   rho_layer1=None, rho_layer2=None):
    """
    Create an unstructured acoustic model with bilayer velocities.
    
    Parameters
    ----------
    n_elements : int
        Number of elements per dimension
    order : int
        Polynomial order (1, 2, or 3)
    on_nodes : bool
        Whether model properties are on nodes (for API compatibility - not used)
    vp_layer1 : float
        Velocity for layer 1
    vp_layer2 : float
        Velocity for layer 2 (if None, use same as layer 1 for single layer)
    z_boundary : int
        Element z-index for interface (if None, use ez/2)
    rho_layer1 : float
        Density for layer 1
    rho_layer2 : float
        Density for layer 2 (if None, use same as layer 1 for single layer)   
    
    Returns
    -------
    tuple
        (model, elem_to_nodes, coord_node)
    """
    
    if vp_layer1 is None:
        vp_layer1 = VP_TRUE_LAYER1
    
    if vp_layer2 is None:
        vp_layer2 = vp_layer1  # Single layer if layer2 not specified

    if rho_layer1 is None:
        rho_layer1 = RHO_TRUE_LAYER1
    
    if rho_layer2 is None:
        rho_layer2 = rho_layer1  # Single layer if layer2 not specified
    
    if z_boundary is None:
        z_boundary = n_elements // 2
    
    # Create bilayer model with face connectivity properly initialized
    model, elem_to_nodes, coord_node, kokkos_refs = create_bilayer_model_data_from_params(
        order, 
        (n_elements, n_elements, n_elements),
        (DOMAIN_SIZE, DOMAIN_SIZE, DOMAIN_SIZE),
        (vp_layer1, vp_layer2),
        (rho_layer1, rho_layer2),
        z_boundary,
        on_nodes=on_nodes
    )
    
    print(f"  Created UNSTRUCTURED bilayer model:")
    print(f"    Mesh: {n_elements}^3 elements, order {order}")
    print(f"    Layer 1 velocity: {vp_layer1:.1f} m/s (z < z_interface={z_boundary})")
    print(f"    Layer 2 velocity: {vp_layer2:.1f} m/s (z >= z_interface={z_boundary})")
    print(f"    Layer 1 density: {rho_layer1:.1f} kg/m³ (z < z_interface={z_boundary})")
    print(f"    Layer 2 density: {rho_layer2:.1f} kg/m³ (z >= z_interface={z_boundary})")
    
    return model, elem_to_nodes, coord_node, kokkos_refs


def find_source_element(src_coord, elem_to_nodes, coord_node, order):
        """Find element containing source coordinate."""
        # RHS location (source position) - find element containing source, using actual mesh coordinates
        n_elem = elem_to_nodes.shape[0]
        n_gll = order + 1
        
        for elem_idx in range(n_elem):
            node_indices = elem_to_nodes[elem_idx, :]
            
            # Get bounding box from corners
            corners = [
                node_indices[0],  # (0,0,0)
                node_indices[n_gll-1],  # (n_gll-1, 0, 0)
                node_indices[n_gll * (n_gll-1)],  # (0, n_gll-1, 0)
                node_indices[n_gll * n_gll - 1],  # (n_gll-1, n_gll-1, 0)
                node_indices[n_gll * n_gll * (n_gll-1)],  # (0, 0, n_gll-1)
                node_indices[-1],  # Last corner
            ]
            
            x_corners = [coord_node[0][c] for c in corners]
            y_corners = [coord_node[1][c] for c in corners]
            z_corners = [coord_node[2][c] for c in corners]
            
            x_min, x_max = min(x_corners), max(x_corners)
            y_min, y_max = min(y_corners), max(y_corners)
            z_min, z_max = min(z_corners), max(z_corners)
            
            tol = 1e-3
            if (x_min - tol <= src_coord[0] <= x_max + tol and
                y_min - tol <= src_coord[1] <= y_max + tol and
                z_min - tol <= src_coord[2] <= z_max + tol):
                return elem_idx
        
        return 0  # Default to first element


def setup_forward_solver(model_type='initial'):
    """
    Set up unstructured solver, layer model, and allocate fields for forward propagation.
    
    Parameters
    ----------
    model_type : str
        'true' for true model (observed data)
        'initial' for initial model (synthetic data) - perturbed version
    """
    print(f"Setting up forward solver ({model_type} model)...")
    
    # Determine model parameters based on model_type
    on_nodes = (MODEL_DISCRETIZATION == "ONNODES")
    
    if model_type == 'true':
        vp_layer1 = VP_TRUE_LAYER1
        vp_layer2 = VP_TRUE_LAYER2
        print(f"  Velocity: Layer 1 = {vp_layer1} m/s, Layer 2 = {vp_layer2} m/s")
        rho_layer1 = RHO_TRUE_LAYER1
        rho_layer2 = RHO_TRUE_LAYER2
        print(f"  Density: Layer 1 = {rho_layer1} kg/m³, Layer 2 = {rho_layer2} kg/m³")
    elif model_type == 'initial':
        vp_layer1 = VP_INITIAL_LAYER1
        vp_layer2 = VP_INITIAL_LAYER2
        print(f"  Velocity: Layer 1 = {vp_layer1} m/s, Layer 2 = {vp_layer2} m/s")
        rho_layer1 = RHO_INITIAL_LAYER1
        rho_layer2 = RHO_INITIAL_LAYER2
        print(f"  Density: Layer 1 = {rho_layer1} kg/m³, Layer 2 = {rho_layer2} kg/m³")
        print(f"  Perturbation L1 (vp): {100*(VP_INITIAL_LAYER1-VP_TRUE_LAYER1)/VP_TRUE_LAYER1:.1f}%")
        print(f"  Perturbation L2 (vp): {100*(VP_INITIAL_LAYER2-VP_TRUE_LAYER2)/VP_TRUE_LAYER2:.1f}%")
        print(f"  Perturbation L1 (rho): {100*(RHO_INITIAL_LAYER1-RHO_TRUE_LAYER1)/RHO_TRUE_LAYER1:.1f}%")
        print(f"  Perturbation L2 (rho): {100*(RHO_INITIAL_LAYER2-RHO_TRUE_LAYER2)/RHO_TRUE_LAYER2:.1f}%")
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    # Create unstructured bilayer model with specified layer velocities
    model, elem_to_nodes, coord_node, kokkos_refs = create_layer_unstruct_model(
        EX, ORDER, on_nodes=on_nodes,
        vp_layer1=vp_layer1, vp_layer2=vp_layer2, 
        z_boundary=LAYER_BOUNDARY, 
        rho_layer1=rho_layer1, rho_layer2=rho_layer2
    )
    discretization_str = "on NODES" if on_nodes else "on ELEMENTS"
    print(f"  Model created: {EX}^3 elements, order {ORDER}, model {discretization_str}")
        
    # Create UNSTRUCTURED solver with chosen discretization
    model_location = (Solver.ModelLocationType.ONNODES if on_nodes 
                      else Solver.ModelLocationType.ONELEMENTS)
    solver = Solver.create_solver(
        Solver.MethodType.SEM,
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.UNSTRUCT,  # Unstructured mesh to match model
        model_location,
        Solver.PhysicType.ACOUSTIC,
        ORDER
    )
    print(f"  Solver created: MAKUTU implementation (UNSTRUCTURED)")
    
    # Initialize solver
    sponge_size = [0.0, 0.0, 0.0]  # No absorbing boundaries for simplicity
    solver.compute_fe_init(model, sponge_size, False, 0.015)
    print(f"  Solver initialized")
    
    # Diagnostic: Check that model Vp was properly set
    print(f"  Model verification (kokkos_refs):")
    kk_vp_elem = kokkos_refs['vp_element']
    np_vp_elem = np.array(kk_vp_elem, copy=False)
    print(f"    Element Vp range: [{np_vp_elem.min():.1f}, {np_vp_elem.max():.1f}] m/s")
    print(f"    Layer 1 count: {np.sum(np_vp_elem == vp_layer1)}")
    print(f"    Layer 2 count: {np.sum(np_vp_elem == vp_layer2)}")
    
    # Allocate fields (wavefields are ALWAYS on nodes, regardless of model discretization)
    # Model discretization only affects where model parameters live, not wavefield resolution
    n_dof = ((EX * ORDER + 1) ** 3)  # Always node-based
    kk_p_prev, np_p_prev = allocate_field(n_dof, 0.0)
    kk_p_curr, np_p_curr = allocate_field(n_dof, 0.0)
    
    # Create RHS arrays (1 source)
    kk_rhs_term = kokkos.array(
        [1, N_TIME_STEPS], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT
    )
    np_rhs_term = np.array(kk_rhs_term, copy=False)
    for t in range(N_TIME_STEPS):
        np_rhs_term[0, t] = source_term_ricker(t * DT, F0)
    
    src_elem_idx = find_source_element(SRC_COORD, elem_to_nodes, coord_node, ORDER)
    kk_rhs_element = kokkos.array([1], dtype=kokkos.int32, space=MEM_SPACE, layout=LAYOUT)
    np_rhs_element = np.array(kk_rhs_element, copy=False)
    np_rhs_element[0] = src_elem_idx
    
    # Diagnostic: Check source element and model at source
    print(f"  Source diagnostics:")
    print(f"    Source element index: {src_elem_idx}")
    print(f"    Source Vp (element): {kokkos_refs['vp_element'][src_elem_idx]:.1f} m/s")
    print(f"    RHS source term range: [{np_rhs_term.min():.6e}, {np_rhs_term.max():.6e}]")
    print(f"    RHS source term at t=0,50,100: [{np_rhs_term[0,0]:.6e}, {np_rhs_term[0,50]:.6e}, {np_rhs_term[0,100]:.6e}]")
    
    # RHS weights
    n_pts_per_elem = (ORDER + 1) ** 3
    kk_rhs_weight = kokkos.array(
        [1, n_pts_per_elem], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT
    )
    np_rhs_weight = np.array(kk_rhs_weight, copy=False)
    np_rhs_weight[0, :] = 1.0 / n_pts_per_elem
    
    # Create wavefield and RHS
    wavefield = Solver.WavefieldAcoustic(kk_p_prev, kk_p_curr)
    rhs = Solver.RhsAcoustic(kk_rhs_term, kk_rhs_element, kk_rhs_weight)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)
    
    return solver, model, data, {
        'n_dof': n_dof,
        'p_prev_kk': kk_p_prev,
        'p_curr_kk': kk_p_curr,
        'wavefield': wavefield,
        'elem_to_nodes': elem_to_nodes,
        'coord_node': coord_node,
        'rhs_kk': kk_rhs_term,  # Keep RHS alive
        'rhs_element_kk': kk_rhs_element,  # Keep RHS element alive
        'rhs_weight_kk': kk_rhs_weight,  # Keep RHS weights alive
        'model_kokkos_refs': kokkos_refs,  # Keep model's Kokkos arrays alive
    }


def forward_propagation(solver, data, fields, output_key, save_snapshots=False):
    """
    Perform forward propagation and extract receiver data.
    
    Parameters
    ----------
    solver : Solver.Solver
        FUnTiDES solver
    data : Solver.SEMsolverDataAcoustic
        Solver data
    fields : dict
        Dictionary containing Kokkos pressure fields and metadata
    output_key : str
        Label for console output
    save_snapshots : bool
        If True, save wavefields snapshots
    
    Returns
    -------
    seismo : ndarray
        Receiver seismogram [n_time_steps]
    snapshots : ndarray or None
        Wavefields at snapshot times [n_snapshots, n_dof] if saved, else None
    """
    print(f"\n{output_key}:")
    print(f"  Running forward propagation for {N_TIME_STEPS} timesteps...")
    
    n_dof = fields['n_dof']
    kk_p_curr = fields['p_curr_kk']
    kk_p_prev = fields['p_prev_kk']
    wavefield = fields['wavefield']
    rcv_idx = get_receiver_index_3d(RCV_COORD, DOMAIN_SIZE, EX, ORDER)
    
    print(f"  Receiver index: {rcv_idx}")
    
    # Diagnostic: Verify wavefields are initialized to zero
    np_p_curr_check = np.array(kk_p_curr, copy=False)
    np_p_prev_check = np.array(kk_p_prev, copy=False)
    print(f"  Initial wavefields (should be zero):")
    print(f"    p_curr: [{np_p_curr_check.min():.6e}, {np_p_curr_check.max():.6e}]")
    print(f"    p_prev: [{np_p_prev_check.min():.6e}, {np_p_prev_check.max():.6e}]")
    del np_p_curr_check, np_p_prev_check
    
    seismo = np.zeros(N_TIME_STEPS, dtype=np.float32)
    snapshots_list = []
    
    start_time = time.time()
    
    for t in range(N_TIME_STEPS):
        # Compute step
        solver.compute_forces(DT, t, data)
        solver.update_solution(DT, data)
        
        # Get CURRENT field reference FRESH from wavefield (not stored once)
        kk_p_curr = wavefield.get_current_field(0)
        np_p_curr = np.array(kk_p_curr, copy=False)
        
        # Extract receiver data from current pressure field
        seismo[t] = np_p_curr[rcv_idx]
        
        # Diagnostic: Print wavefield range at key timesteps
        if t < 5 or t in [10, 50, 100, 200, 500, 999] or (t + 1) % 500 == 0:
            wf_min = np_p_curr.min()
            wf_max = np_p_curr.max()
            wf_rms = np.sqrt(np.mean(np_p_curr**2))
            rcv_val = np_p_curr[rcv_idx]
            print(f"    Step {t:4d}: wavefield=[{wf_min:10.6e}, {wf_max:10.6e}], RMS={wf_rms:10.6e}, receiver={rcv_val:10.6e}")
        
        # Save snapshots if requested
        if save_snapshots and t % SNAPSHOT_INTERVAL == 0:
            snapshots_list.append(np_p_curr.copy())
        
        # Swap wavefields for next iteration
        data.swap_wavefields()
        
        # Delete the local numpy view
        del np_p_curr
    
    elapsed = time.time() - start_time
    print(f"  Completed in {elapsed:.2f} seconds")
    print(f"  Receiver index: {rcv_idx}")
    print(f"  Seismogram range: [{seismo.min():.6e}, {seismo.max():.6e}]")
    
    snapshots = np.array(snapshots_list, dtype=np.float32) if snapshots_list else None
    
    return seismo, snapshots


def setup_adjoint_solver(residual):
    """
    Set up unstructured solver, layer model, and allocate fields for adjoint propagation.
    
    Parameters
    ----------
    residual : ndarray
        Residual seismogram to backpropagate [n_time_steps]
    
    Returns
    -------
    solver : Solver.Solver
        FUnTiDES solver for adjoint
    data : Solver.SEMsolverDataAcoustic
        Adjoint solver data
    fields : dict
        Adjoint wavefield and related arrays
    """
    print("\nSetting up adjoint solver...")
    
    # Create unstructured model (same as forward INITIAL model, with chosen discretization)
    # The adjoint solver needs the INITIAL model since we're computing gradients for the inversion
    on_nodes = (MODEL_DISCRETIZATION == "ONNODES")
    model, elem_to_nodes, coord_node, kokkos_refs = create_layer_unstruct_model(
        EX, ORDER, on_nodes=on_nodes,
        vp_layer1=VP_INITIAL_LAYER1, vp_layer2=VP_INITIAL_LAYER2,
        z_boundary=LAYER_BOUNDARY, 
        rho_layer1=RHO_INITIAL_LAYER1, rho_layer2=RHO_INITIAL_LAYER2
    )
        
    # Create UNSTRUCTURED solver with chosen discretization
    model_location = (Solver.ModelLocationType.ONNODES if on_nodes 
                      else Solver.ModelLocationType.ONELEMENTS)
    solver = Solver.create_solver(
        Solver.MethodType.SEM,
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.UNSTRUCT,  # Unstructured to match model
        model_location,
        Solver.PhysicType.ACOUSTIC,
        ORDER
    )
    print(f"  Solver created: MAKUTU implementation (UNSTRUCTURED)")
    
    # Initialize solver
    sponge_size = [0.0, 0.0, 0.0]
    print(f"  Calling compute_fe_init for adjoint solver...")
    solver.compute_fe_init(model, sponge_size, False, 0.015)
    print(f"  Adjoint solver initialized successfully")
    
    # Allocate adjoint wavefield (wavefields are ALWAYS on nodes, regardless of model discretization)
    # Model discretization only affects where model parameters live, not wavefield resolution
    n_dof = ((EX * ORDER + 1) ** 3)  # Always node-based: 101^3 = 1,030,301
    kk_q_prev, _ = allocate_field(n_dof, 0.0)
    kk_q_curr, _ = allocate_field(n_dof, 0.0)
    
    # Create adjoint RHS (residual as source)
    # Residual is applied at receiver during backward propagation
    kk_rhs_term_adj = kokkos.array(
        [1, N_TIME_STEPS], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT
    )
    np_rhs_term_adj = np.array(kk_rhs_term_adj, copy=False)

    # Compute adjointe receiver index. Corresponds to source index in forward solver.
    rcv_idx = get_receiver_index_3d(SRC_COORD, DOMAIN_SIZE, EX, ORDER)
    print(f"  Adjoint receiver index (source in forward): {rcv_idx}")
    
    # Fill RHS with residual in FORWARD time order
    # The adjoint solver's backward time loop automatically handles the time reversal
    # At each adjoint_step = k, the solver is at physical time (N_TIME_STEPS-1-k),
    # so RHS[k] should contain the residual at that physical time
    for t in range(N_TIME_STEPS):
        # RHS at step t should have residual at physical time (N_TIME_STEPS-1-t)
        np_rhs_term_adj[0, t] = residual[N_TIME_STEPS - 1 - t]

    # RHS location (receiver position - acts as adjoint source)
    # Use grid-based element indexing (same as forward solver source)
    rcv_elem_idx = coord_to_element_index(RCV_COORD, DOMAIN_SIZE, EX, ORDER)
    
    # Validate element index
    if rcv_elem_idx < 0 or rcv_elem_idx >= EX**3:
        raise ValueError(f"Receiver element index {rcv_elem_idx} is out of bounds [0, {EX**3-1}]")
    
    kk_rhs_element_adj = kokkos.array([1], dtype=kokkos.int32, space=MEM_SPACE, layout=LAYOUT)
    np_rhs_element_adj = np.array(kk_rhs_element_adj, copy=False)
    np_rhs_element_adj[0] = rcv_elem_idx
    
    # RHS weights
    n_pts_per_elem = (ORDER + 1) ** 3
    kk_rhs_weight_adj = kokkos.array(
        [1, n_pts_per_elem], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT
    )
    np_rhs_weight_adj = np.array(kk_rhs_weight_adj, copy=False)
    np_rhs_weight_adj[0, :] = 1.0 / n_pts_per_elem
    
    # Validate RHS weights
    if np.any(np.isnan(np_rhs_weight_adj)) or np.any(np.isinf(np_rhs_weight_adj)):
        raise ValueError("RHS weights contain invalid values")
    
    # Create wavefield and RHS
    wavefield_adj = Solver.WavefieldAcoustic(kk_q_prev, kk_q_curr)
    rhs_adj = Solver.RhsAcoustic(kk_rhs_term_adj, kk_rhs_element_adj, kk_rhs_weight_adj)
    data_adj = Solver.SEMsolverDataAcoustic(wavefield_adj, rhs_adj)

    return solver, model, data_adj, {
        'n_dof': n_dof,
        'q_prev_kk': kk_q_prev,
        'q_curr_kk': kk_q_curr,
        'q_prevprev_kk': None,  # Will be allocated in adjoint_propagation
        'wavefield': wavefield_adj,
        'elem_to_nodes': elem_to_nodes,
        'coord_node': coord_node,
        'rcv_idx': rcv_idx, # Source index in forward, receiver index in adjoint
        'rhs_kk': kk_rhs_term_adj,  # Keep RHS alive
        'rhs_element_kk': kk_rhs_element_adj,  # Keep RHS element alive
        'rhs_weight_kk': kk_rhs_weight_adj,  # Keep RHS weights alive
        'model_kokkos_refs': kokkos_refs,  # Keep model's Kokkos arrays alive
    }


def adjoint_propagation(solver_adj, model_adj, data_adj, fields_adj, 
                       forward_snapshots, output_key):
    """
    Perform adjoint propagation and compute gradient.
    
    Uses forward snapshots to compute gradient contributions at each adjoint step.
    Implements the standard adjoint-state method for FWI.
    
    Parameters
    ----------
    solver_adj : Solver.Solver
        FUnTiDES solver for adjoint
    model_adj : Model.ModelStruct
        Acoustic model
    data_adj : Solver.SEMsolverDataAcoustic
        Solver data for adjoint with residual as source
    fields_adj : dict
        Adjoint wavefield arrays and metadata
    forward_snapshots : tuple
        (snapshot_indices, forward_wavefield_array) where:
        - snapshot_indices: list of time step indices corresponding to snapshots
        - forward_wavefield_array: [n_snapshots, n_dof] array
    output_key : str
        Label for console output
    
    Returns
    -------
    gradient : dict
        Gradients: {'kappa': ndarray, 'buoyancy': ndarray}
    """
    print(f"\n{output_key}:")
    print(f"  Running adjoint propagation for {N_TIME_STEPS} timesteps...")
    
    n_dof = fields_adj['n_dof']
    
    # If these are collected, the mesh connectivity may become corrupted
    model_kokkos_refs = fields_adj.get('model_kokkos_refs', {})
    
    # VALIDATION: Verify mesh state before entering adjoint loop
    print(f"  Pre-adjoint-loop mesh validation:")
    print(f"    Wavefield DOF: {n_dof}")
    print(f"    Expected max node index in connectivity: {np.max(fields_adj.get('elem_to_nodes', np.array([0])))}")
    print(f"    Model Kokkos refs dict size: {len(model_kokkos_refs)}")
    
    # Create differentiator for gradient computation (must match solver discretization)
    model_loc_type = (Solver.ModelLocationType.ONNODES 
                      if MODEL_DISCRETIZATION == "ONNODES" 
                      else Solver.ModelLocationType.ONELEMENTS)
    differentiator = Gradient.create_differentiator(
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.UNSTRUCT,  # Unstructured to match forward/adjoint solvers
        model_loc_type,
        Solver.PhysicType.ACOUSTIC,
        ORDER
    )
    
    # Allocate gradient arrays (must match model discretization, not wavefield resolution)
    if MODEL_DISCRETIZATION == "ONNODES":
        grad_dof = n_dof  # Node-based gradient
    else:  # ONELEMENTS
        grad_dof = EX * EY * EZ  # Element-based gradient
    
    kk_grad_kappa, np_grad_kappa = allocate_field(grad_dof, 0.0)
    kk_grad_buoyancy, np_grad_buoyancy = allocate_field(grad_dof, 0.0)
    grad_acoustic = Gradient.GradientAcoustic(kk_grad_kappa, kk_grad_buoyancy)
    
    # Unpack forward snapshots
    snapshot_indices, forward_wavefields = forward_snapshots
    n_snapshots = len(snapshot_indices)
    
    print(f"  Forward snapshots: {n_snapshots} at indices {snapshot_indices[:5]}...") 
    print(f"  Total adjoint time steps: {N_TIME_STEPS}")
    print(f"  Expected gradient computation steps: ~{n_snapshots}")
    
    # Allocate buffers for adjoint wavefield prevprev (for gradient computation)
    kk_q_prevprev, _ = allocate_field(n_dof, 0.0)

    # Pre-allocate forward wavefield array for gradient computation (avoids repeated allocation)
    fwd_kk = kokkos.array([n_dof], dtype=kokkos.float32, space=MEM_SPACE, layout=LAYOUT)
    
    # Allocate adjoint seismogram to diagnose swap behavior
    adjoint_seismo = np.zeros(N_TIME_STEPS, dtype=np.float32)
    rcv_idx = fields_adj['rcv_idx']
    
    start_time = time.time()
    
    # Counter for gradient computation steps
    grad_compute_count = 0
    
    print(f"  Beginning adjoint time loop...")
    
    # Physically starts at time T and propagates back to time 0
    for adjoint_step in range(N_TIME_STEPS):
        # Physical time mapping: adjoint_step=0 corresponds to physical time T
        # adjoint_step=1 corresponds to physical time T-dt, etc.
        physical_time_index = N_TIME_STEPS - 1 - adjoint_step
        
        if adjoint_step % 100 == 0:
            print(f"    Adjoint step {adjoint_step}")
        
        # 1a. BEFORE solver step: save current and previous for use as prevprev after swap
        wavefield_adj = fields_adj['wavefield']
        q_curr_before = wavefield_adj.get_current_field(0)
        q_prev_before = wavefield_adj.get_previous_field(0)
        np_q_curr_before = np.array(q_curr_before, copy=False).copy() #TODO: will be changed once problem with swapping is fixed
        np_q_prev_before = np.array(q_prev_before, copy=False).copy() #TODO: will be changed once problem with swapping is fixed
        
        # 1b. Perform adjoint update step
        solver_adj.compute_forces(DT, adjoint_step, data_adj)
        solver_adj.update_solution(DT, data_adj)
        
        # 1c. Swap wavefields for next iteration (this rotates the leapfrog state)
        data_adj.swap_wavefields()

        # 2. Get adjoint wavefield at current time (AFTER swap, these are the updated fields)
        wavefield_adj = fields_adj['wavefield']
        q_curr_kk = wavefield_adj.get_current_field(0)
        q_prev_kk = wavefield_adj.get_previous_field(0)
        np_q_curr = np.array(q_curr_kk, copy=False)
        np_q_prev = np.array(q_prev_kk, copy=False)
        
        # Extract adjoint seismogram at receiver location
        adjoint_seismo[adjoint_step] = np_q_curr[rcv_idx]
        
        # 3. Find forward snapshot for gradient computation
        # Match forward snapshot at corresponding physical time (backward propagation)
        forward_time = physical_time_index
        
        # Check if forward snapshot exists at this time and snapshot is valid
        should_compute_gradient = False
        snap_idx = -1
        p_fwd = None
        
        if forward_time in snapshot_indices:
            snap_idx = snapshot_indices.index(forward_time)
            if not (snap_idx >= n_snapshots or forward_wavefields is None or len(forward_wavefields) == 0):
                should_compute_gradient = True
                p_fwd = forward_wavefields[snap_idx]
        
        # Only compute gradient if conditions are met
        if should_compute_gradient:
            # 4. Compute gradient via differentiator
            # For gradient computation, use:
            # - Forward wavefield: p_fwd (from snapshot)
            # - Adjoint wavefield: q_curr (current), q_prev (previous), q_prevprev (from before swap)
            
            # NO MANUAL ROTATION NEEDED - the swap_wavefields() already handles leapfrog
            # The saved "before" values become the prevprev after the swap
            
            # Copy forward snapshot into pre-allocated Kokkos array (create local view)
            np_fwd = np.array(fwd_kk, copy=False)
            np_fwd[:] = p_fwd
            del np_fwd
            
            # Copy saved prevprev data (from before swap) into kk_q_prevprev
            np.copyto(np.array(kk_q_prevprev, copy=False), np_q_curr_before) #TODO: will be changed once problem with swapping is fixed
            # Another way to do this while avoiding copies is: np.copyto(np.array(kk_q_prevprev, copy=False), q_prev_before)

            fwd_view = Gradient.WavefieldViewForwardAcoustic(fwd_kk)
            
            # Build backward view from adjoint wavefield
            # After swap: current is at time t, previous is at time t-dt
            # Prevprev (from before swap) is at time t-2dt
            bwd_view = Gradient.WavefieldViewBackwardAcoustic(
                q_curr_kk, q_prev_kk, kk_q_prevprev
            )
            
            # Assemble gradient data
            grad_data = Gradient.GradientDataAcoustic(fwd_view, bwd_view, grad_acoustic)
            
            # 5. Compute and accumulate gradient
            differentiator.compute(model_adj, grad_data, DT)
            grad_compute_count += 1
            
            if grad_compute_count == 1 or adjoint_step % 500 == 0:
                # Check gradient values after computation
                np_grad_k_check = np.array(kk_grad_kappa, copy=False)
                print(f"    Adjoint step {adjoint_step} / {N_TIME_STEPS} (phys_time={physical_time_index}, snap_idx={snap_idx})")
                print(f"      Gradient kappa: [{np_grad_k_check.min():.6e}, {np_grad_k_check.max():.6e}]")
        else:
            # Skip gradient computation - no snapshot or invalid snapshot access
            if adjoint_step % 500 == 0:
                if forward_time not in snapshot_indices:
                    print(f"    Adjoint step {adjoint_step} / {N_TIME_STEPS} (phys time {physical_time_index}, no forward snapshot)")
                else:
                    print(f"    Adjoint step {adjoint_step} / {N_TIME_STEPS} (invalid snapshot access)")
        
        # ALWAYS DELETE local numpy views at end of iteration (even if gradient was skipped)
        del np_q_curr, np_q_prev, np_q_curr_before, np_q_prev_before
    
    elapsed = time.time() - start_time
    print(f"  Completed in {elapsed:.2f} seconds")
    print(f"  Gradient computations performed: {grad_compute_count}")
    
    # Create local numpy views for final output (these will be deleted after return)
    np_grad_kappa_out = np.array(kk_grad_kappa, copy=True)
    np_grad_buoyancy_out = np.array(kk_grad_buoyancy, copy=True)
    
    print(f"  Gradient (kappa) range: [{np_grad_kappa_out.min():.6e}, {np_grad_kappa_out.max():.6e}]")
    print(f"  Gradient (buoyancy) range: [{np_grad_buoyancy_out.min():.6e}, {np_grad_buoyancy_out.max():.6e}]")
    
    # Print adjoint seismogram diagnostics
    print(f"  Adjoint seismogram range: [{adjoint_seismo.min():.6e}, {adjoint_seismo.max():.6e}]")
    
    # Delete temporary Kokkos arrays before function returns
    del fwd_kk, kk_q_prevprev
    del np_grad_kappa, np_grad_buoyancy, grad_acoustic, differentiator
    del kk_grad_kappa, kk_grad_buoyancy  # Also delete gradient Kokkos arrays
    
    return {
        'kappa': np_grad_kappa_out,
        'buoyancy': np_grad_buoyancy_out,
        'adjoint_seismo': adjoint_seismo,
    }


# =============================================================================
# Main FWI Loop
# =============================================================================

def main():
    """Execute single FWI iteration."""
    
    print("\n" + "=" * 70)
    print("FULL WAVEFORM INVERSION - GRADIENT TESTING")
    print("=" * 70)
    
    print("\nConfiguration:")
    print(f"  Domain: {DOMAIN_SIZE}^3 m")
    print(f"  Elements: {EX}^3 (order {ORDER})")
    print(f"  Model discretization: {MODEL_DISCRETIZATION}")
    print(f"  Mesh: UNSTRUCTURED (CartesianUnstructBuilder)")
    print(f"  Wavefield DOF: {(EX * ORDER + 1) ** 3} (always nodes, {EX * ORDER + 1}^3)")
    # Gradient discretization depends on MODEL_DISCRETIZATION
    if MODEL_DISCRETIZATION == "ONNODES":
        print(f"  Gradient DOF: {(EX * ORDER + 1) ** 3} (nodes)")
    else:
        print(f"  Gradient DOF: {EX ** 3} (elements)")
    print(f"  Time steps: {N_TIME_STEPS} x {DT} s = {N_TIME_STEPS * DT} s")
    print(f"  True model Vp (layer 1): {VP_TRUE_LAYER1} m/s, Vp (layer 2): {VP_TRUE_LAYER2} m/s")
    print(f"  Initial model Vp (layer 1): {VP_INITIAL_LAYER1} m/s, Vp (layer 2): {VP_INITIAL_LAYER2} m/s")
    print(f"  True model Rho (layer 1): {RHO_TRUE_LAYER1} kg/m³, Rho (layer 2): {RHO_TRUE_LAYER2} kg/m³")
    print(f"  Initial model Rho (layer 1): {RHO_INITIAL_LAYER1} kg/m³, Rho (layer 2): {RHO_INITIAL_LAYER2} kg/m³")
    print(f"  Source: {SRC_COORD}")
    print(f"  Receiver: {RCV_COORD}")
    print(f"  Output: {OUTPUT_DIR}")
    
    # Initialize Kokkos
    kokkos.initialize()
    print("\nKokkos initialized")
       
    # Initialize variables to None for proper cleanup
    solver_fwd1 = model_true = data_fwd1 = fields_fwd1 = None
    solver_fwd2 = model_init = data_fwd2 = fields_fwd2 = None
    solver_adj = model_adj = data_adj = fields_adj = None
    obs_file = syn_file = snap_file = res_file = grad_k_file = grad_b_file = None
    
    try:
        # =====================================================================
        # STEP 1: Forward propagation with TRUE model (observed data)
        # =====================================================================
        print("\n" + "-" * 70)
        print("STEP 1: FORWARD PROPAGATION (TRUE MODEL)")
        print("-" * 70)
        
        solver_fwd1, model_true, data_fwd1, fields_fwd1 = setup_forward_solver(model_type='true')
        
        obs_seismo, _ = forward_propagation(
            solver_fwd1, data_fwd1, fields_fwd1, "Forward (True model)", save_snapshots=False
        )
        
        # Save observed seismogram
        obs_file = os.path.join(OUTPUT_DIR, "observed_seismo.npy")
        np.save(obs_file, obs_seismo)
        print(f"\n  Observed seismogram saved: {obs_file}")
        
        # CRITICAL: Clean up first forward resources before creating second solver
        # NOTE: Delete objects but do NOT reset_solver_state - let garbage collection clean up naturally
        del solver_fwd1, model_true, data_fwd1, fields_fwd1
        gc.collect()

        # =====================================================================
        # STEP 2: Forward propagation with INITIAL model (synthetic data)
        # =====================================================================
        print("\n" + "-" * 70)
        print("STEP 2: FORWARD PROPAGATION (INITIAL MODEL)")
        print("-" * 70)
        
        solver_fwd2, model_init, data_fwd2, fields_fwd2 = setup_forward_solver(model_type='initial')
        
        syn_seismo, wavefield_snaps = forward_propagation(
            solver_fwd2, data_fwd2, fields_fwd2, "Forward (Initial model)", save_snapshots=True
        )

        # Save synthetic seismogram
        syn_file = os.path.join(OUTPUT_DIR, "synthetic_seismo.npy")
        np.save(syn_file, syn_seismo)
        print(f"\n  Synthetic seismogram saved: {syn_file}")
        
        # Save wavefields
        snap_file = None
        snapshot_indices = []
        if wavefield_snaps is not None:
            snap_file = os.path.join(OUTPUT_DIR, "wavefield_snapshots.npy")
            np.save(snap_file, wavefield_snaps)
            print(f"  Wavefield snapshots ({len(wavefield_snaps)} frames) saved: {snap_file}")
            
            # Record snapshot times (every SNAPSHOT_INTERVAL timesteps)
            snapshot_indices = [t for t in range(0, N_TIME_STEPS, SNAPSHOT_INTERVAL)]
        
        # =====================================================================
        # STEP 3: Compute residuals
        # =====================================================================
        print("\n" + "-" * 70)
        print("STEP 3: COMPUTE RESIDUALS")
        print("-" * 70)
        
        residual = obs_seismo - syn_seismo
        print(f"  Residual range: [{residual.min():.6e}, {residual.max():.6e}]")
        print(f"  Misfit (L2): {np.sum(residual**2):.6e}")
        
        res_file = os.path.join(OUTPUT_DIR, "residual.npy")
        np.save(res_file, residual)
        print(f"  Residuals saved: {res_file}")
        
        # CRITICAL: Clean up forward resources before creating adjoint
        del solver_fwd2, model_init, data_fwd2, fields_fwd2
        gc.collect()
        
        # =====================================================================
        # STEP 4: Adjoint propagation with RESIDUAL as source
        # =====================================================================
        print("\n" + "-" * 70)
        print("STEP 4: ADJOINT PROPAGATION AND GRADIENT")
        print("-" * 70)
        
        # Set up adjoint solver with residual as source
        solver_adj, model_adj, data_adj, fields_adj = setup_adjoint_solver(residual)
        
        # Validate that we have forward snapshots before running adjoint
        if wavefield_snaps is None or len(wavefield_snaps) == 0:
            print("\n  WARNING: No forward snapshots available for gradient computation!")
            print("  Skipping adjoint propagation...")
            gradient = {
                'kappa': np.zeros(fields_adj['n_dof'], dtype=np.float32),
                'buoyancy': np.zeros(fields_adj['n_dof'], dtype=np.float32),
            }
        else:
            # Perform adjoint propagation with gradient computation
            gradient = adjoint_propagation(
                solver_adj, model_adj, data_adj, fields_adj,
                (snapshot_indices, wavefield_snaps),
                "Adjoint propagation"
            )
        
        # Extract gradient arrays
        grad_kappa = gradient['kappa']
        grad_buoyancy = gradient['buoyancy']
        adjoint_seismo = gradient.get('adjoint_seismo', None)
        
        # Save gradients (numpy arrays are copied when saved to disk)
        grad_k_file = os.path.join(OUTPUT_DIR, "gradient_kappa.npy")
        grad_b_file = os.path.join(OUTPUT_DIR, "gradient_buoyancy.npy")
        np.save(grad_k_file, grad_kappa)
        np.save(grad_b_file, grad_buoyancy)
        print(f"\n  Gradient (kappa) saved: {grad_k_file}")
        print(f"  Gradient (buoyancy) saved: {grad_b_file}")
        
        # Save adjoint seismogram if available
        if adjoint_seismo is not None:
            adj_seismo_file = os.path.join(OUTPUT_DIR, "adjoint_seismo.npy")
            np.save(adj_seismo_file, adjoint_seismo)
            print(f"  Adjoint seismogram saved: {adj_seismo_file}")
        
        # Delete gradient dictionary references to Kokkos arrays
        del gradient
        
        # =====================================================================
        # Summary
        # =====================================================================
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print(f"Observed seismogram: {obs_file}")
        print(f"Synthetic seismogram: {syn_file}")
        if snap_file:
            print(f"Wavefield snapshots: {snap_file}")
        print(f"Residuals: {res_file}")
        if adjoint_seismo is not None:
            adj_seismo_file = os.path.join(OUTPUT_DIR, "adjoint_seismo.npy")
            print(f"Adjoint seismogram: {adj_seismo_file}")
        print(f"Gradient (kappa): {grad_k_file}")
        print(f"Gradient (buoyancy): {grad_b_file}")
        print(f"\nFWI iteration completed successfully!")
        
    except Exception as e:
        print(f"\nError during FWI execution: {e}")
        import traceback
        traceback.print_exc()
        
    finally:
        # Clean up Kokkos objects BEFORE finalizing Kokkos
        print("\nCleaning up Kokkos objects...")
        
        # Clear field dictionaries (they contain numpy views of Kokkos arrays)
        try:
            if fields_fwd1:
                fields_fwd1.clear()
            if fields_fwd2:
                fields_fwd2.clear()
            if fields_adj:
                fields_adj.clear()
        except:
            pass
        
        # Delete objects in reverse order of creation
        try:
            del solver_adj, model_adj, data_adj, fields_adj
        except:
            pass
        
        try:
            del solver_fwd2, model_init, data_fwd2, fields_fwd2
        except:
            pass
        
        try:
            del solver_fwd1, model_true, data_fwd1, fields_fwd1
        except:
            pass
        
        # Force garbage collection to ensure all destructors are called
        gc.collect()
        
        # Now finalize Kokkos
        try:
            kokkos.finalize()
            print("Kokkos finalized successfully")
        except Exception as e:
            print(f"Kokkos finalization completed")


if __name__ == "__main__":
    main()
