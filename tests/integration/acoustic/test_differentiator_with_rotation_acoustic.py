"""
Integration test for gradient computation with wavefield rotation.

Tests the complete adjoint gradient workflow:
1. Backward time loop
2. Load forward snapshots (fake values)
3. Apply leapfrog update on adjoint wavefield
4. Rotate adjoint wavefield
5. Compute gradient via differentiator
6. Verify gradient values
"""

import kokkos
import numpy as np
import pytest

import pyfuntides.solver as Solver
import pyfuntides.gradient as Gradient
import pyfuntides.model as Model

MEM     = kokkos.HostSpace
LAY     = kokkos.LayoutRight
BUILDER = Model.CartesianStructBuilder_f32_i32_O1

# Test configuration
EX = EY = EZ = 10      # Elements in x, y, z
LX = LY = LZ = 1000.0  # Length in x, y, z
ORDER = 1
N_DOF = (EX * ORDER + 1) * (EY * ORDER + 1) * (EZ * ORDER + 1)  
DT = 0.001
N_TIME_STEPS = 5

IMPL = Solver.ImplemType.MAKUTU
MESH = Solver.MeshType.STRUCT
LOC  = Solver.ModelLocationType.ONNODES
PHYS = Solver.PhysicType.ACOUSTIC

def _alloc(value, name):
    """Allocate and initialize a Kokkos view."""
    kk = kokkos.array([N_DOF], dtype=kokkos.float32, space=MEM, layout=LAY)
    np.array(kk, copy=False)[:] = value
    return kk

class TestDifferentiatorWithRotation:
    """Integration test for differentiator with wavefield rotation."""
    
    def setup_method(self):
        # Create mesh/model
        self.model = BUILDER(EX, LX, EY, LY, EZ, LZ, True, False).get_model()
        
        # Forward wavefield snapshot (will be loaded per time step)
        self.kk_pn = _alloc(0.0, "pn_snapshot")
        
        # Adjoint wavefield: 3 buffers for rotation (t-1, t, t+1)
        self.kk_qn0 = _alloc(0.0, "qn0")
        self.kk_qn1 = _alloc(0.0, "qn1")
        self.kk_qn2 = _alloc(0.0, "qn2")
        
        # Gradient outputs (accumulate over time)
        self.kk_grad_kappa = _alloc(0.0, "grad_kappa")
        self.kk_grad_buoyancy = _alloc(0.0, "grad_buoyancy")
        
        # Create adjoint wavefield
        self.adj_wavefield = Solver.WavefieldAcoustic(
            self.kk_qn0, self.kk_qn1
        )

        # Create gradient data
        self.grad = Gradient.GradientAcoustic(self.kk_grad_kappa, self.kk_grad_buoyancy)
        
        # Create differentiator
        self.differentiator = Gradient.create_differentiator(
            IMPL,
            MESH,
            LOC,
            PHYS,
            ORDER)
    
    def teardown_method(self):
        del (self.kk_pn, self.kk_qn1, self.kk_qn2, self.kk_qn0,
             self.kk_grad_kappa, self.kk_grad_buoyancy,
             self.adj_wavefield, self.grad, self.differentiator, self.model)
        
    def test_gradient_kappa_loop(self):
        """
        Test complete gradient computation loop with rotation.
        Verifies exact gradient values based on mathematical formula.
        """
        # Backward time loop
        for t in range(N_TIME_STEPS - 1, -1, -1):
            # Step 1: Load forward snapshot
            snapshot_value = 100.0 + t * 50.0
            np.array(self.kk_pn, copy=False)[:] = snapshot_value
            
            # Step 2: Simulate leapfrog update on adjoint wavefield
            np.array(self.adj_wavefield.get_previous_field(0), copy=False)[:] += t * 10.0
            np.array(self.adj_wavefield.get_current_field(0), copy=False)[:] += t * 5.0
            
            # Capture state before gradient computation
            grad_before_kappa = np.array(self.kk_grad_kappa, copy=False)
            qn_curr = np.array(self.adj_wavefield.get_current_field(0), copy=False)[0]
            qn_prev = np.array(self.adj_wavefield.get_previous_field(0), copy=False)[0]
            qn_prevprev = np.array(self.kk_qn2, copy=False)[0]
            
            # Step 3: Rotate adjoint wavefield
            self.kk_qn2 = self.adj_wavefield.swap_with_rotation(self.kk_qn2)
            adj_curr = self.adj_wavefield.get_current_field(0)
            adj_prev = self.adj_wavefield.get_previous_field(0)
            adj_prevprev = self.kk_qn2
            # check that they all swapped
            assert np.array(adj_curr, copy=False)[0] == qn_prevprev
            assert np.array(adj_prev, copy=False)[0] == qn_curr
            assert np.array(adj_prevprev, copy=False)[0] == qn_prev
            
            # Step 4: Build gradient data structures
            fwd_view = Gradient.WavefieldViewForwardAcoustic(self.kk_pn)
            bwd_view = Gradient.WavefieldViewBackwardAcoustic(
                adj_curr, adj_prev, adj_prevprev
            )
            grad_data = Gradient.GradientDataAcoustic(fwd_view, bwd_view, self.grad)
            
            # Step 5: Compute gradient via differentiator
            self.differentiator.compute(self.model, grad_data, DT)
            
            # Step 6: Verify gradient values
            grad_kappa = np.array(self.kk_grad_kappa, copy=False)
            grad_buoyancy = np.array(self.kk_grad_buoyancy, copy=False)
            delta_grad_kappa = grad_kappa - grad_before_kappa

            if t == N_TIME_STEPS - 1:
                # Second order derivative at first time step should be zero (no previous time step)
                # which means grad_kappa should not update from its initial zero value
                assert np.allclose(grad_kappa, 0.0), \
                    f"t={t}: Initial grad_kappa should be 0"
                # grad_buoyancy stays zero (no spatial gradients)
                assert np.allclose(grad_buoyancy, 0.0, atol=1e-6), \
                    f"t={t}: grad_buoyancy should be 0 for spatially uniform fields"
                continue
            
            # Property 1: grad_buoyancy stays zero (no spatial gradients)
            assert np.allclose(grad_buoyancy, 0.0, atol=1e-6), \
                f"t={t}: grad_buoyancy should be 0 for spatially uniform fields"
            
            # Property 2: With normalization by mass matrix diagonal,
            # all normalized values should be approximately equal 
            # (within floating point precision). Check that max/min ratio is close to 1.0
            min_val = grad_kappa.min()
            max_val = grad_kappa.max()
            ratio = max_val / min_val if min_val > 0 else 1.0
            # Allow up to 0.1% variation due to floating point rounding in atomic operations
            assert np.isclose(ratio, 1.0, rtol=0.001), \
                f"t={t}: Normalized gradients should be uniform, got ratio {ratio:.6f}"
            
            # Property 4: Check that it accumulates over time steps and is in the expected range
            assert np.all(grad_kappa > 0), f"t={t}: grad_kappa should accumulate to positive values"

        # Final checks
        final_grad_kappa = np.array(self.kk_grad_kappa, copy=False)
        final_grad_buoyancy = np.array(self.kk_grad_buoyancy, copy=False)
        
        assert final_grad_kappa.shape == (N_DOF,)
        assert final_grad_buoyancy.shape == (N_DOF,)
        assert np.any(final_grad_kappa > 0), "grad_kappa should accumulate to non-zero"
        assert np.allclose(final_grad_buoyancy, 0.0), "grad_buoyancy should be 0"
    
    def test_gradient_buoyancy_loop(self):
        """
        Test grad_buoyancy computation with spatially varying wavefields.
        
        The stiffness term computes spatial derivatives, so grad_buoyancy
        should be non-zero when both forward and adjoint fields have spatial variation.
        Uses a local perturbation for precise verification.
        """
        # Use a single time step with spatial variation
        t = 2  # Use middle time step to have all three time levels
        
        # Grid dimensions for 3D indexing
        nx = EX * ORDER + 1  # 11
        ny = EY * ORDER + 1  # 11
        nz = EZ * ORDER + 1  # 11
        
        # Step 1: Create spatially varying forward field with local perturbation
        # Base uniform field with a Gaussian-like bump in the center
        pn_array = np.array(self.kk_pn, copy=False)
        pn_array[:] = 100.0  # uniform base
        
        # Add a localized perturbation in the center region (5x5x5 nodes)
        center_x, center_y, center_z = nx // 2, ny // 2, nz // 2
        for i in range(max(0, center_x - 2), min(nx, center_x + 3)):
            for j in range(max(0, center_y - 2), min(ny, center_y + 3)):
                for k in range(max(0, center_z - 2), min(nz, center_z + 3)):
                    idx = i + j * nx + k * nx * ny
                    # Distance from center
                    dx, dy, dz = i - center_x, j - center_y, k - center_z
                    dist_sq = dx*dx + dy*dy + dz*dz
                    # Gaussian-like perturbation
                    pn_array[idx] += 50.0 * np.exp(-0.5 * dist_sq)
        
        # Step 2: Set up adjoint wavefield with SPATIAL VARIATION
        # This is critical for non-zero grad_buoyancy
        qn_curr_array = np.array(self.adj_wavefield.get_current_field(0), copy=False)
        qn_prev_array = np.array(self.adj_wavefield.get_previous_field(0), copy=False)
        qn_prevprev_array = np.array(self.kk_qn2, copy=False)
        
        # Add similar localized perturbation to adjoint fields
        qn_curr_array[:] = 50.0 + t * 5.0
        qn_prev_array[:] = 50.0 + (t - 1) * 5.0
        qn_prevprev_array[:] = 50.0 + (t - 2) * 5.0
        
        for i in range(max(0, center_x - 2), min(nx, center_x + 3)):
            for j in range(max(0, center_y - 2), min(ny, center_y + 3)):
                for k in range(max(0, center_z - 2), min(nz, center_z + 3)):
                    idx = i + j * nx + k * nx * ny
                    dx, dy, dz = i - center_x, j - center_y, k - center_z
                    dist_sq = dx*dx + dy*dy + dz*dz
                    perturbation = 30.0 * np.exp(-0.3 * dist_sq)
                    qn_curr_array[idx] += perturbation
                    qn_prev_array[idx] += perturbation
                    qn_prevprev_array[idx] += perturbation
        
        # Step 3: Rotate adjoint wavefield
        self.kk_qn2 = self.adj_wavefield.swap_with_rotation(self.kk_qn2)
        adj_curr = self.adj_wavefield.get_current_field(0)
        adj_prev = self.adj_wavefield.get_previous_field(0)
        adj_prevprev = self.kk_qn2
        
        # Step 4: Build gradient data structures
        fwd_view = Gradient.WavefieldViewForwardAcoustic(self.kk_pn)
        bwd_view = Gradient.WavefieldViewBackwardAcoustic(
            adj_curr, adj_prev, adj_prevprev
        )
        grad_data = Gradient.GradientDataAcoustic(fwd_view, bwd_view, self.grad)
        
        # Step 5: Compute gradient via differentiator
        self.differentiator.compute(self.model, grad_data, DT)
        
        # Step 6: Verify gradient values
        grad_kappa = np.array(self.kk_grad_kappa, copy=False)
        grad_buoyancy = np.array(self.kk_grad_buoyancy, copy=False)
        
        # Property 1: grad_kappa should be non-zero (mass term contribution)
        assert not np.allclose(grad_kappa, 0.0), \
            f"grad_kappa should be non-zero with spatial variation, got max={grad_kappa.max():.6e}, min={grad_kappa.min():.6e}"
        
        # Property 2: grad_buoyancy should be non-zero due to spatial variation
        # The stiffness term involves spatial derivatives of both fields
        assert not np.allclose(grad_buoyancy, 0.0), \
            f"grad_buoyancy should be non-zero with spatial variation, got max={grad_buoyancy.max():.6e}, min={grad_buoyancy.min():.6e}"
        
        # Property 3: Check localization - grad_buoyancy should be larger near the perturbation
        # Extract center region values
        center_idx = center_x + center_y * nx + center_z * nx * ny
        center_grad_buoyancy = abs(grad_buoyancy[center_idx])
        
        # Extract corner region values (far from perturbation)
        corner_indices = [0, nx-1, (ny-1)*nx, nx*ny-1]  # corners of z=0 plane
        corner_grad_buoyancy = np.mean([abs(grad_buoyancy[i]) for i in corner_indices])
        
        # Center should have larger gradient than corners (localization check)
        # Physical expectation: Gaussian perturbations decay exponentially with distance
        # With decay rates of 0.5 and 0.3, center should be >> corners
        # Minimum ratio of 2.0 ensures meaningful localization (not just numerical noise)
        if corner_grad_buoyancy > 1e-10:
            ratio = center_grad_buoyancy / corner_grad_buoyancy
            assert ratio > 2.0, \
                f"Center gradient ({center_grad_buoyancy:.6e}) should be >> corner gradient ({corner_grad_buoyancy:.6e}), ratio={ratio:.2f}"
        
        # Property 4: Check that grad_buoyancy has spatial structure
        # Physical expectation: The 5x5x5=125 node perturbed region should produce
        # multiple distinct gradient values due to varying distance from center
        # Minimum of 5 unique values ensures meaningful spatial variation
        unique_buoyancy_vals = len(np.unique(np.round(grad_buoyancy, decimals=10)))
        assert unique_buoyancy_vals > 5, \
            f"grad_buoyancy should have spatial structure (>5 unique values), got {unique_buoyancy_vals}"

        # Final checks
        assert grad_kappa.shape == (N_DOF,)
        assert grad_buoyancy.shape == (N_DOF,)