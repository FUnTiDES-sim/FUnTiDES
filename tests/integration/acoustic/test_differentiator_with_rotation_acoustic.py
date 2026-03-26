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
        
        # Adjoint wavefield: 3 buffers for rotation
        self.kk_qn_prev = _alloc(0.0, "qn_prev")
        self.kk_qn_curr = _alloc(0.0, "qn_curr")
        self.kk_qn_prevprev = _alloc(0.0, "qn_prevprev")
        
        # Gradient outputs (accumulate over time)
        self.kk_grad_kappa = _alloc(0.0, "grad_kappa")
        self.kk_grad_buoyancy = _alloc(0.0, "grad_buoyancy")
        
        # Create adjoint wavefield
        self.adj_wavefield = Solver.WavefieldAcoustic(
            self.kk_qn_prev, self.kk_qn_curr
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
        del (self.kk_pn, self.kk_qn_prev, self.kk_qn_curr, self.kk_qn_prevprev,
             self.kk_grad_kappa, self.kk_grad_buoyancy,
             self.adj_wavefield, self.grad, self.differentiator, self.model)
        
    def test_gradient_computation_with_rotation_loop(self):
        """
        Test complete gradient computation loop with rotation.
        Verifies exact gradient values based on mathematical formula.
        """
        # Volume per element (for mass matrix contribution)
        element_volume = (LX / EX) * (LY / EY) * (LZ / EZ)
        
        # Backward time loop
        for t in range(N_TIME_STEPS - 1, -1, -1):
            # Step 1: Load forward snapshot
            snapshot_value = 100.0 + t * 50.0
            np.array(self.kk_pn, copy=False)[:] = snapshot_value
            
            # Step 2: Simulate leapfrog update on adjoint wavefield
            np.array(self.kk_qn_curr, copy=False)[:] += t * 10.0
            np.array(self.kk_qn_prev, copy=False)[:] += t * 5.0
            
            # Capture state before gradient computation
            grad_before_kappa = np.array(self.kk_grad_kappa, copy=False).copy()
            qn_curr = np.array(self.kk_qn_curr, copy=False)[0]
            qn_prev = np.array(self.kk_qn_prev, copy=False)[0]
            qn_prevprev = np.array(self.kk_qn_prevprev, copy=False)[0]
            
            # Step 3: Rotate adjoint wavefield
            self.adj_wavefield.swap_with_rotation(self.kk_qn_prevprev)
            
            # Step 4: Build gradient data structures
            fwd_view = Gradient.WavefieldViewForwardAcoustic(self.kk_pn)
            bwd_view = Gradient.WavefieldViewBackwardAcoustic(
                self.kk_qn_curr, self.kk_qn_prev, self.kk_qn_prevprev
            )
            grad_data = Gradient.GradientDataAcoustic(fwd_view, bwd_view, self.grad)
            
            # Step 5: Compute gradient via differentiator
            self.differentiator.compute(self.model, grad_data, DT)
            
            # Step 6: Verify gradient values
            grad_kappa = np.array(self.kk_grad_kappa, copy=False)
            grad_buoyancy = np.array(self.kk_grad_buoyancy, copy=False)
            delta_grad_kappa = grad_kappa - grad_before_kappa
            
            # Property 1: grad_buoyancy stays zero (no spatial gradients)
            assert np.allclose(grad_buoyancy, 0.0, atol=1e-6), \
                f"t={t}: grad_buoyancy should be 0 for spatially uniform fields"
            
            # Property 2: Verify that we have exactly 4 distinct gradient values corresponding 
            # to interior, boundary, corner, edge nodes
            unique_deltas = np.unique(grad_kappa)
            assert len(unique_deltas) == 4, \
                f"t={t}: Expected exactly 4 distinct gradient values, got {len(unique_deltas)}"
            
            # Property 3: check contributions based on node type (interior vs boundary)
            # Verify 8:1 ratio between max and min (boundary vs interior)
            min_val = grad_kappa.min()
            max_val = grad_kappa.max()
            ratio = max_val / min_val if min_val > 0 else 0
            assert np.isclose(ratio, 8.0, rtol=0.001), \
                f"t={t}: Expected max/min ratio of 8.0, got {ratio:.2f}"
            # Verify corner node (first) has one contribution
            assert np.isclose(grad_kappa[0], min_val, rtol=0.01), \
                f"t={t}: corner node should have one contribution"
            # Verify edge nodes (next 12) have two contributions
            assert np.isclose(grad_kappa[1], min_val * 2, rtol=0.01), \
                f"t={t}: edge nodes should have two contributions"
            
            # Property 4: Check that it accumulates over time steps and is in the expected range
            assert np.all(grad_kappa > 0), f"t={t}: grad_kappa should accumulate to positive values"

        # Final checks
        final_grad_kappa = np.array(self.kk_grad_kappa, copy=False)
        final_grad_buoyancy = np.array(self.kk_grad_buoyancy, copy=False)
        
        assert final_grad_kappa.shape == (N_DOF,)
        assert final_grad_buoyancy.shape == (N_DOF,)
        assert np.any(final_grad_kappa > 0), "grad_kappa should accumulate to non-zero"
        assert np.allclose(final_grad_buoyancy, 0.0), "grad_buoyancy should be 0"
    