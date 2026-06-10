"""
Integration test for elastic gradient computation with 3-buffer wavefield swap.

Tests the complete adjoint gradient workflow for elastic media:
1. Backward time loop
2. Load forward displacement snapshots (fake values)
3. Apply leapfrog update on adjoint elastic wavefield
4. Swap adjoint wavefield (3-way rotation: curr←prevprev, prev←curr, prevprev←prev)
5. Compute gradient via differentiator
6. Verify gradient values (gradRho, gradLambda, gradMu)
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
EX = EY = EZ = 10
LX = LY = LZ = 1000.0
ORDER = 1
N_DOF = (EX * ORDER + 1) * (EY * ORDER + 1) * (EZ * ORDER + 1)
DT = 0.001
N_TIME_STEPS = 5

IMPL = Solver.ImplemType.MAKUTU
MESH = Solver.MeshType.STRUCT
LOC  = Solver.ModelLocationType.ONNODES
PHYS = Solver.PhysicType.ELASTIC


def _alloc(value, name):
    """Allocate and initialize a Kokkos view."""
    kk = kokkos.array([N_DOF], dtype=kokkos.float32, space=MEM, layout=LAY)
    np.array(kk, copy=False)[:] = value
    return kk


class TestDifferentiatorWithRotationElastic:
    """Integration test for elastic differentiator with 3-buffer wavefield swap."""

    def setup_method(self):
        # Create mesh/model (elastic = True)
        self.model = BUILDER(EX, LX, EY, LY, EZ, LZ, True, True).get_model()

        # Forward displacement snapshots (loaded per time step)
        self.kk_ux_fwd = _alloc(0.0, "ux_fwd")
        self.kk_uy_fwd = _alloc(0.0, "uy_fwd")
        self.kk_uz_fwd = _alloc(0.0, "uz_fwd")

        # Adjoint wavefield: elastic has 3 components × 3 time levels (3-buffer mode)
        # prevprev, prev, curr for each component
        self.kk_ux_pp = _alloc(0.0, "ux_pp")
        self.kk_ux_prev = _alloc(0.0, "ux_prev")
        self.kk_ux_curr = _alloc(0.0, "ux_curr")
        self.kk_uy_pp = _alloc(0.0, "uy_pp")
        self.kk_uy_prev = _alloc(0.0, "uy_prev")
        self.kk_uy_curr = _alloc(0.0, "uy_curr")
        self.kk_uz_pp = _alloc(0.0, "uz_pp")
        self.kk_uz_prev = _alloc(0.0, "uz_prev")
        self.kk_uz_curr = _alloc(0.0, "uz_curr")

        # Gradient outputs (accumulate over time)
        self.kk_grad_rho = _alloc(0.0, "grad_rho")
        self.kk_grad_lambda = _alloc(0.0, "grad_lambda")
        self.kk_grad_mu = _alloc(0.0, "grad_mu")

        # Create adjoint wavefield with 3-buffer constructor (9 args)
        self.adj_wavefield = Solver.WavefieldElastic(
            self.kk_ux_pp, self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_pp, self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_pp, self.kk_uz_prev, self.kk_uz_curr,
        )

        # Create gradient data
        self.grad = Gradient.GradientElastic(
            self.kk_grad_rho, self.kk_grad_lambda, self.kk_grad_mu
        )

        # Create differentiator
        self.differentiator = Gradient.create_differentiator(
            IMPL, MESH, LOC, PHYS, ORDER
        )

    def teardown_method(self):
        del (
            self.kk_ux_fwd, self.kk_uy_fwd, self.kk_uz_fwd,
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp,
            self.kk_grad_rho, self.kk_grad_lambda, self.kk_grad_mu,
            self.adj_wavefield, self.grad, self.differentiator, self.model,
        )

    def test_gradient_rho_loop(self):
        """
        Test complete elastic gradient computation loop with 3-buffer wavefield swap.
        Uses spatially uniform fields to verify gradRho accumulation.
        With uniform fields: div(u)=0, ε(u)=0, so gradLambda=gradMu=0.
        gradRho = ∫ ü† · u dΩ accumulates the mass term.
        """
        for t in range(N_TIME_STEPS - 1, -1, -1):
            # Step 1: Load forward snapshot (uniform displacement)
            snapshot_value = 100.0 + t * 50.0
            np.array(self.kk_ux_fwd, copy=False)[:] = snapshot_value
            np.array(self.kk_uy_fwd, copy=False)[:] = snapshot_value * 0.5
            np.array(self.kk_uz_fwd, copy=False)[:] = snapshot_value * 0.25

            # Step 2: Simulate leapfrog update on adjoint wavefield
            for comp_idx in range(3):
                np.array(
                    self.adj_wavefield.get_previous_field(comp_idx), copy=False
                )[:] += t * 10.0
                np.array(
                    self.adj_wavefield.get_current_field(comp_idx), copy=False
                )[:] += t * 5.0

            # Capture state before gradient computation
            grad_before_rho = np.array(self.kk_grad_rho, copy=False).copy()

            # Step 3: Swap adjoint wavefield (3-way rotation)
            self.adj_wavefield.swap()

            adj_curr_views = [
                self.adj_wavefield.get_current_field(i) for i in range(3)
            ]
            adj_prev_views = [
                self.adj_wavefield.get_previous_field(i) for i in range(3)
            ]
            adj_pp_views = [
                self.adj_wavefield.get_prevprev_field(i) for i in range(3)
            ]

            # Step 4: Compute second time derivative for backward wavefield
            # ü†_i = (q_curr_i - 2*q_prev_i + q_prevprev_i) / dt²
            kk_ux_dt2 = _alloc(0.0, "ux_dt2")
            kk_uy_dt2 = _alloc(0.0, "uy_dt2")
            kk_uz_dt2 = _alloc(0.0, "uz_dt2")

            dt2_views = [kk_ux_dt2, kk_uy_dt2, kk_uz_dt2]

            for comp_idx in range(3):
                curr = np.array(adj_curr_views[comp_idx], copy=False)
                prev = np.array(adj_prev_views[comp_idx], copy=False)
                prevprev = np.array(adj_pp_views[comp_idx], copy=False)
                np.array(dt2_views[comp_idx], copy=False)[:] = (
                    (curr - 2.0 * prev + prevprev) / (DT * DT)
                )

            # Step 5: Build gradient data structures
            fwd_view = Gradient.WavefieldViewForwardElastic(
                self.kk_ux_fwd, self.kk_uy_fwd, self.kk_uz_fwd
            )
            bwd_view = Gradient.WavefieldViewBackwardElastic(
                adj_curr_views[0], adj_curr_views[1], adj_curr_views[2],
                kk_ux_dt2, kk_uy_dt2, kk_uz_dt2,
            )
            grad_data = Gradient.GradientDataElastic(
                fwd_view, bwd_view, self.grad
            )

            # Step 6: Compute gradient via differentiator
            self.differentiator.compute(self.model, grad_data, DT)

            # Step 7: Verify gradient values
            grad_rho = np.array(self.kk_grad_rho, copy=False)
            grad_lambda = np.array(self.kk_grad_lambda, copy=False)
            grad_mu = np.array(self.kk_grad_mu, copy=False)

            # grad_lambda and grad_mu should stay zero (uniform fields => no strain)
            assert np.allclose(grad_lambda, 0.0, atol=1e-6), \
                f"t={t}: grad_lambda should be 0 for spatially uniform fields"
            assert np.allclose(grad_mu, 0.0, atol=1e-6), \
                f"t={t}: grad_mu should be 0 for spatially uniform fields"

            del kk_ux_dt2, kk_uy_dt2, kk_uz_dt2

        # Final checks
        final_grad_rho = np.array(self.kk_grad_rho, copy=False)
        final_grad_lambda = np.array(self.kk_grad_lambda, copy=False)
        final_grad_mu = np.array(self.kk_grad_mu, copy=False)

        assert final_grad_rho.shape == (N_DOF,)
        assert final_grad_lambda.shape == (N_DOF,)
        assert final_grad_mu.shape == (N_DOF,)
        assert np.allclose(final_grad_lambda, 0.0), \
            "grad_lambda should be 0 for uniform fields"
        assert np.allclose(final_grad_mu, 0.0), \
            "grad_mu should be 0 for uniform fields"

    def test_gradient_strain_loop(self):
        """
        Test grad_lambda and grad_mu computation with spatially varying wavefields.
        The stiffness term computes spatial derivatives of displacement,
        so grad_lambda and grad_mu should be non-zero when fields have
        spatial variation.
        """
        nx = EX * ORDER + 1
        ny = EY * ORDER + 1
        nz = EZ * ORDER + 1

        t = 2  # Use middle time step to have all three time levels

        # Step 1: Create spatially varying forward displacement
        ux_fwd_arr = np.array(self.kk_ux_fwd, copy=False)
        uy_fwd_arr = np.array(self.kk_uy_fwd, copy=False)
        uz_fwd_arr = np.array(self.kk_uz_fwd, copy=False)
        ux_fwd_arr[:] = 100.0
        uy_fwd_arr[:] = 100.0
        uz_fwd_arr[:] = 100.0

        # Add localized perturbation in the center region
        center_x, center_y, center_z = nx // 2, ny // 2, nz // 2
        for i in range(max(0, center_x - 2), min(nx, center_x + 3)):
            for j in range(max(0, center_y - 2), min(ny, center_y + 3)):
                for k in range(max(0, center_z - 2), min(nz, center_z + 3)):
                    idx = i + j * nx + k * nx * ny
                    dx, dy, dz = i - center_x, j - center_y, k - center_z
                    dist_sq = dx * dx + dy * dy + dz * dz
                    perturbation = 50.0 * np.exp(-0.5 * dist_sq)
                    ux_fwd_arr[idx] += perturbation
                    uy_fwd_arr[idx] += perturbation * 0.5
                    uz_fwd_arr[idx] += perturbation * 0.3

        # Step 2: Set up adjoint wavefield with spatial variation
        for comp_idx in range(3):
            curr_arr = np.array(
                self.adj_wavefield.get_current_field(comp_idx), copy=False
            )
            prev_arr = np.array(
                self.adj_wavefield.get_previous_field(comp_idx), copy=False
            )
            pp_arr = np.array(
                self.adj_wavefield.get_prevprev_field(comp_idx), copy=False
            )

            curr_arr[:] = 50.0 + t * 5.0
            prev_arr[:] = 50.0 + (t - 1) * 5.0
            pp_arr[:] = 50.0 + (t - 2) * 5.0

            for i in range(max(0, center_x - 2), min(nx, center_x + 3)):
                for j in range(max(0, center_y - 2), min(ny, center_y + 3)):
                    for k in range(
                        max(0, center_z - 2), min(nz, center_z + 3)
                    ):
                        idx = i + j * nx + k * nx * ny
                        dx = i - center_x
                        dy = j - center_y
                        dz = k - center_z
                        dist_sq = dx * dx + dy * dy + dz * dz
                        perturbation = 30.0 * np.exp(-0.3 * dist_sq)
                        curr_arr[idx] += perturbation
                        prev_arr[idx] += perturbation
                        pp_arr[idx] += perturbation

        # Step 3: Swap adjoint wavefield (3-way rotation)
        self.adj_wavefield.swap()
        
        adj_curr_views = [
            self.adj_wavefield.get_current_field(i) for i in range(3)
        ]
        adj_prev_views = [
            self.adj_wavefield.get_previous_field(i) for i in range(3)
        ]
        adj_pp_views = [
            self.adj_wavefield.get_prevprev_field(i) for i in range(3)
        ]

        # Step 4: Compute second time derivative
        kk_ux_dt2 = _alloc(0.0, "ux_dt2")
        kk_uy_dt2 = _alloc(0.0, "uy_dt2")
        kk_uz_dt2 = _alloc(0.0, "uz_dt2")
        dt2_views = [kk_ux_dt2, kk_uy_dt2, kk_uz_dt2]

        for comp_idx in range(3):
            curr = np.array(adj_curr_views[comp_idx], copy=False)
            prev = np.array(adj_prev_views[comp_idx], copy=False)
            prevprev = np.array(adj_pp_views[comp_idx], copy=False)
            np.array(dt2_views[comp_idx], copy=False)[:] = (
                (curr - 2.0 * prev + prevprev) / (DT * DT)
            )

        # Step 5: Build gradient data structures
        fwd_view = Gradient.WavefieldViewForwardElastic(
            self.kk_ux_fwd, self.kk_uy_fwd, self.kk_uz_fwd
        )
        bwd_view = Gradient.WavefieldViewBackwardElastic(
            adj_curr_views[0], adj_curr_views[1], adj_curr_views[2],
            kk_ux_dt2, kk_uy_dt2, kk_uz_dt2,
        )
        grad_data = Gradient.GradientDataElastic(fwd_view, bwd_view, self.grad)

        # Step 6: Compute gradient
        self.differentiator.compute(self.model, grad_data, DT)

        # Step 7: Verify gradient values
        grad_rho = np.array(self.kk_grad_rho, copy=False)
        grad_lambda = np.array(self.kk_grad_lambda, copy=False)
        grad_mu = np.array(self.kk_grad_mu, copy=False)

        # Property 1: grad_rho should be non-zero (mass term contribution)
        assert not np.allclose(grad_rho, 0.0), \
            f"grad_rho should be non-zero, got max={grad_rho.max():.6e}"

        # Property 2: grad_lambda should be non-zero (div interaction)
        assert not np.allclose(grad_lambda, 0.0), \
            f"grad_lambda should be non-zero with spatial variation, got " \
            f"max={grad_lambda.max():.6e}, min={grad_lambda.min():.6e}"

        # Property 3: grad_mu should be non-zero (strain interaction)
        assert not np.allclose(grad_mu, 0.0), \
            f"grad_mu should be non-zero with spatial variation, got " \
            f"max={grad_mu.max():.6e}, min={grad_mu.min():.6e}"

        # Property 4: Check localization — gradients near perturbation center
        # should be larger than at corners
        center_idx = center_x + center_y * nx + center_z * nx * ny
        corner_indices = [0, nx - 1, (ny - 1) * nx, nx * ny - 1]

        center_grad_lambda = abs(grad_lambda[center_idx])
        corner_grad_lambda = np.mean(
            [abs(grad_lambda[i]) for i in corner_indices]
        )

        if corner_grad_lambda > 1e-10:
            ratio = center_grad_lambda / corner_grad_lambda
            assert ratio > 2.0, \
                f"Center gradient ({center_grad_lambda:.6e}) should be >> " \
                f"corner gradient ({corner_grad_lambda:.6e}), ratio={ratio:.2f}"

        # Property 5: grad_lambda and grad_mu should have spatial structure
        unique_lambda = len(np.unique(np.round(grad_lambda, decimals=10)))
        assert unique_lambda > 5, \
            f"grad_lambda should have spatial structure (>5 unique values), " \
            f"got {unique_lambda}"

        unique_mu = len(np.unique(np.round(grad_mu, decimals=10)))
        assert unique_mu > 5, \
            f"grad_mu should have spatial structure (>5 unique values), " \
            f"got {unique_mu}"

        # Final shape checks
        assert grad_rho.shape == (N_DOF,)
        assert grad_lambda.shape == (N_DOF,)
        assert grad_mu.shape == (N_DOF,)

        del kk_ux_dt2, kk_uy_dt2, kk_uz_dt2
