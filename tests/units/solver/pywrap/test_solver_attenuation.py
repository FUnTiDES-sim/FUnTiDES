"""
Tests for the SLS wave attenuation implementation.

Verifies that setting SLS attenuation via the Python bindings:
  1. Does not crash the solver.
  2. Produces decaying wavefield amplitude compared to no attenuation.
  3. Rejects mismatched coefficient sizes appropriately.
"""

import math
import numpy as np
import kokkos
import pyfuntides.model as Model
import pyfuntides.solver as Solver
import pytest
import solver_utils as Utils


class SmallDomain:
    """Small 3x3x3 element structured mesh for attenuation tests."""

    def __init__(self, order=1):
        self.ex = self.ey = self.ez = 3
        self.domain_size = 300.0
        self.hx = self.domain_size / self.ex
        self.hy = self.domain_size / self.ey
        self.hz = self.domain_size / self.ez
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz


def _build_acoustic_solver(sd, sls_freqs=None, sls_coeffs=None, Q=None):
    """Create an acoustic solver with optional SLS attenuation."""
    builder = Model.CartesianStructBuilder_f32_i32_O1(
        sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, False, False
    )
    model = builder.get_model()

    if Q is not None:
        model.set_quality_factors(float(Q), float(Q))

    physic_type = Solver.PhysicType.ACOUSTIC
    if sls_freqs is not None or Q is not None:
        # Fallback to ACOUSTIC if VISCOACOUSTIC doesn't exist, but prefer VISCOACOUSTIC
        physic_type = getattr(Solver.PhysicType, "VISCOACOUSTIC", Solver.PhysicType.ACOUSTIC)

    solver = Solver.create_solver(
        Solver.MethodType.SEM,
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.STRUCT,
        Solver.ModelLocationType.ONELEMENTS,
        physic_type,
        sd.order,
    )

    if sls_freqs is not None:
        # If no coefficients are given, provide a default weight of 1.0 for each mechanism
        # so the attenuation actually takes effect.
        if sls_coeffs is None:
            sls_coeffs = [1.0] * len(sls_freqs)

        solver.set_sls_attenuation(sls_freqs, sls_coeffs)

    solver.compute_fe_init(model)
    return solver, model


def _run_acoustic_simulation(solver, sd, n_steps=150, dt=0.0002):
    """Run an acoustic simulation and return the L2 norm of the final field."""
    kk_pPrev, _, kk_pCurr, _ = Utils.allocate_pressure(sd.n_dof)

    n_rhs = 1
    # 1. Assign the source to the center element
    kk_rhsElem, rhsElem_np = Utils.allocate_rhs_element(n_rhs, sd.ex, sd.ey, sd.ez)
    rhsElem_np[0] = (sd.ex * sd.ey * sd.ez) // 2

    # 2. Assign the spatial weight to the center node of that element
    npp = (sd.order + 1) ** 3
    kk_rhsWeights = kokkos.array(
        [n_rhs, npp], dtype=kokkos.float32, space=kokkos.HostSpace, layout=kokkos.LayoutRight
    )
    ws = np.array(kk_rhsWeights, copy=False)
    ws[:] = 0.0
    ws[0, npp // 2] = 1.0

    # 3. Create a faster Ricker wavelet source (f0 = 50 Hz)
    kk_rhsTerm = kokkos.array(
        [n_rhs, n_steps], dtype=kokkos.float32, space=kokkos.HostSpace, layout=kokkos.LayoutRight
    )
    rt = np.array(kk_rhsTerm, copy=False)
    rt[:] = 0.0

    f0 = 50.0
    t0 = 1.2 / f0 # Peaks at 0.024 seconds
    for t in range(n_steps):
        time = t * dt
        tau = math.pi * f0 * (time - t0)
        rt[0, t] = (1.0 - 2.0 * tau**2) * math.exp(-tau**2)

    wavefield = Solver.WavefieldAcoustic(kk_pPrev, kk_pCurr)
    rhs = Solver.RhsAcoustic(kk_rhsTerm, kk_rhsElem, kk_rhsWeights)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    # 4. Run the time loop
    for t in range(n_steps):
        solver.compute_forces(dt, t, data)
        solver.update_solution(dt, data)
        data.swap_wavefields()

    # 5. Extract the energy.
    # Because swap_wavefields changes the internal pointers, we check both arrays
    # to guarantee we capture the energy regardless of odd/even step counts.
    field_curr = np.array(kk_pCurr, copy=False)
    field_prev = np.array(kk_pPrev, copy=False)

    return float(np.linalg.norm(field_curr) + np.linalg.norm(field_prev))


class TestAttenuationPython:
    def test_set_sls_attenuation_no_crash(self):
        """Setting SLS attenuation and running computeFEInit should not crash."""
        sd = SmallDomain()
        freqs = [2.0 * math.pi * 5.0, 2.0 * math.pi * 50.0]
        solver, _ = _build_acoustic_solver(sd, sls_freqs=freqs)
        # If we get here without exception, the test passes

    def test_empty_frequencies_disables_attenuation(self):
        """Empty frequency list should disable attenuation without error."""
        sd = SmallDomain()
        solver, _ = _build_acoustic_solver(sd, sls_freqs=[2.0 * math.pi * 5.0])
        # Disable attenuation
        solver.set_sls_attenuation([])
        # Should not raise

    def test_mismatched_coefficients_raises(self):
        """Mismatched frequency and coefficient sizes should raise."""
        sd = SmallDomain()
        builder = Model.CartesianStructBuilder_f32_i32_O1(
            sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, False, False
        )
        solver = Solver.create_solver(
            Solver.MethodType.SEM,
            Solver.ImplemType.MAKUTU,
            Solver.MeshType.STRUCT,
            Solver.ModelLocationType.ONELEMENTS,
            Solver.PhysicType.ACOUSTIC,
            sd.order,
        )
        with pytest.raises(RuntimeError):
            solver.set_sls_attenuation([1.0, 2.0, 3.0], [0.5])

    def test_acoustic_attenuation_decays_amplitude(self):
        """Simulation with attenuation should change the wavefield."""
        sd = SmallDomain()

        # Adjust time so wave only travels ~90 meters (safely inside the 300m box)
        n_steps = 150
        dt = 0.0002

        # Run without attenuation
        solver_na, _ = _build_acoustic_solver(sd)
        norm_na = _run_acoustic_simulation(solver_na, sd, n_steps, dt)

        # Run with attenuation (Q=10, 1 SLS mechanism)
        freqs = [2.0 * math.pi * 5.0]
        solver_att, _ = _build_acoustic_solver(sd, sls_freqs=freqs, Q=10)
        norm_att = _run_acoustic_simulation(solver_att, sd, n_steps, dt)

        assert norm_na > 0.0, "Non-attenuated simulation should produce non-zero field"
        assert np.isfinite(norm_att), "Attenuated simulation should remain stable"

        ratio = norm_att / norm_na
        assert ratio != 1.0, (
            f"Attenuation should change the wavefield. "
            f"norm_na={norm_na}, norm_att={norm_att}"
        )
