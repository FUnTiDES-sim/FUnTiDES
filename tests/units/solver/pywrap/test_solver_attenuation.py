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

    solver = Solver.create_solver(
        Solver.MethodType.SEM,
        Solver.ImplemType.MAKUTU,
        Solver.MeshType.STRUCT,
        Solver.ModelLocationType.ONELEMENTS,
        Solver.PhysicType.ACOUSTIC,
        sd.order,
    )

    if sls_freqs is not None:
        if sls_coeffs is not None:
            solver.set_sls_attenuation(sls_freqs, sls_coeffs)
        else:
            solver.set_sls_attenuation(sls_freqs)

    solver.compute_fe_init(model)
    return solver, model


def _build_rhs(sd, n_steps, dt):
    """Build a Ricker wavelet RHS source at the center element."""
    n_rhs = 1
    kk_rhsElem, rhsElem_np = Utils.allocate_rhs_element(n_rhs, sd.ex, sd.ey, sd.ez)
    rhsElem_np[0] = (sd.ex * sd.ey * sd.ez) // 2

    npp = (sd.order + 1) ** 3
    kk_rhsWeights = kokkos.array(
        [n_rhs, npp], dtype=kokkos.float32, layout=kokkos.LayoutRight
    )
    ws = np.array(kk_rhsWeights, copy=False)
    ws[:] = 0.0
    ws[0, npp // 2] = 1.0

    kk_rhsTerm = kokkos.array(
        [n_rhs, n_steps], dtype=kokkos.float32, layout=kokkos.LayoutRight
    )
    rt = np.array(kk_rhsTerm, copy=False)
    rt[:] = 0.0

    f0 = 5.0
    t0 = 1.2 / f0
    amplitude = 1.0e8
    for t in range(n_steps):
        time = t * dt
        tau = math.pi * f0 * (time - t0)
        rt[0, t] = amplitude * (1.0 - 2.0 * tau**2) * math.exp(-(tau**2))

    return kk_rhsElem, kk_rhsWeights, kk_rhsTerm


def _run_acoustic_simulation(solver, sd, n_steps=500, dt=0.001):
    """Run an acoustic simulation and return the L2 norm of the final pressure field.

    Note: WavefieldAcoustic deep-copies the Python buffers onto the device at
    construction time. Results must therefore be read back via the wavefield's
    pCurr/pPrev properties, which perform an explicit device-to-host copy.

    swap_wavefields() rotates the internal C++ pointers each step, so the buffer
    holding the final written field depends on the parity of n_steps:
      - even: last written field is in pCurr (m_pnGlobalCurr)
      - odd:  last written field has rotated into pPrev (m_pnGlobalPrev)
    """
    kk_pPrev, pPrev, kk_pCurr, pCurr = Utils.allocate_pressure(sd.n_dof)
    kk_rhsElem, kk_rhsWeights, kk_rhsTerm = _build_rhs(sd, n_steps, dt)

    wavefield = Solver.WavefieldAcoustic(kk_pPrev, kk_pCurr)
    rhs = Solver.RhsAcoustic(kk_rhsTerm, kk_rhsElem, kk_rhsWeights)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    for t in range(n_steps):
        solver.compute_forces(dt, t, data)
        solver.update_solution(dt, data)
        data.swap_wavefields()

    field = wavefield.pCurr if n_steps % 2 == 0 else wavefield.pPrev
    return float(np.linalg.norm(field))


class TestAttenuationPython:
    def test_set_sls_attenuation_no_crash(self):
        """Setting SLS attenuation and running compute_fe_init should not crash."""
        sd = SmallDomain()
        freqs = [2.0 * math.pi * 5.0, 2.0 * math.pi * 50.0]
        solver, _ = _build_acoustic_solver(sd, sls_freqs=freqs)

    def test_empty_frequencies_disables_attenuation(self):
        """Empty frequency list should disable attenuation without error."""
        sd = SmallDomain()
        solver, _ = _build_acoustic_solver(sd, sls_freqs=[2.0 * math.pi * 5.0])
        solver.set_sls_attenuation([])

    def test_mismatched_coefficients_raises(self):
        """Mismatched frequency and coefficient sizes should raise a RuntimeError."""
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
        """Simulation with SLS attenuation should produce a smaller L2 norm."""
        sd = SmallDomain()
        n_steps = 500
        dt = 0.001

        solver_na, _ = _build_acoustic_solver(sd)
        norm_na = _run_acoustic_simulation(solver_na, sd, n_steps, dt)

        f0_rad = 2.0 * math.pi * 5.0
        sls_coeff = 1.0 / 10.0  # ≈ 1/Q for a single-mechanism approximation
        solver_att, _ = _build_acoustic_solver(
            sd, sls_freqs=[f0_rad], sls_coeffs=[sls_coeff], Q=10
        )
        norm_att = _run_acoustic_simulation(solver_att, sd, n_steps, dt)

        assert norm_na > 0.0, "Non-attenuated simulation produced zero field."
        assert np.isfinite(norm_att), "Attenuated simulation diverged."
        assert norm_att < norm_na, (
            f"Attenuation should reduce wavefield amplitude. "
            f"norm_na={norm_na:.4g}, norm_att={norm_att:.4g}"
        )
