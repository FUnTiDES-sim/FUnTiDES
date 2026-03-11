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


def _build_acoustic_solver(sd, sls_freqs=None, sls_coeffs=None):
    """Create an acoustic solver with optional SLS attenuation."""
    builder = Model.CartesianStructBuilder_f32_i32_O1(
        sd.ex, sd.hx, sd.ey, sd.hy, sd.ez, sd.hz, False, False
    )
    model = builder.get_model()

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


def _run_acoustic_simulation(solver, sd, n_steps=200, dt=0.001):
    """Run an acoustic simulation and return the L2 norm of the final field."""
    kk_pPrev, pPrev, kk_pCurr, pCurr = Utils.allocate_pressure(sd.n_dof)

    # Impulse at center node
    center = sd.n_dof // 2
    pCurr[center] = 1.0

    n_rhs = 1
    kk_rhsElem, _ = Utils.allocate_rhs_element(n_rhs, sd.ex, sd.ey, sd.ez)
    npp = (sd.order + 1) ** 3
    kk_rhsWeights = kokkos.array(
        [n_rhs, npp], dtype=kokkos.float32, space=kokkos.HostSpace, layout=kokkos.LayoutRight
    )
    ws = np.array(kk_rhsWeights, copy=False)
    ws[:] = 0.0

    kk_rhsTerm = kokkos.array(
        [n_rhs, n_steps], dtype=kokkos.float32, space=kokkos.HostSpace, layout=kokkos.LayoutRight
    )
    rt = np.array(kk_rhsTerm, copy=False)
    rt[:] = 0.0

    wavefield = Solver.WavefieldAcoustic(kk_pPrev, kk_pCurr)
    rhs = Solver.RhsAcoustic(kk_rhsTerm, kk_rhsElem, kk_rhsWeights)
    data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    for t in range(n_steps):
        solver.compute_forces(dt, t, data)
        solver.update_solution(dt, data)
        data.swap_wavefields()

    field = np.array(kk_pCurr, copy=False)
    return np.linalg.norm(field)


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
        """Simulation with attenuation should produce smaller amplitude."""
        sd = SmallDomain()
        n_steps = 200
        dt = 0.001

        # Run without attenuation
        solver_na, _ = _build_acoustic_solver(sd)
        norm_na = _run_acoustic_simulation(solver_na, sd, n_steps, dt)

        # Run with attenuation
        freqs = [2.0 * math.pi * 5.0, 2.0 * math.pi * 50.0]
        solver_att, _ = _build_acoustic_solver(sd, sls_freqs=freqs)
        norm_att = _run_acoustic_simulation(solver_att, sd, n_steps, dt)

        assert norm_na > 0.0, "Non-attenuated simulation should produce non-zero field"
        assert norm_att < norm_na, (
            f"Attenuated amplitude ({norm_att}) should be less than "
            f"non-attenuated ({norm_na})"
        )
