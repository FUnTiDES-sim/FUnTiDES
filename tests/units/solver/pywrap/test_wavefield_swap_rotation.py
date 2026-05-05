"""
Tests for WavefieldAcoustic.swap_with_rotation and
WavefieldElastic.swap_with_rotation.

The rotation performs a copy-free 3-way cyclic shift of Kokkos view handles:

    prevPrevBuffer ← prev ← curr ← prevPrevBuffer

After one call:
    - curr     holds the slot that was prevPrevBuffer (ready for next write)
    - prev     holds what was curr   (the most recently computed field)
    - prevPrev holds what was prev   (the field from two steps ago)
"""

import kokkos
import numpy as np
import pytest

import pyfuntides.solver as Solver

N = 50
MEM = kokkos.HostSpace
LAY = kokkos.LayoutRight


def _alloc(value, name):
    kk = kokkos.array([N], dtype=kokkos.float32, space=MEM, layout=LAY)
    np.array(kk, copy=False)[:] = value
    return kk


# ---------------------------------------------------------------------------
# Acoustic
# ---------------------------------------------------------------------------


class TestWavefieldAcousticSwapWithRotation:
    def setup_method(self):
        # prev = 1.0, curr = 2.0, prevprev = 3.0
        self.kk_prev    = _alloc(1.0, "prev")
        self.kk_curr    = _alloc(2.0, "curr")
        self.kk_prevprev = _alloc(3.0, "prevprev")

        self.wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)

    def teardown_method(self):
        # Explicitly delete views before session finalize to avoid deallocation warnings
        del self.kk_prev, self.kk_curr, self.kk_prevprev
        del self.wavefield

    def test_rotation_moves_prev_to_prevprev(self):
        # After rotation, the prevprev parameter now points to what was prev
        self.wavefield.swap_with_rotation(self.kk_prevprev)
        assert np.allclose(np.array(self.kk_prevprev, copy=False), 1.0)

    def test_three_rotations_restore_state(self):
        # Three rotations on a 3-element ring restores original state
        # After 3 rotations, kk_prevprev points back to its original buffer
        original_val = np.array(self.kk_prevprev, copy=False)[0]
        for _ in range(3):
            self.wavefield.swap_with_rotation(self.kk_prevprev)
        assert np.array(self.kk_prevprev, copy=False)[0] == pytest.approx(original_val)

    def test_no_data_copy_mutation_visible_through_rotated_handle(self):
        # After rotation kk_prevprev aliases the old kk_prev allocation.
        # Writing through one must be visible through the other.
        self.wavefield.swap_with_rotation(self.kk_prevprev)
        np.array(self.kk_prevprev, copy=False)[0] = 999.0
        # kk_prev is the underlying allocation that kk_prevprev now points to
        assert np.array(self.kk_prev, copy=False)[0] == pytest.approx(999.0)

    def test_swap_still_works_after_rotation(self):
        # Verify that swap() can be called after swap_with_rotation()
        # without crashing. We can't directly observe the internal state
        # from Python, but we can verify the operation completes.
        self.wavefield.swap_with_rotation(self.kk_prevprev)
        self.wavefield.swap()  # Should not crash
        # After rotation then swap, do another rotation to verify consistency
        self.wavefield.swap_with_rotation(self.kk_prevprev)
        # If we got here without crashing, the test passes


# ---------------------------------------------------------------------------
# Elastic
# ---------------------------------------------------------------------------


class TestWavefieldElasticSwapWithRotation:
    def setup_method(self):
        # ux: prev=1, curr=2, prevprev=3
        # uy: prev=4, curr=5, prevprev=6
        # uz: prev=7, curr=8, prevprev=9
        self.kk_ux_prev = _alloc(1.0, "ux_prev")
        self.kk_ux_curr = _alloc(2.0, "ux_curr")
        self.kk_uy_prev = _alloc(4.0, "uy_prev")
        self.kk_uy_curr = _alloc(5.0, "uy_curr")
        self.kk_uz_prev = _alloc(7.0, "uz_prev")
        self.kk_uz_curr = _alloc(8.0, "uz_curr")

        self.kk_ux_pp = _alloc(3.0, "ux_pp")
        self.kk_uy_pp = _alloc(6.0, "uy_pp")
        self.kk_uz_pp = _alloc(9.0, "uz_pp")

        self.wavefield = Solver.WavefieldElastic(
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
        )

    def teardown_method(self):
        # Explicitly delete views before session finalize to avoid deallocation warnings
        del self.kk_ux_prev, self.kk_ux_curr, self.kk_ux_pp
        del self.kk_uy_prev, self.kk_uy_curr, self.kk_uy_pp
        del self.kk_uz_prev, self.kk_uz_curr, self.kk_uz_pp
        del self.wavefield

    def test_rotation_correct_for_all_components(self):
        self.wavefield.swap_with_rotation(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )
        # After rotation, prevprev parameters point to what was prev
        assert np.allclose(np.array(self.kk_ux_pp, copy=False), 1.0)
        assert np.allclose(np.array(self.kk_uy_pp, copy=False), 4.0)
        assert np.allclose(np.array(self.kk_uz_pp, copy=False), 7.0)

    def test_three_rotations_restore_state(self):
        # Save original values
        orig_ux = np.array(self.kk_ux_pp, copy=False)[0]
        orig_uy = np.array(self.kk_uy_pp, copy=False)[0]
        orig_uz = np.array(self.kk_uz_pp, copy=False)[0]
        
        for _ in range(3):
            self.wavefield.swap_with_rotation(
                self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
            )
        
        # After 3 rotations, prevprev should point back to original buffers
        assert np.array(self.kk_ux_pp, copy=False)[0] == pytest.approx(orig_ux)
        assert np.array(self.kk_uy_pp, copy=False)[0] == pytest.approx(orig_uy)
        assert np.array(self.kk_uz_pp, copy=False)[0] == pytest.approx(orig_uz)

    def test_no_data_copy_ux(self):
        self.wavefield.swap_with_rotation(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )
        np.array(self.kk_ux_pp, copy=False)[0] = 999.0
        assert np.array(self.kk_ux_prev, copy=False)[0] == pytest.approx(999.0)
