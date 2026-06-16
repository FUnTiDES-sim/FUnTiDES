"""
Tests for 3-buffer wavefield construction and swap() with backward mode.

The 3-buffer constructor enables backward/adjoint time-stepping by providing
prevprev, prev, and curr buffers. The swap() method automatically performs
3-way rotation when hasPrevPrev() returns True:

    curr ← prevPrev, prev ← curr, prevPrev ← prev

After one swap():
    - curr     holds what was prevPrev (ready for next write)
    - prev     holds what was curr (most recently computed field)
    - prevPrev holds what was prev (field from two steps ago)
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
# Acoustic 3-buffer tests
# ---------------------------------------------------------------------------


class TestWavefieldAcousticThreeBuffer:
    def setup_method(self):
        # prevprev = 3.0, prev = 1.0, curr = 2.0
        self.kk_prevprev = _alloc(3.0, "prevprev")
        self.kk_prev = _alloc(1.0, "prev")
        self.kk_curr = _alloc(2.0, "curr")

    def teardown_method(self):
        del self.kk_prevprev, self.kk_prev, self.kk_curr

    def test_three_buffer_constructor(self):
        """3-buffer constructor should set hasPrevPrev to True."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prevprev, self.kk_prev, self.kk_curr)
        assert wavefield.has_prevprev()

    def test_swap_performs_three_way_rotation(self):
        """swap() should perform 3-way rotation when hasPrevPrev is True."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prevprev, self.kk_prev, self.kk_curr)
        wavefield.swap()
        
        # After 3-way rotation: curr←prevprev(3), prev←curr(2), prevprev←prev(1)
        curr_field = wavefield.get_current_field(0)
        prev_field = wavefield.get_previous_field(0)
        prevprev_field = wavefield.get_prevprev_field(0)
        
        assert np.allclose(np.array(curr_field, copy=False), 3.0)
        assert np.allclose(np.array(prev_field, copy=False), 2.0)
        assert np.allclose(np.array(prevprev_field, copy=False), 1.0)

    def test_three_swaps_restore_state(self):
        """Three swaps on 3-buffer ring should restore original state."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prevprev, self.kk_prev, self.kk_curr)
        
        wavefield.swap()
        wavefield.swap()
        wavefield.swap()
        
        # Should be back to original: prevprev=3, prev=1, curr=2
        curr_field = wavefield.get_current_field(0)
        prev_field = wavefield.get_previous_field(0)
        prevprev_field = wavefield.get_prevprev_field(0)
        
        assert np.allclose(np.array(curr_field, copy=False), 2.0)
        assert np.allclose(np.array(prev_field, copy=False), 1.0)
        assert np.allclose(np.array(prevprev_field, copy=False), 3.0)

    def test_no_data_copy(self):
        """swap() should not copy data, only swap view handles."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prevprev, self.kk_prev, self.kk_curr)
        wavefield.swap()
        
        # Modify through prevprev field accessor
        prevprev_field = wavefield.get_prevprev_field(0)
        np.array(prevprev_field, copy=False)[0] = 999.0
        
        # Should be visible through original kk_prev buffer (which is now prevprev)
        assert np.array(self.kk_prev, copy=False)[0] == pytest.approx(999.0)


# ---------------------------------------------------------------------------
# Elastic 3-buffer tests  
# ---------------------------------------------------------------------------


class TestWavefieldElasticThreeBuffer:
    def setup_method(self):
        # ux: prevprev=3, prev=1, curr=2
        # uy: prevprev=6, prev=4, curr=5
        # uz: prevprev=9, prev=7, curr=8
        self.kk_ux_pp = _alloc(3.0, "ux_pp")
        self.kk_ux_prev = _alloc(1.0, "ux_prev")
        self.kk_ux_curr = _alloc(2.0, "ux_curr")
        self.kk_uy_pp = _alloc(6.0, "uy_pp")
        self.kk_uy_prev = _alloc(4.0, "uy_prev")
        self.kk_uy_curr = _alloc(5.0, "uy_curr")
        self.kk_uz_pp = _alloc(9.0, "uz_pp")
        self.kk_uz_prev = _alloc(7.0, "uz_prev")
        self.kk_uz_curr = _alloc(8.0, "uz_curr")

    def teardown_method(self):
        del self.kk_ux_pp, self.kk_ux_prev, self.kk_ux_curr
        del self.kk_uy_pp, self.kk_uy_prev, self.kk_uy_curr
        del self.kk_uz_pp, self.kk_uz_prev, self.kk_uz_curr

    def test_nine_argument_constructor(self):
        """9-arg constructor should enable 3-buffer mode."""
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_pp, self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_pp, self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_pp, self.kk_uz_prev, self.kk_uz_curr,
        )
        assert wavefield.has_prevprev()

    def test_swap_rotates_all_components(self):
        """swap() should rotate all three components correctly."""
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_pp, self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_pp, self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_pp, self.kk_uz_prev, self.kk_uz_curr,
        )
        wavefield.swap()
        
        # After rotation: curr←prevprev, prev←curr, prevprev←prev
        assert np.allclose(np.array(wavefield.get_current_field(0), copy=False), 3.0)  # ux
        assert np.allclose(np.array(wavefield.get_current_field(1), copy=False), 6.0)  # uy
        assert np.allclose(np.array(wavefield.get_current_field(2), copy=False), 9.0)  # uz
        
        assert np.allclose(np.array(wavefield.get_previous_field(0), copy=False), 2.0)  # ux
        assert np.allclose(np.array(wavefield.get_prevprev_field(0), copy=False), 1.0)  # ux

    def test_three_swaps_restore_state(self):
        """Three swaps should restore original state for all components."""
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_pp, self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_pp, self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_pp, self.kk_uz_prev, self.kk_uz_curr,
        )
        
        wavefield.swap()
        wavefield.swap()
        wavefield.swap()
        
        # Back to original state
        assert np.allclose(np.array(wavefield.get_current_field(0), copy=False), 2.0)
        assert np.allclose(np.array(wavefield.get_previous_field(0), copy=False), 1.0)
        assert np.allclose(np.array(wavefield.get_prevprev_field(0), copy=False), 3.0)


# ---------------------------------------------------------------------------
# Acoustic 2-buffer tests (forward mode)
# ---------------------------------------------------------------------------


class TestWavefieldAcousticTwoBuffer:
    def setup_method(self):
        # prev = 1.0, curr = 2.0
        self.kk_prev = _alloc(1.0, "prev")
        self.kk_curr = _alloc(2.0, "curr")

    def teardown_method(self):
        del self.kk_prev, self.kk_curr

    def test_two_buffer_constructor(self):
        """2-buffer constructor should set hasPrevPrev to False."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)
        assert not wavefield.has_prevprev()

    def test_swap_performs_two_way_swap(self):
        """swap() should perform 2-way swap when hasPrevPrev is False."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)
        wavefield.swap()
        
        # After 2-way swap: curr←prev(1), prev←curr(2)
        curr_field = wavefield.get_current_field(0)
        prev_field = wavefield.get_previous_field(0)
        
        assert np.allclose(np.array(curr_field, copy=False), 1.0)
        assert np.allclose(np.array(prev_field, copy=False), 2.0)

    def test_two_swaps_restore_state(self):
        """Two swaps on 2-buffer ring should restore original state."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)
        
        wavefield.swap()
        wavefield.swap()
        
        # Should be back to original: prev=1, curr=2
        curr_field = wavefield.get_current_field(0)
        prev_field = wavefield.get_previous_field(0)
        
        assert np.allclose(np.array(curr_field, copy=False), 2.0)
        assert np.allclose(np.array(prev_field, copy=False), 1.0)

    def test_no_data_copy_two_buffer(self):
        """swap() should not copy data in 2-buffer mode."""
        wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)
        wavefield.swap()
        
        # Modify through prev field accessor
        prev_field = wavefield.get_previous_field(0)
        np.array(prev_field, copy=False)[0] = 999.0
        
        # Should be visible through original kk_curr buffer (which is now prev)
        assert np.array(self.kk_curr, copy=False)[0] == pytest.approx(999.0)


# ---------------------------------------------------------------------------
# Elastic 2-buffer tests (forward mode)
# ---------------------------------------------------------------------------


class TestWavefieldElasticTwoBuffer:
    def setup_method(self):
        # ux: prev=1, curr=2
        # uy: prev=4, curr=5
        # uz: prev=7, curr=8
        self.kk_ux_prev = _alloc(1.0, "ux_prev")
        self.kk_ux_curr = _alloc(2.0, "ux_curr")
        self.kk_uy_prev = _alloc(4.0, "uy_prev")
        self.kk_uy_curr = _alloc(5.0, "uy_curr")
        self.kk_uz_prev = _alloc(7.0, "uz_prev")
        self.kk_uz_curr = _alloc(8.0, "uz_curr")

    def teardown_method(self):
        del self.kk_ux_prev, self.kk_ux_curr
        del self.kk_uy_prev, self.kk_uy_curr
        del self.kk_uz_prev, self.kk_uz_curr

    def test_six_argument_constructor(self):
        """6-arg constructor should enable 2-buffer mode."""
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
        )
        assert not wavefield.has_prevprev()

    def test_swap_swaps_all_components(self):
        """swap() should swap all three components correctly in 2-buffer mode."""
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
        )
        wavefield.swap()
        
        # After swap: curr←prev, prev←curr
        assert np.allclose(np.array(wavefield.get_current_field(0), copy=False), 1.0)  # ux
        assert np.allclose(np.array(wavefield.get_current_field(1), copy=False), 4.0)  # uy
        assert np.allclose(np.array(wavefield.get_current_field(2), copy=False), 7.0)  # uz
        
        assert np.allclose(np.array(wavefield.get_previous_field(0), copy=False), 2.0)  # ux
        assert np.allclose(np.array(wavefield.get_previous_field(1), copy=False), 5.0)  # uy
        assert np.allclose(np.array(wavefield.get_previous_field(2), copy=False), 8.0)  # uz

    def test_two_swaps_restore_state(self):
        """Two swaps should restore original state for all components."""
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
        )
        
        wavefield.swap()
        wavefield.swap()
        
        # Back to original state
        assert np.allclose(np.array(wavefield.get_current_field(0), copy=False), 2.0)
        assert np.allclose(np.array(wavefield.get_previous_field(0), copy=False), 1.0)
        assert np.allclose(np.array(wavefield.get_current_field(1), copy=False), 5.0)
        assert np.allclose(np.array(wavefield.get_previous_field(1), copy=False), 4.0)

