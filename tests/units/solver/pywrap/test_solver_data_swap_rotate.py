"""
Tests for SEMsolverDataAcoustic and SEMsolverDataElastic swap_wavefields()
and rotate_wavefields() methods.

swap_wavefields() exchanges the prev/curr view handles inside the solver data.

rotate_wavefields() performs a copy-free 3-way cyclic shift:

    prevPrevBuffer ← prev ← curr ← prevPrevBuffer

After rotate_wavefields(pp):
    - curr     holds the slot that was prevPrevBuffer  (ready for next write)
    - prev     holds what was curr   (the most recently computed field)
    - pp (returned) holds what was prev  (the field from two steps ago)

Since get_current_field / get_previous_field are not exposed on the Python
solver-data types, internal state is observed indirectly:

  - rotate_wavefields() returns the ejected prev handle, which aliases the
    original kk_prev allocation. Writing through it is visible via kk_prev.
  - After swap_wavefields(), the "prev" inside the solver data is the old
    curr allocation; a subsequent rotate therefore ejects a handle that
    aliases kk_curr.
"""

import kokkos
import numpy as np
import pytest

import pyfuntides.solver as Solver

N = 50
MEM = kokkos.HostSpace
LAY = kokkos.LayoutRight


def _vec(value):
    """Allocate a 1-D float32 Kokkos view filled with *value*."""
    kk = kokkos.array([N], dtype=kokkos.float32, space=MEM, layout=LAY)
    np.array(kk, copy=False)[:] = value
    return kk


def _array2d(rows, cols, value=0.0):
    """Allocate a 2-D float32 Kokkos view filled with *value*."""
    kk = kokkos.array([rows, cols], dtype=kokkos.float32, space=MEM, layout=LAY)
    np.array(kk, copy=False)[:] = value
    return kk


def _vec_int(n, value=0):
    """Allocate a 1-D int32 Kokkos view filled with *value*."""
    kk = kokkos.array([n], dtype=kokkos.int32, space=MEM, layout=LAY)
    np.array(kk, copy=False)[:] = value
    return kk


def _make_acoustic_rhs():
    """Minimal RhsAcoustic for use as a stub in solver data tests."""
    term = _array2d(1, 1)
    element = _vec_int(1)
    weights = _array2d(1, 1)
    return Solver.RhsAcoustic(term, element, weights)


def _make_elastic_rhs():
    """Minimal RhsElastic for use as a stub in solver data tests."""
    term = _array2d(1, 1)
    element = _vec_int(1)
    weights = _array2d(1, 1)
    return Solver.RhsElastic(term, term, term, element, weights)


# ---------------------------------------------------------------------------
# Acoustic — swap_wavefields
# ---------------------------------------------------------------------------


class TestSEMsolverDataAcousticSwap:
    def setup_method(self):
        # prev=1.0, curr=2.0, prevprev=3.0
        self.kk_prev = _vec(1.0)
        self.kk_curr = _vec(2.0)
        self.kk_pp   = _vec(3.0)
        wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)
        rhs = _make_acoustic_rhs()
        self.data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    def teardown_method(self):
        del self.kk_prev, self.kk_curr, self.kk_pp
        del self.data

    def test_swap_once_ejects_original_curr_on_rotate(self):
        # After swap: internal prev = kk_curr (2.0), internal curr = kk_prev (1.0).
        # A subsequent rotate ejects internal prev → returned handle aliases kk_curr.
        self.data.swap_wavefields()
        self.kk_pp = self.data.rotate_wavefields(self.kk_pp)
        assert np.allclose(np.array(self.kk_pp, copy=False), 2.0)

    def test_double_swap_ejects_original_prev_on_rotate(self):
        # Two swaps restore original state → rotate ejects kk_prev (1.0).
        self.data.swap_wavefields()
        self.data.swap_wavefields()
        self.kk_pp = self.data.rotate_wavefields(self.kk_pp)
        assert np.allclose(np.array(self.kk_pp, copy=False), 1.0)


# ---------------------------------------------------------------------------
# Acoustic — rotate_wavefields
# ---------------------------------------------------------------------------


class TestSEMsolverDataAcousticRotate:
    def setup_method(self):
        # prev=1.0, curr=2.0, prevprev=3.0
        self.kk_prev = _vec(1.0)
        self.kk_curr = _vec(2.0)
        self.kk_pp   = _vec(3.0)
        wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)
        rhs = _make_acoustic_rhs()
        self.data = Solver.SEMsolverDataAcoustic(wavefield, rhs)

    def teardown_method(self):
        del self.kk_prev, self.kk_curr, self.kk_pp
        del self.data

    def test_rotate_ejects_old_prev(self):
        # rotate returns the old prev handle (value=1.0, aliases kk_prev).
        self.kk_pp = self.data.rotate_wavefields(self.kk_pp)
        assert np.allclose(np.array(self.kk_pp, copy=False), 1.0)

    def test_three_rotations_restore_state(self):
        # Three rotations on the 3-element ring must restore the original
        # prevprev handle assignment.
        original_pp_val = np.array(self.kk_pp, copy=False)[0]  # 3.0
        for _ in range(3):
            self.kk_pp = self.data.rotate_wavefields(self.kk_pp)
        assert np.array(self.kk_pp, copy=False)[0] == pytest.approx(original_pp_val)

    def test_no_data_copy_mutation_visible_through_returned_handle(self):
        # After rotate, kk_pp aliases the old kk_prev allocation.
        # Writing through kk_pp must be visible via kk_prev.
        self.kk_pp = self.data.rotate_wavefields(self.kk_pp)
        np.array(self.kk_pp, copy=False)[0] = 999.0
        assert np.array(self.kk_prev, copy=False)[0] == pytest.approx(999.0)

    def test_swap_still_works_after_rotate(self):
        # Verify swap_wavefields() does not crash when called after rotate.
        self.kk_pp = self.data.rotate_wavefields(self.kk_pp)
        self.data.swap_wavefields()
        self.kk_pp = self.data.rotate_wavefields(self.kk_pp)


# ---------------------------------------------------------------------------
# Elastic — swap_wavefields
# ---------------------------------------------------------------------------


class TestSEMsolverDataElasticSwap:
    def setup_method(self):
        # ux: prev=1.0, curr=2.0  |  uy: prev=4.0, curr=5.0  |  uz: prev=7.0, curr=8.0
        self.kk_ux_prev = _vec(1.0)
        self.kk_ux_curr = _vec(2.0)
        self.kk_uy_prev = _vec(4.0)
        self.kk_uy_curr = _vec(5.0)
        self.kk_uz_prev = _vec(7.0)
        self.kk_uz_curr = _vec(8.0)
        self.kk_ux_pp = _vec(3.0)
        self.kk_uy_pp = _vec(6.0)
        self.kk_uz_pp = _vec(9.0)
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
        )
        rhs = _make_elastic_rhs()
        self.data = Solver.SEMsolverDataElastic(wavefield, rhs)

    def teardown_method(self):
        del self.kk_ux_prev, self.kk_ux_curr, self.kk_ux_pp
        del self.kk_uy_prev, self.kk_uy_curr, self.kk_uy_pp
        del self.kk_uz_prev, self.kk_uz_curr, self.kk_uz_pp
        del self.data

    def test_swap_once_ejects_original_curr_for_all_components(self):
        # After swap, internal prev = old curr for each component.
        # rotate therefore ejects the old curr handles.
        self.data.swap_wavefields()
        self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp = self.data.rotate_wavefields(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )
        assert np.allclose(np.array(self.kk_ux_pp, copy=False), 2.0)
        assert np.allclose(np.array(self.kk_uy_pp, copy=False), 5.0)
        assert np.allclose(np.array(self.kk_uz_pp, copy=False), 8.0)

    def test_double_swap_ejects_original_prev_for_all_components(self):
        # Two swaps restore original state → rotate ejects the old prev handles.
        self.data.swap_wavefields()
        self.data.swap_wavefields()
        self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp = self.data.rotate_wavefields(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )
        assert np.allclose(np.array(self.kk_ux_pp, copy=False), 1.0)
        assert np.allclose(np.array(self.kk_uy_pp, copy=False), 4.0)
        assert np.allclose(np.array(self.kk_uz_pp, copy=False), 7.0)


# ---------------------------------------------------------------------------
# Elastic — rotate_wavefields
# ---------------------------------------------------------------------------


class TestSEMsolverDataElasticRotate:
    def setup_method(self):
        # ux: prev=1.0, curr=2.0  |  uy: prev=4.0, curr=5.0  |  uz: prev=7.0, curr=8.0
        self.kk_ux_prev = _vec(1.0)
        self.kk_ux_curr = _vec(2.0)
        self.kk_uy_prev = _vec(4.0)
        self.kk_uy_curr = _vec(5.0)
        self.kk_uz_prev = _vec(7.0)
        self.kk_uz_curr = _vec(8.0)
        self.kk_ux_pp = _vec(3.0)
        self.kk_uy_pp = _vec(6.0)
        self.kk_uz_pp = _vec(9.0)
        wavefield = Solver.WavefieldElastic(
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
        )
        rhs = _make_elastic_rhs()
        self.data = Solver.SEMsolverDataElastic(wavefield, rhs)

    def teardown_method(self):
        del self.kk_ux_prev, self.kk_ux_curr, self.kk_ux_pp
        del self.kk_uy_prev, self.kk_uy_curr, self.kk_uy_pp
        del self.kk_uz_prev, self.kk_uz_curr, self.kk_uz_pp
        del self.data

    def test_rotate_ejects_old_prev_for_all_components(self):
        # rotate returns the old prev handle for each component.
        self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp = self.data.rotate_wavefields(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )
        assert np.allclose(np.array(self.kk_ux_pp, copy=False), 1.0)
        assert np.allclose(np.array(self.kk_uy_pp, copy=False), 4.0)
        assert np.allclose(np.array(self.kk_uz_pp, copy=False), 7.0)

    def test_three_rotations_restore_state(self):
        # Three rotations on the 3-element ring must restore the original
        # prevprev handle assignment for each component.
        orig = [
            np.array(self.kk_ux_pp, copy=False)[0],  # 3.0
            np.array(self.kk_uy_pp, copy=False)[0],  # 6.0
            np.array(self.kk_uz_pp, copy=False)[0],  # 9.0
        ]
        for _ in range(3):
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp = self.data.rotate_wavefields(
                self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
            )
        assert np.array(self.kk_ux_pp, copy=False)[0] == pytest.approx(orig[0])
        assert np.array(self.kk_uy_pp, copy=False)[0] == pytest.approx(orig[1])
        assert np.array(self.kk_uz_pp, copy=False)[0] == pytest.approx(orig[2])

    def test_no_data_copy_ux_mutation_visible_through_returned_handle(self):
        # After rotate, kk_ux_pp aliases the old kk_ux_prev allocation.
        # Writing through kk_ux_pp must be visible via kk_ux_prev.
        self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp = self.data.rotate_wavefields(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )
        np.array(self.kk_ux_pp, copy=False)[0] = 999.0
        assert np.array(self.kk_ux_prev, copy=False)[0] == pytest.approx(999.0)

    def test_swap_still_works_after_rotate(self):
        # Verify swap_wavefields() does not crash when called after rotate.
        self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp = self.data.rotate_wavefields(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )
        self.data.swap_wavefields()
        self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp = self.data.rotate_wavefields(
            self.kk_ux_pp, self.kk_uy_pp, self.kk_uz_pp
        )


# ---------------------------------------------------------------------------
# Isolation: rotating the wavefield object does NOT update the solver data
# ---------------------------------------------------------------------------
#
# SEMsolverData stores its wavefield BY VALUE (C++ copy at construction time).
# Rotating the Python-side WavefieldAcoustic / WavefieldElastic object after
# the solver data has been built only updates the Python object's internal
# handles — the solver data's own copy is left untouched.
#
# Observable consequence:
#   - wavefield.rotate(pp1) ejects the old kk_prev handle into pp1 and
#     replaces the wavefield's curr with kk_pp1_old.
#   - solver_data.rotate_wavefields(pp2) thereafter still ejects kk_prev
#     (value 1.0), NOT kk_curr (value 2.0) — proving the solver data was
#     blind to the wavefield-level rotation.
#   - Had the rotation leaked into the solver data the ejected value would
#     have been 2.0 (the old curr), not 1.0 (the old prev).


class TestAcousticWavefieldRotateDoesNotAffectSolverData:
    def setup_method(self):
        # prev=1.0, curr=2.0 — distinct values so we can tell which was ejected
        self.kk_prev = _vec(1.0)
        self.kk_curr = _vec(2.0)
        self.kk_pp1  = _vec(3.0)   # used for the wavefield-level rotate
        self.kk_pp2  = _vec(10.0)  # used for the solver-data-level rotate
        self.wavefield = Solver.WavefieldAcoustic(self.kk_prev, self.kk_curr)
        rhs = _make_acoustic_rhs()
        self.data = Solver.SEMsolverDataAcoustic(self.wavefield, rhs)

    def teardown_method(self):
        del self.kk_prev, self.kk_curr, self.kk_pp1, self.kk_pp2
        del self.wavefield, self.data

    def test_wavefield_rotate_does_not_update_solver_data_prev_pointer(self):
        # Rotate via the *wavefield* object — solver data must be unaffected.
        self.kk_pp1 = self.wavefield.rotate(self.kk_pp1)
        # kk_pp1 now aliases kk_prev (value 1.0); wavefield's curr → old kk_pp1 (3.0)

        # Now rotate via the *solver data* — it must still eject kk_prev (1.0),
        # not kk_curr (2.0), because the solver data's copy was never changed.
        self.kk_pp2 = self.data.rotate_wavefields(self.kk_pp2)
        assert np.allclose(np.array(self.kk_pp2, copy=False), 1.0), (
            "solver data ejected the wrong prev: wavefield.rotate() must not "
            "update the solver data's internal view handles"
        )

    def test_wavefield_rotate_changes_wavefield_but_not_solver_data(self):
        # After wavefield.rotate(), the wavefield's curr now points to what
        # was kk_pp1 (value 3.0).  If the solver data had been updated, its
        # next rotate_wavefields would eject kk_curr (2.0) instead of kk_prev (1.0).
        self.kk_pp1 = self.wavefield.rotate(self.kk_pp1)
        # Confirm the wavefield DID change: pp1 now has the old prev value (1.0).
        assert np.allclose(np.array(self.kk_pp1, copy=False), 1.0)

        # Confirm the solver data did NOT change: it still ejects kk_prev (1.0).
        self.kk_pp2 = self.data.rotate_wavefields(self.kk_pp2)
        assert np.allclose(np.array(self.kk_pp2, copy=False), 1.0)
        # The solver data's ejected value is 1.0 (kk_prev), not 2.0 (kk_curr).
        assert not np.allclose(np.array(self.kk_pp2, copy=False), 2.0)


class TestElasticWavefieldRotateDoesNotAffectSolverData:
    def setup_method(self):
        # ux: prev=1.0, curr=2.0  |  uy: prev=4.0, curr=5.0  |  uz: prev=7.0, curr=8.0
        self.kk_ux_prev = _vec(1.0)
        self.kk_ux_curr = _vec(2.0)
        self.kk_uy_prev = _vec(4.0)
        self.kk_uy_curr = _vec(5.0)
        self.kk_uz_prev = _vec(7.0)
        self.kk_uz_curr = _vec(8.0)
        # pp1 used for wavefield-level rotate, pp2 used for solver-data-level rotate
        self.kk_ux_pp1 = _vec(3.0)
        self.kk_uy_pp1 = _vec(6.0)
        self.kk_uz_pp1 = _vec(9.0)
        self.kk_ux_pp2 = _vec(10.0)
        self.kk_uy_pp2 = _vec(20.0)
        self.kk_uz_pp2 = _vec(30.0)
        self.wavefield = Solver.WavefieldElastic(
            self.kk_ux_prev, self.kk_ux_curr,
            self.kk_uy_prev, self.kk_uy_curr,
            self.kk_uz_prev, self.kk_uz_curr,
        )
        rhs = _make_elastic_rhs()
        self.data = Solver.SEMsolverDataElastic(self.wavefield, rhs)

    def teardown_method(self):
        del self.kk_ux_prev, self.kk_ux_curr
        del self.kk_uy_prev, self.kk_uy_curr
        del self.kk_uz_prev, self.kk_uz_curr
        del self.kk_ux_pp1, self.kk_uy_pp1, self.kk_uz_pp1
        del self.kk_ux_pp2, self.kk_uy_pp2, self.kk_uz_pp2
        del self.wavefield, self.data

    def test_wavefield_rotate_does_not_update_solver_data_prev_pointers(self):
        # Rotate via the *wavefield* object — solver data must be unaffected.
        self.kk_ux_pp1, self.kk_uy_pp1, self.kk_uz_pp1 = self.wavefield.rotate(
            self.kk_ux_pp1, self.kk_uy_pp1, self.kk_uz_pp1
        )
        # pp1 handles now alias the old prev allocations (1.0, 4.0, 7.0).
        # The solver data should still eject the old prev handles (same values).
        self.kk_ux_pp2, self.kk_uy_pp2, self.kk_uz_pp2 = self.data.rotate_wavefields(
            self.kk_ux_pp2, self.kk_uy_pp2, self.kk_uz_pp2
        )
        assert np.allclose(np.array(self.kk_ux_pp2, copy=False), 1.0), (
            "ux: solver data ejected the wrong prev"
        )
        assert np.allclose(np.array(self.kk_uy_pp2, copy=False), 4.0), (
            "uy: solver data ejected the wrong prev"
        )
        assert np.allclose(np.array(self.kk_uz_pp2, copy=False), 7.0), (
            "uz: solver data ejected the wrong prev"
        )

    def test_wavefield_rotate_changes_wavefield_but_not_solver_data(self):
        # Confirm wavefield DID rotate (pp1 has old prev values).
        self.kk_ux_pp1, self.kk_uy_pp1, self.kk_uz_pp1 = self.wavefield.rotate(
            self.kk_ux_pp1, self.kk_uy_pp1, self.kk_uz_pp1
        )
        assert np.allclose(np.array(self.kk_ux_pp1, copy=False), 1.0)
        assert np.allclose(np.array(self.kk_uy_pp1, copy=False), 4.0)
        assert np.allclose(np.array(self.kk_uz_pp1, copy=False), 7.0)

        # Confirm solver data did NOT change: it still ejects old prev (not old curr).
        self.kk_ux_pp2, self.kk_uy_pp2, self.kk_uz_pp2 = self.data.rotate_wavefields(
            self.kk_ux_pp2, self.kk_uy_pp2, self.kk_uz_pp2
        )
        # Must be old prev values (1.0, 4.0, 7.0), NOT old curr values (2.0, 5.0, 8.0).
        assert np.allclose(np.array(self.kk_ux_pp2, copy=False), 1.0)
        assert np.allclose(np.array(self.kk_uy_pp2, copy=False), 4.0)
        assert np.allclose(np.array(self.kk_uz_pp2, copy=False), 7.0)
        assert not np.allclose(np.array(self.kk_ux_pp2, copy=False), 2.0)
        assert not np.allclose(np.array(self.kk_uy_pp2, copy=False), 5.0)
        assert not np.allclose(np.array(self.kk_uz_pp2, copy=False), 8.0)
