"""
pytest entry point of the validation against closed-form solutions: one test per
case of cases.py. Full simulations, marked slow: run them with --run-slow, e.g.

    pytest tests/validation --run-slow -k "acoustic and order2"

Kokkos is initialized by tests/conftest.py. run_validation.py runs the same cases
without pytest.
"""

import pytest

pytest.importorskip("pyfuntides.solver")

from cases import all_cases
from run_validation import validate_case


@pytest.mark.slow
@pytest.mark.parametrize("case", all_cases(), ids=lambda c: c.name)
def test_validation(case, tmp_path):
    result = validate_case(case, tmp_path)
    assert result.total_error < case.tolerance, (
        f"{case.name}: relative L2 error {result.total_error:.3e} >= {case.tolerance:g} "
        f"(per component: {result.errors}), traces in {tmp_path / case.name}"
    )
