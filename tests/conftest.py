import os
import kokkos
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--kokkos-num-threads",
        type=int,
        default=6,
        help="Number of Kokkos threads"
    )
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run tests marked slow (full simulations, skipped by default)",
    )

    # Options of tests/integration/test_self_adjoint.py
    group = parser.getgroup("self-adjoint")
    group.addoption(
        "--sa-debug",
        action="store_true",
        default=False,
        help="Activate debug to have more prints and snapshots (default: False)",
    )
    group.addoption(
        "--mem",
        choices=["CPU", "GPU"],
        default="default_mem",
        help="Choose Kokkos memspace: CPU, GPU (default: [auto-detected])",
    )
    group.addoption(
        "--model",
        choices=["STRUCTURED", "UNSTRUCTURED"],
        default="STRUCTURED",
        help="Choose model type: STRUCTURED, UNSTRUCTURED (default: STRUCTURED)",
    )
    group.addoption(
        "--impl",
        choices=["CLASSIC", "MAKUTU", "OPTIM", "SHIVA"],
        default="MAKUTU",
        help="Choose implementation type: CLASSIC, MAKUTU, OPTIM, SHIVA (default: MAKUTU)",
    )
    group.addoption(
        "--order",
        type=int,
        default=2,
        choices=range(1, 4),
        help="Polynomial order of the elements (default: 2, max 3)",
    )
    group.addoption(
        "--domain_size",
        type=float,
        default=1500.0,
        help="Size of the cubic domain (default: 1500.0)",
    )
    group.addoption(
        "--ex",
        type=int,
        default=50,
        help="Number of elements in x-direction (default: 50)",
    )
    group.addoption(
        "--ey",
        type=int,
        default=50,
        help="Number of elements in y-direction (default: 50)",
    )
    group.addoption(
        "--ez",
        type=int,
        default=50,
        help="Number of elements in z-direction (default: 50)",
    )
    group.addoption(
        "--f0",
        type=float,
        default=5.0,
        help="Peak frequency for the Ricker source term (default: 5.0)",
    )
    group.addoption(
        "--dt",
        type=float,
        default=0.001,
        help="Time step size (default: 0.001)",
    )
    group.addoption(
        "--n_time_steps",
        type=int,
        default=1500,
        help="Number of time steps to run (default: 1500)",
    )
    group.addoption(
        "--n_rhs",
        type=int,
        default=1,
        help="Number of right-hand side sources (default: 1)",
    )
    group.addoption(
        "--on_nodes",
        action="store_true",
        default=False,
        help="Whether to apply model on nodes (default: False)",
    )
    group.addoption(
        "--n_rcv",
        type=int,
        default=1,
        help="Number of receivers (default: 1)",
    )
    group.addoption(
        "--is_elastic",
        action="store_true",
        default=False,
        help="Solving Elastic wave equation (True) or Acoustic wave equation (False) (default: False)",
    )
    group.addoption(
        "--is_backward",
        action="store_true",
        default=False,
        help="Reverse time step by taking -dt (default: False)",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: full simulation, run only with --run-slow"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip_slow = pytest.mark.skip(reason="slow test, use --run-slow to run it")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


def _set_thread_env(n: int):
    v = str(n)
    for var in ("KOKKOS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT"):
        os.environ[var] = v


def pytest_sessionstart(session):
    n = session.config.getoption("--kokkos-num-threads")
    _set_thread_env(n)
    kokkos.initialize()


def pytest_sessionfinish(session, exitstatus):
    kokkos.finalize()
