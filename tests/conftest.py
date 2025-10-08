import os
import kokkos


def pytest_addoption(parser):
    parser.addoption(
        "--threads",
        type=int,
        default=6,
        help="Number of Kokkos threads (sets KOKKOS_NUM_THREADS / OMP_NUM_THREADS / OMP_THREAD_LIMIT before kokkos.initialize())"
    )


def _set_thread_env(n: int):
    v = str(n)
    for var in ("KOKKOS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT"):
        os.environ[var] = v


def pytest_sessionstart(session):
    n = session.config.getoption("--threads")
    _set_thread_env(n)
    kokkos.initialize()


def pytest_sessionfinish(session, exitstatus):
    kokkos.finalize()
