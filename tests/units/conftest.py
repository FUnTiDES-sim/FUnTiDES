import kokkos
import os


# Avoid using the entire developer node for tests
os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OMP_THREAD_LIMIT", "6")
os.environ.setdefault("KOKKOS_NUM_THREADS", "6")


def pytest_sessionstart(session):
    kokkos.initialize()


def pytest_sessionfinish(session, exitstatus):
    kokkos.finalize()
