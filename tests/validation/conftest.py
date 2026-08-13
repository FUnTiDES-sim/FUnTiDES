import os
import kokkos
import pytest


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


def pytest_addoption(parser):
    """
    Parses command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """

    parser.addoption(
        "--debug",
        action="store_true",
        default=False,
        help=f"Activate debug to have more prints (default: False)",
    )
    parser.addoption(
        "--mem",
        choices=["CPU", "GPU"],
        default="default_mem",
        help=f"Choose Kokkos memspace: CPU, GPU (default: [auto-detected])",
    )
    parser.addoption(
        "--model",
        choices=["STRUCTURED", "UNSTRUCTURED"],
        default="STRUCTURED",
        help=f"Choose model type: STRUCTURED, UNSTRUCTURED (default: STRUCTURED)",
    )
    parser.addoption(
        "--impl",
        choices=["CLASSIC", "MAKUTU", "OPTIM", "SHIVA"],
        default="MAKUTU",
        help=f"Choose implementation type: CLASSIC, MAKUTU, OPTIM, SHIVA (default: MAKUTU)",
    )
    parser.addoption(
        "--order",
        type=int,
        default=2,
        choices=range(1, 4),
        help="Polynomial order of the elements (default: 2, max 3)",
    )
    parser.addoption(
        "--domain_size",
        type=float,
        default=1500.0,
        help="Size of the cubic domain (default: 1500.0)",
    )
    parser.addoption(
        "--ex",
        type=int,
        default=50,
        help="Number of elements in x-direction (default: 50)",
    )
    parser.addoption(
        "--ey",
        type=int,
        default=50,
        help="Number of elements in y-direction (default: 50)",
    )
    parser.addoption(
        "--ez",
        type=int,
        default=50,
        help="Number of elements in z-direction (default: 50)",
    )
    parser.addoption(
        "--f0",
        type=float,
        default=5.0,
        help="Peak frequency for the Ricker source term (default: 5.0)",
    )
    parser.addoption(
        "--dt",
        type=float,
        default=0.001,
        help="Time step size (default: 0.001)",
    )
    parser.addoption(
        "--n_time_steps",
        type=int,
        default=1500,
        help="Number of time steps to run (default: 1500)",
    )
    parser.addoption(
        "--n_rhs",
        type=int,
        default=1,
        help="Number of right-hand side sources (default: 1)",
    )
    parser.addoption(
        "--on_nodes",
        action="store_true",
        default=False,
        help="Whether to apply model on nodes (default: False)",
    )
    parser.addoption(
        "--n_rcv",
        type=int,
        default=1,
        help="Number of receivers (default: 1)",
    )
    parser.addoption(
        "--is_elastic",
        action="store_true",
        default=False,
        help="Solving Elastic wave equation (True) or Acoustic wave equation (False) (default: False)",
    )
    parser.addoption(
        "--is_backward",
        action="store_true",
        default=False,
        help="Reverse time step by taking -dt (default: False)",
    )


@pytest.fixture
def on_nodes(request):
    return request.config.getoption("--on_nodes")

@pytest.fixture
def is_elastic(request):
    return request.config.getoption("--is_elastic")

@pytest.fixture
def is_backward(request):
    return request.config.getoption("--is_backward")

@pytest.fixture
def dt(request):
    return request.config.getoption("--dt")

@pytest.fixture
def f0(request):
    return request.config.getoption("--f0")

@pytest.fixture
def n_time_steps(request):
    return request.config.getoption("--n_time_steps")

@pytest.fixture
def n_rhs(request):
    return request.config.getoption("--n_rhs")

@pytest.fixture
def n_rcv(request):
    return request.config.getoption("--n_rcv")

@pytest.fixture
def order(request):
    return request.config.getoption("--order")

@pytest.fixture
def domain_size(request):
    return request.config.getoption("--domain_size")

@pytest.fixture
def ex(request):
    return request.config.getoption("--ex")

@pytest.fixture
def ey(request):
    return request.config.getoption("--ey")

@pytest.fixture
def ez(request):
    return request.config.getoption("--ez")

@pytest.fixture
def mem(request):
    return request.config.getoption("--mem")

@pytest.fixture
def impl(request):
    return request.config.getoption("--impl")

@pytest.fixture
def model(request):
    return request.config.getoption("--model")
