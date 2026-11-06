from importlib.metadata import version

__version__ = version("pyfuntides")

from pyfuntides import model,solver,gradient,io

from pyfuntides._build_info import METHODS,MAX_SOLVER_ORDER,MAX_DIFFERENTIATOR_ORDER

import kokkos
from pyfuntides import _build_info

_built_with_uvm = "CUDA_UVM" in _build_info.KOKKOS_OPTIONS
_running_with_uvm = kokkos.DefaultMemorySpace == kokkos.CudaUVMSpace

if _built_with_uvm != _running_with_uvm:
    raise ImportError(
        f"pyfuntides was built against a Kokkos with CUDA UVM={_built_with_uvm}, "
        f"but the loaded kokkos module runs with CUDA UVM={_running_with_uvm}. "
        "Reinstall pyfuntides with Kokkos_ROOT pointing to the Kokkos used by pykokkos-base."
    )

