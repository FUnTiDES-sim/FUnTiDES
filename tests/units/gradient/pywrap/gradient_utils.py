import kokkos
import numpy as np

default_memspace = kokkos.HostSpace
default_layout = kokkos.LayoutRight


def allocate_field(size, memspace=default_memspace, layout=default_layout):
    """Allocate a float32 Kokkos 1D array initialised to zero."""
    kk_arr = kokkos.array(
        [size], dtype=kokkos.float32, space=memspace, layout=layout
    )
    arr = np.array(kk_arr, copy=False)
    arr[:] = 0.0
    return kk_arr, arr
