import kokkos
import numpy as np

# Dynamically determine the correct memory space and layout based on how Kokkos was compiled.
# GPU builds (CUDA) expect LayoutLeft and UVM space to share memory with NumPy without segfaulting.
# CPU builds expect LayoutRight and HostSpace.
if hasattr(kokkos, "CudaUVMSpace"):
    default_memspace = kokkos.CudaUVMSpace
    default_layout = kokkos.LayoutLeft
elif hasattr(kokkos, "CudaSpace"):
    default_memspace = kokkos.CudaSpace
    default_layout = kokkos.LayoutLeft
elif hasattr(kokkos, "HIPManagedSpace"):
    default_memspace = kokkos.HIPManagedSpace
    default_layout = kokkos.LayoutLeft
else:
    default_memspace = kokkos.HostSpace
    default_layout = kokkos.LayoutRight


def allocate_pressure(n_dof, memspace=default_memspace, layout=default_layout):
    kk_pnGlobalPrev = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pnGlobalPrev = np.array(kk_pnGlobalPrev, copy=False)
    pnGlobalPrev[:] = 0.0

    kk_pnGlobalCurr = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pnGlobalCurr = np.array(kk_pnGlobalCurr, copy=False)
    pnGlobalCurr[:] = 0.0

    return kk_pnGlobalPrev, pnGlobalPrev, kk_pnGlobalCurr, pnGlobalCurr


def allocate_displacementx(n_dof, memspace=default_memspace, layout=default_layout):
    kk_uxnGlobalPrev = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uxnGlobalPrev = np.array(kk_uxnGlobalPrev, copy=False)
    uxnGlobalPrev[:] = 0.0

    kk_uxnGlobalCurr = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uxnGlobalCurr = np.array(kk_uxnGlobalCurr, copy=False)
    uxnGlobalCurr[:] = 0.0

    return kk_uxnGlobalPrev, uxnGlobalPrev, kk_uxnGlobalCurr, uxnGlobalCurr


def allocate_displacementy(n_dof, memspace=default_memspace, layout=default_layout):
    kk_uynGlobalPrev = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uynGlobalPrev = np.array(kk_uynGlobalPrev, copy=False)
    uynGlobalPrev[:] = 0.0

    kk_uynGlobalCurr = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uynGlobalCurr = np.array(kk_uynGlobalCurr, copy=False)
    uynGlobalCurr[:] = 0.0

    return kk_uynGlobalPrev, uynGlobalPrev, kk_uynGlobalCurr, uynGlobalCurr


def allocate_displacementz(n_dof, memspace=default_memspace, layout=default_layout):
    kk_uznGlobalPrev = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uznGlobalPrev = np.array(kk_uznGlobalPrev, copy=False)
    uznGlobalPrev[:] = 0.0

    kk_uznGlobalCurr = kokkos.array(
        [n_dof], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uznGlobalCurr = np.array(kk_uznGlobalCurr, copy=False)
    uznGlobalCurr[:] = 0.0

    return kk_uznGlobalPrev, uznGlobalPrev, kk_uznGlobalCurr, uznGlobalCurr


def allocate_rhs_term(
    n_rhs, n_time_steps, dt, f0, memspace=default_memspace, layout=default_layout
):
    kk_RHSTerm = kokkos.array(
        [n_rhs, n_time_steps], dtype=kokkos.float32, space=memspace, layout=layout
    )
    RHSTerm = np.array(kk_RHSTerm, copy=False)
    for i in range(n_time_steps):
        RHSTerm[0, i] = source_term(i * dt, f0)
        RHSTerm[1, i] = source_term(i * dt, f0)
    return kk_RHSTerm, RHSTerm


def allocate_rhs_weight(n_rhs, model, memspace=default_memspace, layout=default_layout):
    nb_points = model.get_number_of_points_per_element()
    kk_RHSWeights = kokkos.array(
        [n_rhs, nb_points],
        dtype=kokkos.float32,
        space=memspace,
        layout=layout,
    )
    RHSWeights = np.array(kk_RHSWeights, copy=False)
    for i in range(n_rhs):
        for j in range(model.get_number_of_points_per_element()):
            RHSWeights[i, j] = 1 / model.get_number_of_points_per_element()
    return kk_RHSWeights, RHSWeights


def allocate_rhs_element(
    n_rhs, ex, ey, ez, memspace=default_memspace, layout=default_layout
):
    kk_RHSElement = kokkos.array(
        [n_rhs], dtype=kokkos.int32, space=memspace, layout=layout
    )
    RHSElement = np.array(kk_RHSElement, copy=False)
    # Wrapped in int() to prevent float-to-int NumPy casting errors
    RHSElement[0] = int(ex / 2 + ey / 2 * ex + ez / 2 * ey * ex)  # one half of slice
    RHSElement[1] = int(ex / 3 + ey / 2 * ex + ez / 2 * ey * ex)  # one third of slice
    return kk_RHSElement, RHSElement


def source_term(time_n, f0):
    o_tpeak = 1.0 / f0
    pulse = 0.0
    if time_n <= -0.9 * o_tpeak or time_n >= 2.9 * o_tpeak:
        return pulse

    pi = 3.14157
    lam = (f0 * pi) * (f0 * pi)
    pulse = (
        2.0
        * lam
        * (2.0 * lam * (time_n - o_tpeak) * (time_n - o_tpeak) - 1.0)
        * np.exp(-lam * (time_n - o_tpeak) * (time_n - o_tpeak))
    )
    return pulse
