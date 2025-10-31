import kokkos
import numpy as np

default_memspace = kokkos.HostSpace
default_layout = kokkos.LayoutRight


def allocate_pressure(n_dof, memspace=default_memspace, layout=default_layout):
    kk_pnGlobal = kokkos.array(
        [n_dof, 2], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pnGlobal = np.array(kk_pnGlobal, copy=False)
    pnGlobal[:] = 0.0

    return kk_pnGlobal, pnGlobal

def allocate_displacementx(n_dof, memspace=default_memspace, layout=default_layout):
    kk_uxnGlobal = kokkos.array(
        [n_dof, 2], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uxnGlobal = np.array(kk_uxnGlobal, copy=False)
    uxnGlobal[:] = 0.0

    return kk_uxnGlobal, uxnGlobal

def allocate_displacementy(n_dof, memspace=default_memspace, layout=default_layout):
    kk_uynGlobal = kokkos.array(
        [n_dof, 2], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uynGlobal = np.array(kk_uynGlobal, copy=False)
    uynGlobal[:] = 0.0

    return kk_uynGlobal, uynGlobal

def allocate_displacementz(n_dof, memspace=default_memspace, layout=default_layout):
    kk_uznGlobal = kokkos.array(
        [n_dof, 2], dtype=kokkos.float32, space=memspace, layout=layout
    )
    uznGlobal = np.array(kk_uznGlobal, copy=False)
    uznGlobal[:] = 0.0

    return kk_uznGlobal, uznGlobal

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
    RHSElement[0] = ex / 2 + ey / 2 * ex + ez / 2 * ey * ex  # one half of slice
    RHSElement[1] = ex / 3 + ey / 2 * ex + ez / 2 * ey * ex  # one third of slice
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
