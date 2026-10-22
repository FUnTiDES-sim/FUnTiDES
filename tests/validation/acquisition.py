"""
Point source and receiver on a cartesian mesh: element lookup and Lagrange
interpolation weights (Python port of SourceAndReceiverUtils::ComputeRHSWeights).

Elements are axis-aligned boxes, so the reference coordinates follow from the
two opposite corners of the element.
"""

import math

import numpy as np

""" Gauss-Lobatto-Legendre points for orders 1, 2, and 3. """
GLL_POINTS = {
    1: np.array([-1.0, 1.0]),
    2: np.array([-1.0, 0.0, 1.0]),
    3: np.array([-1.0, -1.0 / math.sqrt(5.0), 1.0 / math.sqrt(5.0), 1.0]),
}


def lagrange_1d(order, xi):
    """Values at xi of the order+1 Lagrange polynomials built on the GLL points."""
    nodes = GLL_POINTS[order]
    values = np.ones(order + 1)
    for a in range(order + 1):
        for b in range(order + 1):
            if b != a:
                values[a] *= (xi - nodes[b]) / (nodes[a] - nodes[b])
    return values


def element_index(point, n_elements, domain_size):
    """Index i + j*ex + k*ex*ey of the element holding the point, same numbering as the cartesian builders."""
    ijk = []
    for x, n in zip(point, n_elements):
        if not 0.0 <= x < domain_size:
            raise ValueError(f"point {point} lies outside the [0, {domain_size}) domain")
        ijk.append(int(math.floor(x * n / domain_size)))
    ex, ey, _ = n_elements
    return ijk[0] + ijk[1] * ex + ijk[2] * ex * ey


def element_nodes(model, element, order):
    """Global node indices of an element, local index i + j*(order+1) + k*(order+1)^2."""
    n = order + 1
    return np.array(
        [model.global_node_index(element, i, j, k) for k in range(n) for j in range(n) for i in range(n)],
        dtype=np.int64,
    )


def interpolation_weights(model, element, order, point, basis_order=None):
    """
    Basis function values at the point, local index i + j*(n) + k*(n)^2 with n = basis_order + 1.

    order is the order of the mesh, used to locate the element corners; basis_order, which
    defaults to it, is the order of the interpolating basis (lower in the pMin domain of the
    p-adaptive DG solver).
    """
    basis_order = order if basis_order is None else basis_order
    corners = (model.global_node_index(element, 0, 0, 0), model.global_node_index(element, order, order, order))
    lo, hi = (np.array([model.node_coord(int(n), d) for d in range(3)], dtype=np.float64) for n in corners)
    xi = -1.0 + 2.0 * (np.asarray(point, dtype=np.float64) - lo) / (hi - lo)
    if np.any(np.abs(xi) > 1.0 + 1e-6):
        raise ValueError(f"point {point} is not inside element {element} (reference coordinates {xi})")
    lx, ly, lz = (lagrange_1d(basis_order, x) for x in xi)
    return np.einsum("k,j,i->kji", lz, ly, lx).ravel()
