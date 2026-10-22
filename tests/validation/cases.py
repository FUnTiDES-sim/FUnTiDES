"""
Validation cases: solution at one receiver vs the closed-form homogeneous
full-space solution.

SEM: the 24 configurations of the former funtides-sem driver
(tests/validation/validate_solution.cc of the application):
  2 meshes (cartesian, ucartesian) x 3 orders x 2 model locations x 2 physics.
DG: acoustic only, model on elements only (the DG kernels read the model per
element), 2 meshes x 3 orders.
DG-SEM and p-adaptive DG: same restrictions, 2 meshes x 3 orders (DG-SEM) or
3 (pMin, pMax) pairs. The medium is still homogeneous: these cases check that
the wave crosses the interface between the two sub-domains unchanged.

The parameters not set by that driver are the funtides-sem defaults: 2000 m cube,
source at (1010, 1010, 1010), receiver at (1310, 1310, 1310), Ricker of order 2,
f0 = 5 Hz, tpeak = 0.2 s, no absorbing layer. The first boundary reflection
reaches the receiver after ~1.1 s, past the 0.7 s window.
"""

import itertools
import math
from dataclasses import dataclass

DOMAIN_SIZE = 2000.0
SRC_POSITION = (1010.0, 1010.0, 1010.0)
RCV_POSITION = (1310.0, 1310.0, 1310.0)
# The coupled solvers split the mesh at INTERFACE_Z, their default: DG (or pMin) below, SEM (or pMax)
# above. Their source lies 150 m below it and their receiver 165 m above (off the element faces of
# every mesh, where the DG fields are discontinuous), at about the same distance
# as in the other cases (530 m); the first boundary reflection still arrives after 1.3 s.
INTERFACE_Z = 1000.0
SRC_POSITION_COUPLED_SOLVER = (1010.0, 1010.0, 850.0)
RCV_POSITION_COUPLED_SOLVER = (1310.0, 1310.0, 1165.0)
F0 = 5.0
RICKER_ORDER = 2
TPEAK = 0.2
DT = 0.001
# The explicit DG scheme is stable only for dt of the order of h / (vp * order^2) (SIPG penalty): its time
# step is DG_CFL times that, floored to a multiple of 1e-5 s. The SEM step of 1 ms (about 0.5 of that
# scale at order 3 with 27 m elements) blows up. VP is the Vp of the default homogeneous model.
# Measured on the cartesian mesh: order 1 passes at 0.2 and blows up at 0.3, orders 2 and 3 pass at 0.3.
# The coupled solvers use the value of their highest order.
DG_CFL = {1: 0.15, 2: 0.2, 3: 0.2}
VP = 1500.0
TIME_MAX = 0.7
# The element count per direction is chosen so that every order keeps the same number of nodes per
# direction. The DG fields are stored per element, without shared nodes, hence the coarser mesh.
NODES_PER_DIRECTION = {"sem": 400, "dg": 300, "dgsem": 300, "dgpadaptive": 300}

METHODS = ("sem", "dg", "dgsem", "dgpadaptive")
PHYSICS = ("acoustic", "elastic")
MESHES = ("cartesian", "ucartesian")
ORDERS = (1, 2, 3)
# (pMin, pMax) pairs of the p-adaptive DG solver, which needs pMin < pMax.
PADAPTIVE_ORDERS = ((1, 2), (1, 3), (2, 3))

# Relative L2 thresholds per (physics, order), about 2 to 3 times the measured error. The mesh type and
# the model location do not change the error. Orders 2 and 3 reach the same ~1e-3 floor, set by the
# second-order time scheme: the DG runs give the same 1e-3 at order 2 with dt = 1 ms, and 5.5e-4 at
# dt = 0.33 ms, the time part scaling as dt^2.
TOLERANCE = {
    ("acoustic", 1): 3e-2,
    ("acoustic", 2): 3e-3,
    ("acoustic", 3): 3e-3,
    ("elastic", 1): 7e-2,
    ("elastic", 2): 3e-3,
    ("elastic", 3): 3e-3,
}
# DG and DG-SEM thresholds, about 2 to 4 times the measured error: DG 2.5e-2, 6.5e-4 and 4.3e-4 at
# orders 1 to 3, DG-SEM 2.6e-2, 4.6e-4 and 3.9e-4 (the interface adds no visible error).
DG_TOLERANCE = {1: 5e-2, 2: 2e-3, 3: 2e-3}
# p-adaptive thresholds per (pMin, pMax), about twice the measured error. The mesh is sized for pMax,
# so the pMin side is under-resolved and dominates the error: 3.3e-2 for (1, 2) with 20 m elements,
# 5.9e-2 for (1, 3) with 27 m elements, 1.3e-3 for (2, 3).
PADAPTIVE_TOLERANCE = {(1, 2): 7e-2, (1, 3): 1.2e-1, (2, 3): 3e-3}


@dataclass(frozen=True)
class ValidationCase:
    method: str
    physics: str
    mesh: str
    order: int
    on_nodes: bool
    # Lower order of the p-adaptive DG solver (pMin), 0 for the other methods. order is then pMax.
    order_min: int = 0

    @property
    def is_elastic(self):
        return self.physics == "elastic"

    @property
    def is_coupled(self):
        return self.method in ("dgsem", "dgpadaptive")

    @property
    def source(self):
        return SRC_POSITION_COUPLED_SOLVER if self.is_coupled else SRC_POSITION

    @property
    def receiver(self):
        return RCV_POSITION_COUPLED_SOLVER if self.is_coupled else RCV_POSITION

    def in_lower_domain(self, point):
        """Whether the point lies in the DG (DG-SEM) or pMin (p-adaptive) sub-domain."""
        return point[2] < INTERFACE_Z

    def basis_order(self, point):
        """Order of the basis interpolating the field at the point."""
        if self.method == "dgpadaptive" and self.in_lower_domain(point):
            return self.order_min
        return self.order

    @property
    def elements_per_direction(self):
        return NODES_PER_DIRECTION[self.method] // (self.order + 1)

    @property
    def dt(self):
        if self.method == "sem":
            return DT
        h = DOMAIN_SIZE / self.elements_per_direction
        return round(math.floor(DG_CFL[self.order] * h / (VP * self.order**2) / 1e-5) * 1e-5, 5)

    @property
    def n_samples(self):
        return int(round(TIME_MAX / self.dt))

    @property
    def components(self):
        # Column names of the trace files written by the analytical scripts.
        return ("ux", "uy", "uz") if self.is_elastic else ("pressure_at_receiver",)

    @property
    def tolerance(self):
        if self.method == "dgpadaptive":
            return PADAPTIVE_TOLERANCE[(self.order_min, self.order)]
        if self.method in ("dg", "dgsem"):
            return DG_TOLERANCE[self.order]
        return TOLERANCE[(self.physics, self.order)]

    @property
    def name(self):
        suffix = "_model_on_nodes" if self.on_nodes else ""
        prefix = "" if self.method == "sem" else f"{self.method}_"
        order = f"order{self.order_min}to{self.order}" if self.order_min else f"order{self.order}"
        return f"{prefix}{self.physics}_{order}_{self.mesh}{suffix}"


def all_cases():
    sem = [
        ValidationCase("sem", physics, mesh, order, on_nodes)
        for physics, mesh, order, on_nodes in itertools.product(PHYSICS, MESHES, ORDERS, (False, True))
    ]
    dg = [
        ValidationCase(method, "acoustic", mesh, order, False)
        for method, mesh, order in itertools.product(("dg", "dgsem"), MESHES, ORDERS)
    ]
    padaptive = [
        ValidationCase("dgpadaptive", "acoustic", mesh, order_max, False, order_min)
        for mesh, (order_min, order_max) in itertools.product(MESHES, PADAPTIVE_ORDERS)
    ]
    return sem + dg + padaptive
