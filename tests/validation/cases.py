"""
Validation cases: solution at one receiver vs the closed-form homogeneous
full-space solution.

SEM: the 24 configurations of the former funtides-sem driver
(tests/validation/validate_solution.cc of the application):
  2 meshes (cartesian, ucartesian) x 3 orders x 2 model locations x 2 physics.
DG: acoustic only, model on elements only (the DG kernels read the model per
element), 2 meshes x 3 orders.

The parameters not set by that driver are the funtides-sem defaults: 2000 m cube,
source at (1010, 1010, 1010), receiver at (1310, 1310, 1310), Ricker of order 2,
f0 = 5 Hz, tpeak = 0.2 s, no absorbing layer. The first boundary reflection
reaches the receiver after ~1.1 s, past the 0.7 s window.
"""

import itertools
import math
from dataclasses import dataclass

DOMAIN_SIZE = 2000.0
SOURCE = (1010.0, 1010.0, 1010.0)
RECEIVER = (1310.0, 1310.0, 1310.0)
F0 = 5.0
RICKER_ORDER = 2
TPEAK = 0.2
DT = 0.001
# The explicit DG scheme is stable only for dt of the order of h / (vp * order^2) (SIPG penalty): its time
# step is DG_CFL times that, floored to a multiple of 1e-5 s. The SEM step of 1 ms (about 0.5 of that
# scale at order 3 with 27 m elements) blows up. VP is the Vp of the default homogeneous model.
# Measured on the cartesian mesh: order 1 passes at 0.2 and blows up at 0.3, orders 2 and 3 pass at 0.3.
DG_CFL = {1: 0.15, 2: 0.2, 3: 0.2}
VP = 1500.0
TIME_MAX = 0.7
# The element count per direction is chosen so that every order keeps the same number of nodes per
# direction. The DG fields are stored per element, without shared nodes, hence the coarser mesh.
NODES_PER_DIRECTION = {"sem": 400, "dg": 300}

METHODS = ("sem", "dg")
PHYSICS = ("acoustic", "elastic")
MESHES = ("cartesian", "ucartesian")
ORDERS = (1, 2, 3)

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
# DG thresholds, about 2 to 3 times the measured error: 2.4e-2 at order 1, 6.5e-4 at order 2 and
# 5.3e-4 at order 3.
DG_TOLERANCE = {1: 5e-2, 2: 2e-3, 3: 2e-3}


@dataclass(frozen=True)
class ValidationCase:
    method: str
    physics: str
    mesh: str
    order: int
    on_nodes: bool

    @property
    def is_elastic(self):
        return self.physics == "elastic"

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
        return DG_TOLERANCE[self.order] if self.method == "dg" else TOLERANCE[(self.physics, self.order)]

    @property
    def name(self):
        suffix = "_model_on_nodes" if self.on_nodes else ""
        prefix = "" if self.method == "sem" else f"{self.method}_"
        return f"{prefix}{self.physics}_order{self.order}_{self.mesh}{suffix}"


def all_cases():
    sem = [
        ValidationCase("sem", physics, mesh, order, on_nodes)
        for physics, mesh, order, on_nodes in itertools.product(PHYSICS, MESHES, ORDERS, (False, True))
    ]
    dg = [ValidationCase("dg", "acoustic", mesh, order, False) for mesh, order in itertools.product(MESHES, ORDERS)]
    return sem + dg
