"""
Validation cases: SEM solution at one receiver vs the closed-form homogeneous
full-space solution.

Same 24 configurations as the former funtides-sem driver
(tests/validation/validate_solution.cc of the application):
  2 meshes (cartesian, ucartesian) x 3 orders x 2 model locations x 2 physics.

The parameters not set by that driver are the funtides-sem defaults: 2000 m cube,
source at (1010, 1010, 1010), receiver at (1310, 1310, 1310), Ricker of order 2,
f0 = 5 Hz, tpeak = 0.2 s, no absorbing layer. The first boundary reflection
reaches the receiver after ~1.1 s, past the 0.7 s window.
"""

import itertools
from dataclasses import dataclass

DOMAIN_SIZE = 2000.0
SOURCE = (1010.0, 1010.0, 1010.0)
RECEIVER = (1310.0, 1310.0, 1310.0)
F0 = 5.0
RICKER_ORDER = 2
TPEAK = 0.2
DT = 0.001
TIME_MAX = 0.7
# The element count per direction is chosen so that every order keeps ~400 nodes per direction.
NODES_PER_DIRECTION = 400

PHYSICS = ("acoustic", "elastic")
MESHES = ("cartesian", "ucartesian")
ORDERS = (1, 2, 3)

# Relative L2 thresholds per (physics, order), about 2 to 3 times the measured error. The mesh type and
# the model location do not change the error. Orders 2 and 3 reach the same ~1e-3 floor, likely set by the
# second-order time scheme rather than by the space discretization.
TOLERANCE = {
    ("acoustic", 1): 3e-2,
    ("acoustic", 2): 3e-3,
    ("acoustic", 3): 3e-3,
    ("elastic", 1): 7e-2,
    ("elastic", 2): 3e-3,
    ("elastic", 3): 3e-3,
}


@dataclass(frozen=True)
class ValidationCase:
    physics: str
    mesh: str
    order: int
    on_nodes: bool

    @property
    def is_elastic(self):
        return self.physics == "elastic"

    @property
    def elements_per_direction(self):
        return NODES_PER_DIRECTION // (self.order + 1)

    @property
    def n_samples(self):
        return int(round(TIME_MAX / DT))

    @property
    def components(self):
        # Column names of the trace files written by the analytical scripts.
        return ("ux", "uy", "uz") if self.is_elastic else ("pressure_at_receiver",)

    @property
    def tolerance(self):
        return TOLERANCE[(self.physics, self.order)]

    @property
    def name(self):
        suffix = "_model_on_nodes" if self.on_nodes else ""
        return f"{self.physics}_order{self.order}_{self.mesh}{suffix}"


def all_cases():
    return [
        ValidationCase(physics, mesh, order, on_nodes)
        for physics, mesh, order, on_nodes in itertools.product(PHYSICS, MESHES, ORDERS, (False, True))
    ]
