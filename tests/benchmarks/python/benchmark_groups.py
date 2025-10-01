from enum import Enum, auto


class BenchmarkGroup(Enum):
    COMPUTE_FE_INIT = auto()
    COMPUTE_ONE_STEP = auto()
