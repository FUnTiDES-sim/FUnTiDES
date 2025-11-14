class StructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 100
        self.domain_size = 2000
        self.hx = self.domain_size / self.ex
        self.hy = self.domain_size / self.ey
        self.hz = self.domain_size / self.ez
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz


class UnstructData:
    def __init__(self, order):
        self.ex = self.ey = self.ez = 100
        self.lx = self.ly = self.lz = 2000
        self.order = order
        self.nx = self.ex * self.order + 1
        self.ny = self.ey * self.order + 1
        self.nz = self.ez * self.order + 1
        self.n_dof = self.nx * self.ny * self.nz
