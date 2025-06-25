"""
Mesher for SEM.
This mesher is used to create a hexahedric mesh from a SEP with projected vp.
Uses pysabl, pydw and pydiva_mesh modules and runs in parallel using MPI.
Outputs vtk files for visualization of the mesh and optionally GLL nodes.
Usage: mpirun -np 2 python mesher.py par=<mesher_parameter_file>
"""

import numpy as np
import time
from datetime import datetime

# SEM
import pysabl.parser as parser
from pysabl.m_parallel_timer import timing, timer_print
from pysabl.m_parallelism_context import Parallelism
from pydw.m_earth_model import Earth_model
from pysabl.m_database import Database
from pydw.m_propagator_param import Propagator_param
from pydiva_mesh.core.m_hexa_mesh import Hexa_mesh
from pydiva_mesh.core.m_hexa_model import Mesh_model
from pydiva_mesh.core.m_hexa_vtk_export import Hexa_vtk_export
from pydiva_mesh.core.m_wrap_sem_proxy_unstruct import Wrap_sem_proxy_unstruct
from pydiva_mesh.core.m_hexa_ghost_unstruct import Mesh_ghost_unstruct

# PROXY
import libpykokkos as kokkos
import pysolver as Solver
import pysem as Sem

ArrayInt = kokkos.KokkosView_int32_HostSpace_LayoutRight_2
ArrayReal = kokkos.KokkosView_float32_HostSpace_LayoutRight_2
VectorInt = kokkos.KokkosView_int32_HostSpace_LayoutRight_1
VectorReal = kokkos.KokkosView_float32_HostSpace_LayoutRight_1


class Mesher:
    def __init__(self, paral: Parallelism, param: parser.Param):
        self.paral = paral
        self.prop_param = Propagator_param()
        self.db = Database()
        self.mesh = Hexa_mesh()
        self.earth_model = Earth_model()
        self.mesh_model = Mesh_model()
        self.vtk_export = Hexa_vtk_export()
        self.wrap = Wrap_sem_proxy_unstruct()
        self.ghost = Mesh_ghost_unstruct()

        self.outdir = param.frompar(
            "outdir",
            default="/tmp",
            help="Directory where to store temporary propagation results (if any)",
        )
        self.vtk_name = param.frompar(
            "vtk_out_file",
            default="hexamesh",
            help="Basename of output vtk files. If gll nodes requested, suffix -gll will be added.",
        )
        self.out_gll_vtk = param.frompar(
            "out_gll_vtk",
            default=False,
            help="Whether gll nodes should be exported into a vtk file.",
        )
        self.out_ghost_vtk = param.frompar(
            "out_ghost_vtk",
            default=False,
            help="Whether ghosts should be exported into a vtk file.",
        )

    @timing
    def init_propagator_param(self):
        self.prop_param.init(self.outdir)

    @timing
    def init_mesh_from_param(self):
        self.mesh.get_parameters(self.prop_param)

    @timing
    def init_earth_model(self):
        self.db.init("identity")
        self.earth_model.init(self.db)
        self.earth_model.info(self.paral.get_rank() == 0)

    @timing
    def read_earth_model(self):
        self.earth_model.read(self.paral)

    @timing
    def build_mesh(self):
        self.mesh.init(self.paral, self.earth_model.vp_model, self.earth_model)
        self.mesh.info()

    @timing
    def load_model_on_mesh(self):
        self.mesh_model.init_model_gll(self.mesh, self.earth_model)
        self.mesh_model.load_model_gll(self.mesh, self.earth_model)
        self.mesh_model.info_model_gll(self.paral.get_rank() == 0)
        self.mesh_model.init_model_elem(self.mesh, self.earth_model)
        self.mesh_model.extrapolate_model_elem(self.mesh, self.earth_model)
        self.mesh_model.info_model_elem(self.paral.get_rank() == 0)

    @timing
    def com_ghost(self):
        self.ghost.init(self.mesh, o_model_elem=True, o_model_node=True)
        self.ghost.compute(self.mesh, self.mesh_model)
        self.ghost.info(self.paral.get_rank() == 0)

    @timing
    def wrap_mesh(self):
        self.wrap.init(self.mesh, self.mesh_model, self.ghost)

        ngllx = self.wrap.ngllx
        nglly = self.wrap.nglly
        ngllz = self.wrap.ngllz
        min_node = self.wrap.min_node
        max_node = self.wrap.max_node
        min_elem = self.wrap.min_elem
        max_elem = self.wrap.max_elem
        elem_to_node = self.wrap.elem_to_node
        coord_node = self.wrap.coord_node
        rho_elem = self.wrap.rho_elem
        vp_elem = self.wrap.vp_elem
        vs_elem = self.wrap.vs_elem
        rho_node = self.wrap.rho_node
        vp_node = self.wrap.vp_node
        vs_node = self.wrap.vs_node

        if self.paral.get_rank() == 0:
            print(f"{'-' * 65}")
            print("Hexahedral mesh wrap")
            print(f"{'-' * 65}")
            print(f"{'nglls':<15} {[ngllx, nglly, ngllz]}")
            print(f"{'min/max node':<15} {str(min_node):<10} / {str(max_node):<10}")
            print(f"{'min/max elem':<15} {str(min_elem):<10} / {str(max_elem):<10}")
            print(f"{'-' * 65}")
            print(f"{'Array':<15} {'Shape':<20} {'Min':<15} {'Max':<15}")
            print(f"{'-' * 65}")
            print(
                f"{'elem_to_node':<15} {str(elem_to_node.shape):<20} {int(elem_to_node.min()):<15} {int(elem_to_node.max()):<15}"
            )
            print(
                f"{'coord_node':<15} {str(coord_node.shape):<20} {int(coord_node.min()):<15} {int(coord_node.max()):<15}"
            )
            print(
                f"{'rho_elem':<15} {str(rho_elem.shape):<20} {rho_elem.min():<15.6f} {rho_elem.max():<15.6f}"
            )
            print(
                f"{'vp_elem':<15} {str(vp_elem.shape):<20} {vp_elem.min():<15.6f} {vp_elem.max():<15.6f}"
            )
            print(
                f"{'vs_elem':<15} {str(vs_elem.shape):<20} {vs_elem.min():<15.6f} {vs_elem.max():<15.6f}"
            )
            print(
                f"{'rho_node':<15} {str(rho_node.shape):<20} {rho_node.min():<15.6f} {rho_node.max():<15.6f}"
            )
            print(
                f"{'vp_node':<15} {str(vp_node.shape):<20} {vp_node.min():<15.6f} {vp_node.max():<15.6f}"
            )
            print(
                f"{'vs_node':<15} {str(vs_node.shape):<20} {vs_node.min():<15.6f} {vs_node.max():<15.6f}"
            )

    @timing
    def export_ghost(self):
        import vtk

        n_ghost_node = self.wrap.n_ghost_node
        ghost_coord = self.wrap.ghost_coord_node
        ghost_coord = ghost_coord.reshape((3, n_ghost_node), order="F")

        # Write ghost nodes as a VTP (vtkPolyData) file of points
        points = vtk.vtkPoints()
        for i in range(n_ghost_node):
            points.InsertNextPoint(ghost_coord[0, i], ghost_coord[1, i], ghost_coord[2, i])

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)

        # Add a vertex cell for each point so they are visible in a preview
        vertices = vtk.vtkCellArray()
        for i in range(n_ghost_node):
            vertices.InsertNextCell(1)
            vertices.InsertCellPoint(i)
        polydata.SetVerts(vertices)

        writer = vtk.vtkXMLPolyDataWriter()
        filename = f"ghost_points{self.paral.get_rank()}.vtp"
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        writer.Write()

    @timing
    def export_gll(self):
        self.vtk_export.export_gll(
            self.paral, self.mesh, self.mesh_model, self.vtk_name
        )

    @timing
    def export_elements(self):
        self.vtk_export.export_elements(
            self.paral, self.mesh, self.mesh_model, self.vtk_name
        )

    @timing
    def free(self):
        self.mesh_model.free_model_elem()
        self.mesh_model.free_model_gll()
        self.ghost.free()
        self.wrap.free()
        self.mesh.free()
        self.earth_model.free()

    def run(self):
        self.init_propagator_param()
        self.init_mesh_from_param()
        self.init_earth_model()
        self.read_earth_model()
        self.build_mesh()
        self.load_model_on_mesh()
        self.com_ghost()
        self.wrap_mesh()
        if self.out_ghost_vtk:
            self.export_ghost()
        if self.out_gll_vtk:
            self.export_gll()
        self.export_elements()
        self.free()
        timer_print()


class SEM:

    def __init__(self, paral: Parallelism, wrap: Wrap_sem_proxy_unstruct):
        # from wrapper
        self.paral = paral
        self.ngllx = self.wrap.ngllx
        self.nglly = self.wrap.nglly
        self.ngllz = self.wrap.ngllz
        self.min_node = self.wrap.min_node
        self.max_node = self.wrap.max_node
        self.min_elem = self.wrap.min_elem
        self.max_elem = self.wrap.max_elem
        self.elem_to_node = self.wrap.elem_to_node
        self.coord_node = self.wrap.coord_node
        self.rho_elem = self.wrap.rho_elem
        self.vp_elem = self.wrap.vp_elem
        self.vs_elem = self.wrap.vs_elem
        self.rho_node = self.wrap.rho_node
        self.vp_node = self.wrap.vp_node
        self.vs_node = self.wrap.vs_node

        # for solver
        self.order = self.wrap.ngllx - 1
        self.ngll = self.wrap.ngllx * self.wrap.nglly * self.wrap.ngllz
        self.n_elem = self.wrap.max_elem - self.wrap.max_elem + 1
        self.n_node = self.wrap.max_node - self.wrap.min_node + 1
        kokkos.initialize()

    @timing
    def init_fe_arrays(self):
        print("Initializing FE arrays...")
        self.interiorNodes = np.linspace(0, self.n_node - 1, self.n_node, dtype=int)  # TODO so far only interior nodes

        # TODO very inneficient
        # (elem,dof_local) -> dof_global
        self.nodesList = np.zeros((self.n_elem, self.ngll), dtype=int)
        # (elem,dof_local) -> coord
        self.nodesCoordsX = np.zeros((self.n_elem, self.ngll), dtype=np.float32)
        self.nodesCoordsY = np.zeros((self.n_elem, self.ngll), dtype=np.float32)
        self.nodesCoordsZ = np.zeros((self.n_elem, self.ngll), dtype=np.float32)
        for e in range(self.n_elem):
            for k in range(self.ngllz):
                for j in range(self.nglly):
                    for i in range(self.ngllx):
                        dof_local = i + j * self.ngllx + k * self.ngllx * self.nglly
                        dof_global = self.elem_to_node[i, j, k, e] - self.min_node
                        self.nodesList[e, dof_local] = dof_global
                        self.nodesCoordsX[e, dof_local] = self.coord_node[0, dof_global]
                        self.nodesCoordsY[e, dof_local] = self.coord_node[1, dof_global]
                        self.nodesCoordsZ[e, dof_local] = self.coord_node[2, dof_global]

        # wrap in kokkos views
        self.kk_interiorNodes = VectorInt(self.interiorNodes, (self.n_node,))
        self.kk_nodesList = ArrayInt(self.nodesList, (self.n_elem, self.ngll))
        self.kk_nodesCoordsX = ArrayReal(self.nodesCoordsX, (self.n_elem, self.ngll))
        self.kk_nodesCoordsY = ArrayReal(self.nodesCoordsY, (self.n_elem, self.ngll))
        self.kk_nodesCoordsZ = ArrayReal(self.nodesCoordsZ, (self.n_elem, self.ngll))
        self.kk_model = VectorReal(self.vp_elem, (self.n_elem,))

    @timing
    def init_info(self):
        print("Initializing SEMinfo...")
        self.myInfo = Sem.SEMinfo()
        self.myInfo.numberOfNodes = self.n_node
        self.myInfo.numberOfElements = self.n_elem
        self.myInfo.numberOfPointsPerElement = self.ngll
        self.myInfo.numberOfInteriorNodes = self.n_node   # TODO so far only interior nodes

    @timing
    def init_pressure(self):
        print("Initializing Pressure...")
        self.pnGlobal = np.zeros((self.n_node, 2), dtype=np.float32)
        self.kk_pnGlobal = ArrayReal(self.pnGlobal, (self.n_node, 2))

    @timing
    def init_solver(self):
        print("Initializing Solver...")
        self.solver = Solver.SEMsolver()
        self.solver.computeFEInitWithoutMesh(self.myInfo, self.kk_nodesList, self.kk_nodesCoordsX, self.kk_nodesCoordsY, self.kk_nodesCoordsZ, self.kk_interiorNodes, self.kk_model)

    @timing
    def inti_source(self):
        print("Initializing Source...")
        self.f0 = 5
        self.timeStep = 0.001
        self.nTimeSteps = 300
        self.numberOfRHS = 1

        self.RHSElement = np.array([self.numberOfRHS],dtype=int)
        self.RHSElement[0] = self.max_elem - self.min_elem + 1 # TODO just take mid element for now
        self.kk_RHSElement = VectorInt(self.RHSElement, (self.numberOfRHS,))
        print("   - RHS element number ", self.RHSElement[0])

        self.RHSTerm = np.zeros((self.numberOfRHS, self.nTimeSteps), dtype=np.float32)
        for i in range(self.nTimeSteps):
            self.RHSTerm[0, i] = self.sourceTerm(i*self.timeStep, self.f0)
        self.kk_RHSTerm = ArrayReal(self.RHSTerm, (self.numberOfRHS, self.nTimeSteps))

    def sourceTerm(self, time_n, f0):
        o_tpeak = 1.0/f0
        pulse = 0.0
        if time_n <= -0.9*o_tpeak or time_n >= 2.9*o_tpeak:
            return pulse

        pi = 3.14157
        lam = (f0*pi)*(f0*pi)
        pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*np.exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak))
        return pulse

    @timing
    def propagate(self):
        # Add timing variables
        start_time = time.time()
        simulation_start = datetime.now()
        iteration_times = []

        print(f"Simulation started at: {simulation_start}")

        # propagate one step forward in time
        i1 = 0
        i2 = 1
        maxval = 1

        for timeSample in range(self.nTimeSteps):
            iter_start = time.time()

            if timeSample % 100 == 0:
                print("sum pnGlobal[:, i1]", np.sum(self.pnGlobal[:, i1]))

            self.solver.computeOneStep(timeSample, self.order, self.ngll, i1, i2, self.myInfo, self.kk_RHSTerm, self.kk_pnGlobal, self.kk_RHSElement)

            iter_time = time.time() - iter_start
            iteration_times.append(iter_time)

            if timeSample % 100 == 0:
                print("sum pnGlobal[:, i1]", np.sum(self.pnGlobal[:, i1]))
                print(f"RHSElement[0] = {self.RHSElement[0]}")
                maxval = np.max(np.abs(self.pnGlobal))/5
                print(f"Pressure max val after iteration {timeSample}: {maxval}")
                print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
                elementSource = self.nodesList[self.RHSElement[0], 0]
                elapsed_time = time.time() - start_time
                print(f"Time iteration {timeSample}/{self.nTimeSteps}")
                print(f"Elapsed time: {elapsed_time:.2f} seconds")
                print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
                print(f"Pressure={self.pnGlobal[elementSource,0]}")

            print(f"Percentage of zeros in pnGlobal = {np.count_nonzero(self.pnGlobal[i1] == 0) / self.pnGlobal[i1].size * 100}")

            tmp = i1
            i1 = i2
            i2 = tmp

        # Print final timing statistics
        end_time = time.time()
        simulation_end = datetime.now()
        total_time = end_time - start_time

        print("\nSimulation Statistics:")
        print(f"Start time: {simulation_start}")
        print(f"End time: {simulation_end}")
        print(f"Total runtime: {total_time:.2f} seconds")
        print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
        print(f"Min iteration time: {np.min(iteration_times):.4f} seconds")
        print(f"Max iteration time: {np.max(iteration_times):.4f} seconds")

    @timing
    def free(self):
        del self.solver
        del self.myInfo
        kokkos.finalize()

    @timing
    def run(self):
        self.init_info()
        self.init_fe_arrays()
        self.init_pressure()
        self.init_solver()
        self.inti_source()
        self.propagate()
        self.free()


def main():
    # Initialize MPI parallelism
    paral = Parallelism()
    paral.initialize()

    # Initialize parametrization
    param = parser.Param(
        name="SEM GLL mesher",
        purpose="Create Hexahedric mesh from SEP with projected vp",
        parallelism="MPI parallelism on vertical axis",
        author="G. Fuss",
    )

    # Mesh
    mesher = Mesher(paral, param)
    mesher.run()

    # SEM
    sem = SEM(paral, mesher.wrap)
    sem.run()


if __name__ == "__main__":
    main()
