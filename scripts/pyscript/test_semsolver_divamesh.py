"""
Test SEM solver with DIVA hexahedral mesh.
Uses pysabl, pydw and pydiva_mesh modules on top of regular SEM solver modules (kokkos, pysolver, pysem).
Outputs vtk files for visualization of the mesh and optionally GLL nodes.
Usage: mpirun -np 2 python test_semsolver_divamesh.py par=<mesher_parameter_file>
"""

# BASE
import numpy as np
import time
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# SABL
import pysabl.parser as parser
from pysabl.m_database import Database
from pysabl.m_grid import Grid
from pysabl.m_parallel_timer import timing, timer_print
from pysabl.m_parallelism_context import Parallelism

# DIVA-WAVE-EQUATION
from pydw.m_earth_model import Earth_model
from pydw.m_propagator_param import Propagator_param

# DIVA-MESH
from pydiva_mesh.core.m_hexa_mesh import Hexa_mesh
from pydiva_mesh.core.m_hexa_model import Mesh_model
from pydiva_mesh.core.m_hexa_vtk_export import Hexa_vtk_export
from pydiva_mesh.core.m_wrap_sem_proxy_unstruct import Wrap_sem_proxy_unstruct

# SEM-PROXY
import libpykokkos as kokkos
import pysem as Sem
import pysolver as Solver

ArrayInt = kokkos.KokkosView_int32_HostSpace_LayoutRight_2
ArrayReal = kokkos.KokkosView_float32_HostSpace_LayoutRight_2
VectorInt = kokkos.KokkosView_int32_HostSpace_LayoutRight_1
VectorReal = kokkos.KokkosView_float32_HostSpace_LayoutRight_1


class Mesher:
    """
    DIVA Mesher for SEM (without ghosting).
    """

    def __init__(self, paral: Parallelism, param: parser.Param):
        self.paral = paral
        self.prop_param = Propagator_param()
        self.db = Database()
        self.mesh = Hexa_mesh()
        self.grid = Grid()
        self.earth_model = Earth_model()
        self.mesh_model = Mesh_model()
        self.vtk_export = Hexa_vtk_export()
        self.wrap = Wrap_sem_proxy_unstruct()

        # Fix the grid size for a small case
        # Warning : make sure these are defined the same in the paramter file
        #  -> current limitation of DIVA
        self.ex = self.ey = self.ez = 100
        self.dx = self.dy = self.dz = 15
        self.ox = self.oy = self.oz = 0

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
        # Use the fixed grid size for the mesh and not the model grid
        grid_n = np.array([self.ex, self.ey, self.ez], dtype=int)
        grid_d = np.array([self.dx, self.dy, self.dz], dtype=float)
        grid_o = np.array([self.ox, self.oy, self.oz], dtype=float)
        self.grid.init_list(grid_n, grid_d, grid_o)
        self.mesh.init(self.paral, self.grid, self.earth_model)
        # Or use the grid from the earth model
        #self.mesh.init(self.paral, self.earth_model.vp_model, self.earth_model)
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
    def wrap_mesh(self):
        self.wrap.init(self.mesh, self.mesh_model)

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
            print(f"{'ex/ey/ez':<15} {[self.ex, self.ey, self.ez]}")
            print(f"{'ox/oy/oz':<15} {[self.ox, self.oy, self.oz]}")
            print(f"{'dx/dy/dz':<15} {[self.dx, self.dy, self.dz]}")
            print(f"{'-' * 65}")
            print(f"{'nglls':<15} {[ngllx, nglly, ngllz]}")
            print(f"{'min/max node':<15} {str(min_node):<10} / {str(max_node):<10}")
            print(f"{'min/max elem':<15} {str(min_elem):<10} / {str(max_elem):<10}")
            print(f"{'-' * 65}")
            print(f"{'Array':<15} {'Shape':<20} {'Min':<15} {'Max':<15}")
            print(f"{'-' * 65}")
            print(f"{'elem_to_node':<15} {str(elem_to_node.shape):<20} {int(elem_to_node.min()):<15} {int(elem_to_node.max()):<15}")
            print(f"{'coord_node X':<15} {str(coord_node[0].shape):<20} {int(coord_node[0].min()):<15} {int(coord_node[0].max()):<15}")
            print(f"{'coord_node Y':<15} {str(coord_node[1].shape):<20} {int(coord_node[1].min()):<15} {int(coord_node[1].max()):<15}")
            print(f"{'coord_node Z':<15} {str(coord_node[2].shape):<20} {int(coord_node[2].min()):<15} {int(coord_node[2].max()):<15}")
            print(f"{'rho_elem':<15} {str(rho_elem.shape):<20} {rho_elem.min():<15.6f} {rho_elem.max():<15.6f}")
            print(f"{'vp_elem':<15} {str(vp_elem.shape):<20} {vp_elem.min():<15.6f} {vp_elem.max():<15.6f}")
            print(f"{'vs_elem':<15} {str(vs_elem.shape):<20} {vs_elem.min():<15.6f} {vs_elem.max():<15.6f}")
            print(f"{'rho_node':<15} {str(rho_node.shape):<20} {rho_node.min():<15.6f} {rho_node.max():<15.6f}")
            print(f"{'vp_node':<15} {str(vp_node.shape):<20} {vp_node.min():<15.6f} {vp_node.max():<15.6f}")
            print(f"{'vs_node':<15} {str(vs_node.shape):<20} {vs_node.min():<15.6f} {vs_node.max():<15.6f}")
            print(f"{'-' * 65}", flush=True)

    @timing
    def export_gll(self):
        self.vtk_export.export_gll(self.paral, self.mesh, self.mesh_model, self.vtk_name)

    @timing
    def export_elements(self):
        self.vtk_export.export_elements(self.paral, self.mesh, self.mesh_model, self.vtk_name)

    @timing
    def free(self):
        self.mesh_model.free_model_elem()
        self.mesh_model.free_model_gll()
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
        self.wrap_mesh()
        if self.out_gll_vtk:
            self.export_gll()
        self.export_elements()


class SEM:
    """
    SEM solver using the DIVA hexahedral mesh.
    """

    def __init__(self, paral: Parallelism, wrap: Wrap_sem_proxy_unstruct):
        # from wrapper
        self.paral = paral
        self.ngllx = wrap.ngllx
        self.nglly = wrap.nglly
        self.ngllz = wrap.ngllz
        self.min_node = wrap.min_node
        self.max_node = wrap.max_node
        self.min_elem = wrap.min_elem
        self.max_elem = wrap.max_elem
        self.elem_to_node = wrap.elem_to_node
        self.coord_node = wrap.coord_node
        self.rho_elem = wrap.rho_elem
        self.vp_elem = wrap.vp_elem
        self.vs_elem = wrap.vs_elem
        self.rho_node = wrap.rho_node
        self.vp_node = wrap.vp_node
        self.vs_node = wrap.vs_node
        self.order = wrap.ngllx - 1
        self.ngll = wrap.ngllx * wrap.nglly * wrap.ngllz
        self.n_elem = wrap.max_elem - wrap.min_elem + 1
        self.n_node = wrap.max_node - wrap.min_node + 1

        # source parameters
        self.f0 = 5
        self.timeStep = 0.001
        self.nTimeSteps = 300
        self.numberOfRHS = 1

        # kokkos start
        kokkos.initialize()

    @timing
    def init_fe_arrays(self):
        print("Initializing FE arrays...", flush=True)
        #  id -> dof_global | TODO so far only interior nodes
        self.interiorNodes = np.linspace(0, self.n_node - 1, self.n_node, dtype=np.int32)
        # (elem,dof_local) -> dof_global
        self.nodesList = np.empty((self.n_elem, self.ngll), dtype=np.int32)
        # (elem,dof_local) -> coord
        self.nodesCoordsX = np.empty((self.n_elem, self.ngll), dtype=np.float32)
        self.nodesCoordsY = np.empty((self.n_elem, self.ngll), dtype=np.float32)
        self.nodesCoordsZ = np.empty((self.n_elem, self.ngll), dtype=np.float32)

        # TODO very inefficient
        for e in range(self.n_elem):
            for k in range(self.ngllz):
                for j in range(self.nglly):
                    for i in range(self.ngllx):
                        # dof_local_diva = i + j * self.ngllx + k * self.nglly * self.ngllx
                        dof_local_sem = i + k * self.ngllx + j * self.ngllz * self.ngllx
                        dof_global_diva = self.elem_to_node[i, j, k, e] - self.min_node
                        self.nodesList[e, dof_local_sem] = dof_global_diva

        # Use slicing to assign coordinates for all elements at once (assuming n_vertices == 8 for hexahedral elements)
        # The indices for the 8 corners of a hexahedron in (i, j, k) order:
        corner_indices = [
            (0, 0, 0),
            (-1, 0, 0),
            (0, 0, -1),
            (-1, 0, -1),
            (0, -1, 0),
            (-1, -1, 0),
            (0, -1, -1),
            (-1, -1, -1),
        ]
        for idx, (i, j, k) in enumerate(corner_indices):
            node_indices = self.elem_to_node[i, j, k, :] - self.min_node  # shape: (n_elem,)
            self.nodesCoordsX[:, idx] = self.coord_node[0, node_indices]
            self.nodesCoordsY[:, idx] = self.coord_node[1, node_indices]
            self.nodesCoordsZ[:, idx] = self.coord_node[2, node_indices]

        # wrap in kokkos views
        self.kk_interiorNodes = VectorInt(self.interiorNodes, (self.n_node,))
        self.kk_nodesList = ArrayInt(self.nodesList, (self.n_elem, self.ngll))
        self.kk_nodesCoordsX = ArrayReal(self.nodesCoordsX, (self.n_elem, self.ngll))
        self.kk_nodesCoordsY = ArrayReal(self.nodesCoordsY, (self.n_elem, self.ngll))
        self.kk_nodesCoordsZ = ArrayReal(self.nodesCoordsZ, (self.n_elem, self.ngll))
        self.kk_model = VectorReal(self.vp_elem, (self.n_elem,))

    @timing
    def init_info(self):
        print("Initializing SEMinfo...", flush=True)
        self.myInfo = Sem.SEMinfo()
        self.myInfo.numberOfNodes = self.n_node
        self.myInfo.numberOfElements = self.n_elem
        self.myInfo.numberOfPointsPerElement = self.ngll
        # TODO so far only interior nodes
        self.myInfo.numberOfInteriorNodes = self.n_node

    @timing
    def init_pressure(self):
        print("Initializing Pressure...", flush=True)
        self.pnGlobal = np.zeros((self.n_node, 2), dtype=np.float32)
        self.kk_pnGlobal = ArrayReal(self.pnGlobal, (self.n_node, 2))

    @timing
    def init_solver(self):
        print("Initializing Solver...", flush=True)
        self.solver = Solver.SEMsolver()
        self.solver.computeFEInitWithoutMesh(self.myInfo, self.kk_nodesList, self.kk_nodesCoordsX, self.kk_nodesCoordsY, self.kk_nodesCoordsZ, self.kk_interiorNodes, self.kk_model)

    @timing
    def init_source(self):
        print("Initializing Source...", flush=True)

        self.RHSElement = np.array([self.numberOfRHS], dtype=np.int32)

        # Find the element at minimum z and middle of x and y coordinates
        elem_centers_x = np.zeros(self.n_elem)
        elem_centers_y = np.zeros(self.n_elem)
        elem_centers_z = np.zeros(self.n_elem)
        for e in range(self.n_elem):
            node_indices = self.nodesList[e, :]
            x = self.coord_node[0, node_indices]
            y = self.coord_node[1, node_indices]
            z = self.coord_node[2, node_indices]
            elem_centers_x[e] = np.mean(x)
            elem_centers_y[e] = np.mean(y)
            elem_centers_z[e] = np.mean(z)
        mid_z = 0.5 * (np.min(elem_centers_z) + np.max(elem_centers_z))
        z_tol = 50
        mid_z_indices = np.where(np.abs(elem_centers_z - mid_z) < z_tol)[0]
        mid_x = 0.5 * (np.min(elem_centers_x) + np.max(elem_centers_x))
        mid_y = 0.5 * (np.min(elem_centers_y) + np.max(elem_centers_y))
        # Find the element among min_z_indices closest to (mid_x, mid_y)
        dists = (elem_centers_x[mid_z_indices] - mid_x) ** 2 + (elem_centers_y[mid_z_indices] - mid_y) ** 2
        min_idx = mid_z_indices[np.argmin(dists)]
        self.RHSElement[0] = min_idx
        self.kk_RHSElement = VectorInt(self.RHSElement, (self.numberOfRHS,))

        # Store coordinates of the RHS element's first node
        self.rhs_elem = self.RHSElement[0]
        self.rhs_node = self.nodesList[self.rhs_elem, 0]
        self.x_rhs = self.coord_node[0, self.rhs_node]
        self.y_rhs = self.coord_node[1, self.rhs_node]
        self.z_rhs = self.coord_node[2, self.rhs_node]

        # Find all nodes with the same X coordinate (within a tolerance) to use for plotting
        tol = 1e-6
        x_coords = self.coord_node[0, :]
        mask = np.abs(x_coords - self.x_rhs) < tol
        self.rhs_node_indices = np.where(mask)[0]
        self.nodes_rhs_slice_y = self.coord_node[1, self.rhs_node_indices]
        self.nodes_rhs_slice_z = self.coord_node[2, self.rhs_node_indices]

        # Apply source term
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

        print("Propagating...", flush=True)
        print(f"Simulation started at: {simulation_start}", flush=True)

        # Propagate one step forward in time
        i1 = 0
        i2 = 1
        self.maxval = 1

        for timeSample in range(self.nTimeSteps):
            print(f"   - time iteration {timeSample}/{self.nTimeSteps}", flush=True)
            iter_start = time.time()

            self.solver.computeOneStep(timeSample, self.order, self.ngll, i1, i2, self.myInfo, self.kk_RHSTerm, self.kk_pnGlobal, self.kk_RHSElement)

            iter_time = time.time() - iter_start
            iteration_times.append(iter_time)

            self.export_snapshot_png(timeSample, i1)

            if timeSample % 100 == 0:
                print(f"     sum pnGlobal[:, i1] = {np.sum(self.pnGlobal[:, i1])}")
                print(f"     RHSElement[0] = {self.RHSElement[0]}")
                self.maxval = np.max(np.abs(self.pnGlobal))/5
                print(f"     Pressure max val after iteration {timeSample}: {self.maxval}")
                print(f"     Average iteration time: {np.mean(iteration_times):.4f} seconds")
                elementSource = self.nodesList[self.RHSElement[0], 0]
                elapsed_time = time.time() - start_time
                print(f"     Time iteration {timeSample}/{self.nTimeSteps}")
                print(f"     Elapsed time: {elapsed_time:.2f} seconds")
                print(f"     Average iteration time: {np.mean(iteration_times):.4f} seconds")
                print(f"     Pressure={self.pnGlobal[elementSource,0]}")
                print(f"     Percentage of zeros in pnGlobal = {np.count_nonzero(self.pnGlobal[i1] == 0) / self.pnGlobal[i1].size * 100}", flush=True)

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
        print(f"Max iteration time: {np.max(iteration_times):.4f} seconds", flush=True)

    @timing
    def info(self):
        print(f"{'-' * 65}")
        print("SEM Solver Info")
        print(f"{'-' * 65}")
        print(f"{'RHS Element':<20} {self.rhs_elem}")
        print(f"{'RHS Node':<20} {self.rhs_node}")
        print(f"{'RHS Coordinates':<20} X={self.x_rhs:.6f}, Y={self.y_rhs:.6f}, Z={self.z_rhs:.6f}")
        print(f"{'f0':<20} {self.f0}")
        print(f"{'timeStep':<20} {self.timeStep}")
        print(f"{'nTimeSteps':<20} {self.nTimeSteps}")
        print(f"{'numberOfRHS':<20} {self.numberOfRHS}")
        print(f"{'Array':<20} {'Shape':<20} {'Min':<15} {'Max':<15}")
        print(f"{'-' * 65}")
        print(f"{'interiorNodes':<20} {str(self.interiorNodes.shape):<20} {self.interiorNodes.min():<15} {self.interiorNodes.max():<15}")
        print(f"{'nodesList':<20} {str(self.nodesList.shape):<20} {self.nodesList.min():<15} {self.nodesList.max():<15}")
        print(f"{'nodesCoordsX':<20} {str(self.nodesCoordsX.shape):<20} {self.nodesCoordsX.min():<15.6f} {self.nodesCoordsX.max():<15.6f}")
        print(f"{'nodesCoordsY':<20} {str(self.nodesCoordsY.shape):<20} {self.nodesCoordsY.min():<15.6f} {self.nodesCoordsY.max():<15.6f}")
        print(f"{'nodesCoordsZ':<20} {str(self.nodesCoordsZ.shape):<20} {self.nodesCoordsZ.min():<15.6f} {self.nodesCoordsZ.max():<15.6f}")
        print(f"{'RHSElement':<20} {str(self.RHSElement.shape):<20} {self.RHSElement.min():<15} {self.RHSElement.max():<15}")
        print(f"{'RHSTerm':<20} {str(self.RHSTerm.shape):<20} {self.RHSTerm.min():<15.6f} {self.RHSTerm.max():<15.6f}")
        print(f"{'model':<20} {str(self.vp_elem.shape):<20} {self.vp_elem.min():<15.6f} {self.vp_elem.max():<15.6f}")
        print(f"{'-' * 65}", flush=True)

    @timing
    def export_source_plot(self):
        """
        Export a plot of the source time function.
        """
        plt.figure(figsize=(8, 4))
        plt.plot(np.arange(self.nTimeSteps) * self.timeStep, self.RHSTerm[0, :], label="Source Time Function")
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude")
        plt.title("Source Time Function (RHSTerm[0, :])")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig("source_time_function.png", dpi=150)
        plt.close()
        print("Source time function plot saved as source_time_function.png", flush=True)

    @timing
    def export_snapshot_png(self, timeSample, i1):
        """
        Export a PNG snapshot of the pressure field for all nodes at the same X coordinate as the source element.
        """
        pressure_vals = self.pnGlobal[self.rhs_node_indices, i1].copy()
        # Normalize pressure values for better visualization
        if self.maxval  > 0:
            pressure_vals = pressure_vals / self.maxval

        # Interpolate using bilinear method (linear)
        grid_x, grid_y = np.mgrid[min(self.nodes_rhs_slice_y):max(self.nodes_rhs_slice_y):200j, min(self.nodes_rhs_slice_z):max(self.nodes_rhs_slice_z):200j]
        grid_z = griddata((self.nodes_rhs_slice_y, self.nodes_rhs_slice_z), pressure_vals, (grid_x, grid_y), method='nearest')

        plt.figure(figsize=(8, 6))
        plt.imshow(grid_z.T, origin='lower', cmap='viridis', interpolation='nearest', vmin=-0.1, vmax=0.1)
        plt.colorbar(label='Pressure')
        plt.xlabel('Y')
        plt.ylabel('Z')
        plt.title(f'Pressure field at X={self.x_rhs:.3f}, timeSample={timeSample}')
        plt.tight_layout()
        file_name = f"sem_snapshot_{timeSample}.png"
        plt.savefig(file_name, dpi=150)
        plt.close()
        print(f"Snapshot PNG file written as {file_name}", flush=True)

    @timing
    def free(self):
        del self.solver
        del self.myInfo
        kokkos.finalize()

    def run(self):
        self.init_info()
        self.init_fe_arrays()
        self.init_pressure()
        self.init_source()
        self.export_source_plot()
        self.info()
        self.init_solver()
        self.propagate()


def main():
    # Initialize MPI parallelism
    paral = Parallelism()
    paral.initialize()

    # Initialize parametrization
    param = parser.Param(
        name="Solver SEM with DIVA hexahedral mesh",
        purpose="Test SEM solver with DIVA hexahedral mesh (requires DIVA par file)",
        parallelism="SEM solver uses Kokkos for parallelism, Mesher uses MPI and OMP/OACC",
        author="G. Fuss",
    )

    # Mesh
    mesher = Mesher(paral, param)
    mesher.run()

    # SEM
    sem = SEM(paral, mesher.wrap)
    sem.run()

    # Free
    sem.free()
    mesher.free()

    timer_print()


if __name__ == "__main__":
    main()
