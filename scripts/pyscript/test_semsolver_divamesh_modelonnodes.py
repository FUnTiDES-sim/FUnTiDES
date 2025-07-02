"""
Mesher for SEM.
This mesher is used to create a hexahedric mesh from a SEP with projected vp.
Uses pysabl, pydw and pydiva_mesh modules and runs in parallel using MPI.
Outputs vtk files for visualization of the mesh and optionally GLL nodes.
Usage: mpirun -np 2 python mesher.py par=<mesher_parameter_file>
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

import time
from datetime import datetime

# SABL
import pysabl.parser as parser
from pysabl.m_parallel_timer import timing, timer_print
from pysabl.m_parallelism_context import Parallelism
from pysabl.m_database import Database
from pysabl.m_grid import Grid

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
import pysolver as Solver
import pysem as Sem

Array4DInt = kokkos.KokkosView_int64_HostSpace_LayoutRight_4
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
        self.grid = Grid()
        self.earth_model = Earth_model()
        self.mesh_model = Mesh_model()
        self.vtk_export = Hexa_vtk_export()
        self.wrap = Wrap_sem_proxy_unstruct()

        # Fix the grid size for a small case
        # Warning : make sure these that the values used here are the same in parameter file 
        # -> Current limitation of DIVA
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
        grid_n = np.array([self.ex+1, self.ey+1, self.ez+1], dtype=int)
        grid_d = np.array([self.dx, self.dy, self.dz], dtype=float)
        grid_o = np.array([self.ox, self.oy, self.oz], dtype=float)
        self.grid.init_list(grid_n, grid_d, grid_o)
        self.mesh.init(self.paral, self.grid, self.earth_model)
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
        #rho_elem = self.wrap.rho_elem
        #vp_elem = self.wrap.vp_elem
        #vs_elem = self.wrap.vs_elem
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
            print(f"{'elem_to_node':<15} {str(elem_to_node.shape):<20} {int(elem_to_node.min()):<15} {int(elem_to_node.max()):<15}")
            print(f"{'coord_node X':<15} {str(coord_node[0].shape):<20} {int(coord_node[0].min()):<15} {int(coord_node[0].max()):<15}")
            print(f"{'coord_node Y':<15} {str(coord_node[1].shape):<20} {int(coord_node[1].min()):<15} {int(coord_node[1].max()):<15}")
            print(f"{'coord_node Z':<15} {str(coord_node[2].shape):<20} {int(coord_node[2].min()):<15} {int(coord_node[2].max()):<15}")
            #print(f"{'rho_elem':<15} {str(rho_elem.shape):<20} {rho_elem.min():<15.6f} {rho_elem.max():<15.6f}")
            #print(f"{'vp_elem':<15} {str(vp_elem.shape):<20} {vp_elem.min():<15.6f} {vp_elem.max():<15.6f}")
            #print(f"{'vs_elem':<15} {str(vs_elem.shape):<20} {vs_elem.min():<15.6f} {vs_elem.max():<15.6f}")
            print(f"{'rho_node':<15} {str(rho_node.shape):<20} {rho_node.min():<15.6f} {rho_node.max():<15.6f}")
            print(f"{'vp_node':<15} {str(vp_node.shape):<20} {vp_node.min():<15.6f} {vp_node.max():<15.6f}")
            print(f"{'vs_node':<15} {str(vs_node.shape):<20} {vs_node.min():<15.6f} {vs_node.max():<15.6f}")
            print(f"{'-' * 65}", flush=True)

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
        self.elem_to_node = wrap.elem_to_node.T - self.min_node # make it zero-based
        self.coord_node = wrap.coord_node.T
        self.elem_to_node = np.asfortranarray(self.elem_to_node)
        self.coord_node = np.asfortranarray(self.coord_node)
        print("DBG: type of elem_to_node", self.elem_to_node.dtype, flush=True)
        print(self.elem_to_node.flags)  # Should show F_CONTIGUOUS = True
        print(wrap.elem_to_node.flags)  # Should show F_CONTIGUOUS = True
        print("DBG: elem_to_node shape", self.elem_to_node.shape, flush=True)
        print("DBG: type of coord_node", self.coord_node.dtype, flush=True)
        print(self.coord_node.flags)  # Should show F_CONTIGUOUS = True
        print(wrap.coord_node.flags)  # Should show F_CONTIGUOUS = True
        print("DBG: coord_node shape", self.coord_node.shape, flush=True)
        self.rho_elem = wrap.rho_elem
        self.vp_elem = wrap.vp_elem
        self.vs_elem = wrap.vs_elem
        self.rho_node = wrap.rho_node
        self.vp_node = wrap.vp_node
        self.vs_node = wrap.vs_node

        # for solver
        self.order = wrap.ngllx - 1 # ngllx = nglly = ngllz
        self.ngll = wrap.ngllx * wrap.nglly * wrap.ngllz
        self.n_elem = wrap.max_elem - wrap.min_elem + 1
        self.n_node = wrap.max_node - wrap.min_node + 1
        kokkos.initialize()

        # for plotting only
        self.ex = self.ey = self.ez = 100 # to create the snapshots in python

    @timing
    def init_fe_arrays(self):
        print("Initializing FE arrays...", flush=True)
        self.interiorNodes = np.linspace(0, self.n_node - 1, self.n_node, dtype=np.int32)  # TODO so far only interior nodes

        # wrap in kokkos views
        self.kk_interiorNodes = VectorInt(self.interiorNodes[:], (self.n_node,))
        print(f"Python code. interiorNodes original : {self.interiorNodes[10]}")
        print(f"Python code. interiorNodes kokkos : {self.kk_interiorNodes[10]}")
        self.kk_elem_to_node = Array4DInt(self.elem_to_node.T, (self.ngllx, self.nglly, self.ngllz, self.n_elem))
        print(f"Python elem_to_node.shape before reshape: {self.elem_to_node.shape}", flush=True)
        self.elem_to_node = self.elem_to_node.T # Important (but is it really?) to ensure they are the same in python and kokkos 
        print(f"Python elem_to_node.shape after reshape: {self.elem_to_node.shape}", flush=True)
        self.kk_coord_node = ArrayReal(self.coord_node.T, (3, self.n_node))
        self.coord_node = self.coord_node.T # Important to ensure they are the same in python and kokkos

        corner_indices = [
            (0, 0, 0),
            (self.order, 0, 0),
            (0, 0, self.order),
            (self.order, 0, self.order),
            (0, self.order, 0),
            (self.order, self.order, 0),
            (0, self.order, self.order),
            (self.order, self.order, self.order),
        ]
        for idx, (i, j, k) in enumerate(corner_indices):

            print(f"Corner {idx}: e2n[{i},{j},{k},0] = {self.elem_to_node[i, j, k, 0]}")
            print(f"Kokkos Corner {idx}: e2n[{i},{j},{k},0] = {self.kk_elem_to_node[i, j, k, 0]}")

        
            print("Original coord x, y, z:", self.coord_node[:,self.elem_to_node[i, j, k, 0]])
            print("Kokkos coord:", np.array(self.kk_coord_node)[:,self.kk_elem_to_node[i, j, k, 0]])

        self.kk_vpmodel = VectorReal(self.vp_node[:], (self.n_node,))
        self.kk_rhomodel = VectorReal(self.rho_node[:], (self.n_node,))

        print(self.vp_node)
        print(self.rho_node)

    @timing
    def init_info(self):
        print("Initializing SEMinfo...", flush=True)
        self.myInfo = Sem.SEMinfo()
        self.myInfo.numberOfNodes = self.n_node
        self.myInfo.numberOfElements = self.n_elem
        self.myInfo.numberOfPointsPerElement = self.ngll
        self.myInfo.numberOfInteriorNodes = self.n_node   # TODO so far only interior nodes

    @timing
    def init_pressure(self):
        print("Initializing Pressure...", flush=True)
        self.pnGlobal = np.zeros((self.n_node, 2), dtype=np.float32)
        self.kk_pnGlobal = ArrayReal(self.pnGlobal, (self.n_node, 2))

    @timing
    def init_solver(self):
        print("Initializing Solver...", flush=True)
        self.solver = Solver.SEMsolver()
        self.solver.allocateSolverDIVA(self.myInfo)
        self.solver.FEInitDIVA(self.myInfo,self.kk_elem_to_node, self.kk_coord_node,self.kk_interiorNodes,
                                self.kk_rhomodel, self.kk_vpmodel)
        self.solver.set_rhomodel(self.kk_rhomodel)
        self.solver.set_vpmodel(self.kk_vpmodel)
      
    @timing
    def inti_source(self):
        print("Initializing Source...", flush=True)
        self.f0 = 5
        self.timeStep = 0.001
        self.nTimeSteps = 300
        self.numberOfRHS = 1
        self.maxAmplitude = 1000

        self.RHSElement = np.array([self.numberOfRHS], dtype=np.int32)

        # Find the element at minimum z and middle of x and y coordinates
        elem_centers_x = np.zeros(self.n_elem)
        elem_centers_y = np.zeros(self.n_elem)
        elem_centers_z = np.zeros(self.n_elem)
        for e in range(self.n_elem):
            i = [0, self.ngllx - 1]
            j = [0, self.nglly - 1]
            k = [0, self.ngllz - 1]
            node_indices = [self.elem_to_node[ii, jj, kk, e] for ii in i for jj in j for kk in k]
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
        print("   - RHS element number ", self.RHSElement[0], flush=True)

        print("elementtonode orig, kk",self.elem_to_node[0,0,0,self.RHSElement[0]], self.kk_elem_to_node[0,0,0,self.RHSElement[0]], flush=True)
        tmpcoordX = self.kk_coord_node[0, self.kk_elem_to_node[0, 0, 0, self.RHSElement[0]]]
        tmpcoordY = self.kk_coord_node[1, self.kk_elem_to_node[0, 0, 0, self.RHSElement[0]]]
        tmpcoordZ = self.kk_coord_node[2, self.kk_elem_to_node[0, 0, 0, self.RHSElement[0]]]
        print(f"   - Source element coordinates: X={tmpcoordX}, Y={tmpcoordY}, Z={tmpcoordZ}", flush=True)

        # Store coordinates of the RHS element's first node
        self.x_rhs = tmpcoordX
        self.y_rhs = tmpcoordY
        self.z_rhs = tmpcoordZ

        # Find all nodes with the same X coordinate (within a tolerance) to use for plotting
        tol = 1e-6
        x_coords = self.coord_node[0, :]
        mask = np.abs(x_coords - self.x_rhs) < tol
        self.rhs_node_indices = np.where(mask)[0]
        self.nodes_rhs_slice_y = self.coord_node[1, self.rhs_node_indices]
        self.nodes_rhs_slice_z = self.coord_node[2, self.rhs_node_indices]

        self.RHSTerm = np.zeros((self.numberOfRHS, self.nTimeSteps), dtype=np.float32)
        for i in range(self.nTimeSteps):
            self.RHSTerm[0, i] = self.sourceTerm(i*self.timeStep, self.f0)
        self.RHSTerm = self.RHSTerm/np.max(np.abs(self.RHSTerm))
        self.RHSTerm *= self.maxAmplitude
        self.kk_RHSTerm = ArrayReal(self.RHSTerm, (self.numberOfRHS, self.nTimeSteps))

    def sourceTerm(self, time_n, f0):
        o_tpeak = 1.2/f0 # golden rule (it just works this way)
        # This is close to 2.0*np.sqrt(np.pi) / (3*f0), where a cut-off frequency of 3*f0 is assumed
        pulse = 0.0

        # Standard truncation for Ricker wavelet
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
       
        # propagate one step forward in time
        i1 = 0
        i2 = 1
        maxval = 1

        for timeSample in range(self.nTimeSteps):
            print(f"   - time iteration {timeSample}/{self.nTimeSteps}", flush=True)
            print(f"i1 = {i1}, i2 = {i2}", flush=True)
            iter_start = time.time()

            self.solver.computeOneStepDIVA(timeSample, self.order, self.ngll, i1, i2, self.myInfo, self.kk_RHSTerm, self.kk_pnGlobal, self.kk_RHSElement)

            iter_time = time.time() - iter_start
            iteration_times.append(iter_time)

            if timeSample % 10 == 0:
                print(f"     sum pnGlobal[:, i1] = {np.sum(self.pnGlobal[:, i1])}")
                print(f"     Percentage of zeros in pnGlobal = {np.count_nonzero(self.pnGlobal[:,i1] == 0) / self.pnGlobal[:,i1].size * 100}", flush=True)
                print(f"     RHSElement[0] = {self.RHSElement[0]}")
                maxval = np.max(np.abs(self.pnGlobal[:,i1]))              
                self.getSnapshot(maxval, timeSample, i1)
                print(f"     Pressure max val after iteration {timeSample}: {maxval}")
                print(f"     Average iteration time: {np.mean(iteration_times):.4f} seconds")
                elapsed_time = time.time() - start_time
                print(f"     Time iteration {timeSample}/{self.nTimeSteps}")
                print(f"     Elapsed time: {elapsed_time:.2f} seconds")
                print(f"     Average iteration time: {np.mean(iteration_times):.4f} seconds", flush=True)

            i1, i2 = i2, i1

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
        print(f"{'Array':<20} {'Shape':<20} {'Min':<15} {'Max':<15}")
        print(f"{'-' * 65}")
        print(f"{'interiorNodes':<20} {str(self.interiorNodes.shape):<20} {self.interiorNodes.min():<15} {self.interiorNodes.max():<15}")
        print(f"{'RHSElement':<20} {str(self.RHSElement.shape):<20} {self.RHSElement.min():<15} {self.RHSElement.max():<15}")
        print(f"{'RHSTerm':<20} {str(self.RHSTerm.shape):<20} {self.RHSTerm.min():<15.6f} {self.RHSTerm.max():<15.6f}")
        print(f"{'-' * 65}", flush=True)

    @timing
    def getSnapshot(self,maxval, timeSample, i1):
        
        pressure_vals = self.pnGlobal[self.rhs_node_indices, i1].copy()
        # Normalize pressure values for better visualization
        if maxval  > 0:
            pressure_vals = pressure_vals / maxval

        pressure_vals = np.reshape(pressure_vals, (self.ey*self.order+1, self.ez*self.order+1))
       
        maxval = np.max(np.abs(pressure_vals))/50
        pressure_vals = np.clip(pressure_vals, -maxval, maxval)

        plt.figure(figsize=(8, 6))
        plt.imshow(pressure_vals.T , origin='lower', cmap='viridis', interpolation='nearest')
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
    def free(self):
        del self.solver
        del self.myInfo
        kokkos.finalize()

    def run(self):
        self.init_info()
        self.init_fe_arrays()
        self.init_pressure()
        self.inti_source()
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

    # Free
    sem.free()
    mesher.free()

    timer_print()


if __name__ == "__main__":
    main()
