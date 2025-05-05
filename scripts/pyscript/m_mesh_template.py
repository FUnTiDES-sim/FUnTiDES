"""Module defining a Propagator_mesh object describing a propagation mesh"""
from pysabl.m_grid import Grid
from pysabl.m_sep_handler import Sep_handler
from pysabl.m_array_3d import Array_3d
from pysabl.m_grid import Grid
from pysabl.m_parallelism_context import Parallelism
from pysabl.m_dd_cart_grid_ctx import Dd_cart_grid
from pydw.m_earth_model import Earth_model
from pydw.m_propagator_param import Propagator_param
from pydw.m_propagator_mesh import Propagator_mesh

class Mesh_template(Propagator_mesh):
    def __init__(self, handle=None):
        """
        Constructor
        """
        pass


    def get_parameters(self, param):
        """
        Get user parameters from command line and/or Propagator_param structure

        Parameters
        ----------
        param : Propagator_param
            Context which contains common wave equation parameters [in]
        """
        pass


    def init(self, cluster, grid, earth_model):
        """
        Initialize propagation mesh.

        Parameters
        ----------
        cluster : T_Parallelism
            Parallelism object containing MPI communicator [in]
        grid : Grid
            Cartesian grid that propagation mesh needs to cover [in]
        earth_model : Earth_model
            Earth model object containing all physical models (Vp, density, etc.) [in]
        """
        pass


    def project_from_grid(self, grid, dd_grid, name, array_in, array_out):
        """
        Project an array allocated on a cartesian grid to the propagation mesh.

        Parameters
        ----------
        grid : Grid
            Grid context of cartesian grid [in]
        dd_grid : Dd_cart_grid
            Domain Decomposition context of cartesian grid [in]
        name : str
            Wavefield name (if any) [in]
        array_in : float array
            Input real array (on cartesian grid) [in]
        array_out : Array_3d
            Output real array (on propagation mesh)
        """
        pass


    def project_to_grid(self, grid, dd_grid, name, array_in, array_out):
        """
        Project an array allocated on propagation mesh to a cartesian grid.

        Parameters
        ----------
        grid : Grid
            Grid context of cartesian grid [in]
        dd_grid : Dd_cart_grid
            Domain Decomposition context of cartesian grid [in]
        name : str
            Waefield name (if any) [in]
        array_in : Array_3d
            Input real array (on propagation mesh) [in]
        array_out : float array
            Output real array (on cartesian grid) [out]
        """
        pass


    def info(self):
        """
        Print information about propagation mesh
        """
        pass


    def free(self):
        """
        Free propagator mesh
        """
        pass

    def __del__(self):
        """
        Destructor
        """
        pass
