"""Module defining a Projected_data_tw object handling projection of seismic traces on propagation mesh"""
from pydw.m_projected_data_tw import Projected_data_tw
from pydw.m_field_tw import Field_tw
from pydw.m_propagator_mesh import Propagator_mesh
from pysabl.m_dd_cart_grid_ctx import DD_cart_grid
from pydw.m_seismic_trace import Seismic_trace


class Projected_data_fd(Projected_data_tw):
    def __init__(self):
        """
        Constructor
        """
        pass


    def init(self, mesh, field, traces, ldim2, lx, ly, lz, o_tag=None, o_ireg=None, \
        o_dirichlet=None):
        """
        Initialize of Projected_data_tw object. Called before propagation loop.

        Parameters
        ----------
        mesh : Propagator_mesh
            propagation mesh [in]
        field : Field_tw
            field object that will be propagated [in]
        traces : Seismic_trace
            seismic traces to project on propagation mesh [in]
        ldim2 : bool
            flag for 2D grid [in]
        lx : int32
            Ghost size in x direction, if > 0 takes account of injection in ghost zones [in]
        ly : int32
            Ghost size in y direction, if > 0 takes account of injection in ghost zones [in]
        lz : int32
            Ghost size in z direction, if > 0 takes account of injection in ghost zones [in]
        o_tag : str
            optional header message to print in projected_data_info [in]
        o_ireg : int32
            optional half length of sinc window [in]
        o_dirichlet : bool
            optional flag to inject data as a dirichlet condition instead of a source term \
                [in]
        """
        pass


    def set_dirichlet(self, dirichlet):
        """
        Specify if data needs to be injected with Dirichlet boundary condition.
        This method might be called before or after propagation loop.
        """
        pass


    def get(self, n_comp, data, mesh, field, name):
        """
        Extract wavefield value at seismic traces positions. Called during propagation.

        Parameters
        ----------
        n_comp : int32
            number of seismic traces components
        data : float32 2D numpy "F" array of size [n_comp, self.n_data_impact_domain]
            output array where we want to put extracted data [in,out]
        mesh : Propagator_mesh
            propagation mesh context [in]
        field : Field_tw
            wavefields context [in,out]
        name : str array
            name of wavefields to extract [in]
        """
        pass


    def info(self):
        """
        Print information about projected_data object
        """
        pass


    def free(self):
        """
        Free Projected_data object
        """
        pass


    def __del__(self):
        """
        Destructor
        """
        pass
