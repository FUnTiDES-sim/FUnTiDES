"""Module defining a Propagator_tw object handling two-way propagation"""
from pysabl.m_mem_size import Mem_size
from pysabl.m_grid import Grid
from pysabl.m_list_array_3d import List_array_3d
from pysabl.m_array_3d import Array_3d
from pysabl.m_parallelism_context import Parallelism
from pydw.m_projected_data_tw import Projected_data_tw
from pydw.m_earth_model import Earth_model
from pydw.m_propagator_param import Propagator_param
from pydw.m_field_tw import Field_tw
from pydw.m_propagator_mesh import Propagator_mesh
from pydw.m_propagator_tw import Propagator_tw
from pydw.m_seismic_trace import Seismic_trace

class Propagator_ac_cpml(Propagator_tw):
    def __init__(self):
        """
        Constructor
        """
        pass

    def get_parameters(self, param, o_flag_born=None, o_grad_list=None, \
        o_grad_type=None):
        """
        Get user parameters from command line and/or Propagator_param structure

        Parameters
        ----------
        param : Propagator_param
            Context which contains common wave equation parameters [in]
        o_flag_born : bool
            True if propagator will be used for Born modeling
        o_grad_list : str array
            List of names of physical parameters wrt which we will compute gradients
        o_grad_type : int32
            Flag to specify which type of gradient will be computed (adjoint, true amplitude, etc.)
        """
        pass


    def precompute(self, cluster, global_grid, earth_model, param):
        """
        Pre-compute data common to all the shots(e.g. global mesh, stencils, etc.).
        This method is called before looping over the shots.

        Parameters
        ----------
        cluster : Parallelism
            Parallelism context [in]
        global_grid : Grid
            Cartesian grid containing the whole area to process [in]
        earth_model : Earth_model
            Earth model object containing all physical models [in,out]
        param : Propagator_param
            Propagation parameters [in]
        """
        pass


    def get_dt(self, mesh, earth_model):
        """
        Function returns the maximum time step allowed by the numerical scheme

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        earth_model : Earth_model
            Earth model object containing all physical models [in]

        Returns
        -------
        propagator_ac_cpml_get_dt : float32
            Propagation time step
        """
        pass


    def set_dt(self, dt):
        """
        Set propagator numerical time step (rounded by the application).

        Parameters
        ----------
        dt : float32
            Time step [in]
        """
        pass


    def allocate(self, mesh, o_mem=None):
        """
        Allocate the required memory for the propagator.
        If `o_mem` is provided, only update memory counters.


        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        o_mem : Mem_size
            Optional t_mem_size, if provided update memory estimation instead of allocating \
                arrays [in,out]
        """
        pass


    def init(self, mesh, earth_model):
        """
        Initialize the propagator : load earth models, initialize damping, etc.
        This routine is called after `allocate`.

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        earth_model : Earth_model
            Earth model object containing all physical models [in,out]
        """
        pass


    def scale_data(self, prop_dir, mesh, data_xyz, data_t):
        """
        Apply any time-independant scaling to seismic traces that will be injected.
        Example: rho*vp² for 2nd order acoustic wave equation.
        This method is called before propagation loop.

        Parameters
        ----------
        prop_dir : int32
            Direction and type of propagation (forward/backward in time, adjoint, etc.) [in]
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        data_xyz : Projected_data_tw
            Spatial projection of seismic traces that will be injected on propagation mesh [in,out]
        data_t : Seismic_trace
            Seismic traces to scale [in,out]
        """
        pass


    def postprocess_data(self, data_xyz, trace):
        """
        Apply any time-independant scaling to data that has been recorded.
        This method is called after propagation loop.

        Parameters
        ----------
        data_xyz : Projected_data_tw
            Spatial projection on propagation mesh of seismic traces that have been recorded [in,out]
        trace : Seismic_trace
            Seismic traces to update [in,out]
        """
        pass


    def one_step(self, t, prop_dir, mesh, field, data_xyz, data_t):
        """
        Propagation: move one time step wavefields, including source injection, communications, boundary conditions.

        Parameters
        ----------
        t : float32
            Time of current iteration [in]
        prop_dir : int32
            Direction and type of propagation (forward/backward in time, adjoint, etc.) [in]
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        field : Field_tw
            Context of the field to propagate [in,out]
        data_xyz : Projected_data_tw
            Spatial projection on propagator mesh of seismic traces to inject [in,out]
        data_t : Seismic_trace
            Seismic traces to inject [in,out]
        """
        pass


    def one_step_linear(self, t, prop_dir, prop_dir_linear, mesh, field, \
        field_linear, data_xyz, data_t, n_pert, list_pert, image, flag_first_call):
        """
        Linear propagation (e.g. born modeling)

        Parameters
        ----------
        t : float32
            Time of current iteration [in]
        prop_dir : int32
            Direction and type of propagation of background wavefields (forward/backward in time, adjoint, etc.) [in]
        prop_dir_linear : int32
            Direction and type of propagation of scatter wavefields (forward/backward in time, adjoint, etc.) [in]
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        field : Field_tw
            Context of the source field to propagate [in,out]
        field_linear : Field_tw
            Context of the linear field to propagate [in,out]
        data_xyz : Projected_data_tw
            Spatial projection on propagation mesh of seismic traces to inject [in,out]
        data_t : Seismic_trace
            Seismic traces to inject [in,out]
        n_pert : int32
            Number of physical parameters perturbations [in]
        list_pert : str array
            List of strings describing which physical parameters perturbations are provided \
                [in]
        image : List_array_3d
            Linear sources [in]
        flag_first_call : bool
        """
        pass


    def one_step_zo(self, t, prop_dir, mesh, field, data_xyz, data_t, flag_t0, \
        image_zo):
        """
        Zero-offset migration/modeling (exploding reflector)

        Parameters
        ----------
        prop : T_Propagator_Ac_Cpml
            Context of the propagator [in,out]
        t : float32
            Time of current iteration [in]
        prop_dir : int32
            Direction and type of propagation (forward/backward in time, adjoint, etc.) [in]
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        field : Field_tw
            Wavefields to propagate [in,out]
        data_xyz : Projected_data_tw
            Spatial projection of data to inject [in,out]
        data_t : Seismic_trace
            Temporal projection of data to inject [in,out]
        flag_t0 : bool
            True at t=0(inject or construct zero-offset image) [in]
        image_zo : Array_3d
            Zero-Offset image [in,out]
        """
        pass


    def compute_gradient(self, mesh, field_fwd, field_bwd, n_grad, list_grad, \
        gradient, grad_type, flag_last_call):
        """
        Compute gradient of wavefield wrt physical parameters.
        This method is called during backward propagation of adjoint source.

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        field_fwd : Field_tw
            Forward wavefields [in,out]
        field_bwd : Field_tw
            Backward wavefields [in,out]
        n_grad : int32
            Number of physical parameters gradients to compute [in]
        list_grad : str array
            Lit of strings describing which physical parameters gradients need to be \
                computed [in]
        gradient : List_array_3d
            Gradient arrays to update [in,out]
        grad_type : int32
            Type of gradient computation (basic, true amplitude, adjoint, etc.) [in]
        flag_last_call : bool
            True if it is last call to gradient computation [in]
        """
        pass


    def compute_born_gradient(self, mesh, field_fwd, field_fwd_linear, field_bwd, \
        field_bwd_linear, n_pert, list_pert, pert, n_grad, list_grad, gradient, \
        grad_type, flag_last_call):
        """
        Compute gradient of Born wavefield wrt physical (background) parameters

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh context [in,out]
        field_fwd : Field_tw
            Forward wavefields [in,out]
        field_fwd_linear : Field_tw
            Forward Born wavefields [in,out]
        field_bwd : Field_tw
            Backward wavefields [in,out]
        field_bwd_linear : Field_tw
            Backward Born wavefields [in,out]
        n_pert : int32
        list_pert : str array
        pert : List_array_3d
            Perturbations array [in]
        n_grad : int32
            Number of physical parameters gradients to compute [in]
        list_grad : str array
            List of strings describing which physical parameters gradients need to be \
                computed [in]
        gradient : List_array_3d
            Gradient array to update [in,out]
        grad_type : int32
            Type of gradient computation(basic, true amplitude, adjoint, etc.) [in]
        flag_last_call : bool
            True if it is last call to gradient computation [in]
        """
        pass


    def compute_true_amp(self, mesh, field_fwd, field_bwd, gradient, sgn):
        """
        Compute true amplitude imaging condition

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh context [in]
        field_fwd : Field_tw
            Forward wavefields [in]
        field_bwd : Field_tw
            Backward wavefields [in]
        gradient : Array_3d
            Gradient array to update [in,out]
        sgn : int32
            Coeff to apply to gradient correlation terms(1 for classical true amplitude) \
                [in]
        """
        pass


    def info(self):
        """
        Print information about propagator
        """
        pass


    def formulation_order(self):
        """
        Returns formulation order of wave equation solved.

        Returns
        -------
        propagator_ac_cpml_formulation_order : int32
        """
        pass


    def source_delay(self):
        """
        Returns delay (number of time steps) between external force and resulting wavefield.

        Returns
        -------
        propagator_ac_cpml_source_delay : int32
        """
        pass


    def gradient_delay(self, grad_type):
        """
        Returns delay (number of time steps) between external force and resulting wavefield.

        Parameters
        ----------
        grad_type : int32
            Type of gradient computation (adjoint, true amplitude, etc.) [in]

        Returns
        -------
        propagator_ac_cpml_gradient_delay : int32
        """
        pass


    def reorder_forward_gradient_field(self, field, mesh, i):
        """
        Reorder (if needed) the ith wavefield used for gradient computation before returning it to the application.
        This method is called just before method `get_gradient_field` of `Field_tw`

        Parameters
        ----------
        field : Field_tw
            Field context [inout]
        mesh : Propagator_mesh
            Propagation mesh
        i : int32
            Index of wavefield to reorder [in]
        """
        pass


    def reorder_backward_gradient_field(self, field, mesh, i):
        """
        Reorder back (if needed) the ith wavefield used for gradient computation after it has been passed to the application.
        This method is called just after method `get_gradient_field` of `Field_tw`

        Parameters
        ----------
        field : Field_tw
            Field context [inout]
        mesh : Propagator_mesh
            Propagation mesh
        i : int32
            Index of wavefield to retorder [in]
        """
        pass


    def allocate_array(self, mesh, array, o_mem=None):
        """
        Allocate an array (gradient, illumination, etc.) on propagation mesh.

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagator mesh context [in]
        array : Array_3d
            Array to allocate [in,out]
        o_mem : Mem_size
            Memory counter, if provided do not allocate array, just update counters
        """
        pass


    def free(self):
        """
        Free propagator
        """
        pass


    def __del__(self):
        """
        Destructor
        """
        pass
