"""Module defining a Field_tw object storing wavefields of a two-way propagator"""
from pydw.m_propagator_param import Propagator_param
from pydw.m_field_tw import Field_tw
from pydw.m_propagator_mesh import Propagator_mesh
from pysabl.m_array_3d import Array_3d
from pysabl.m_mem_size import Mem_size
from pysabl.m_parallelism_context import Parallelism

class Field_template(Field_tw):

    def __init__(self):
        """
        Constructor
        """
        pass


    def get_parameters(self, param, o_grad_only=None, o_grad_type=None):
        """
        Get user parameters from command line and/or Propagator_param structure

        Parameters
        ----------
        param : Propagator_param
            Context which contains common wave equation parameters [in]
        o_grad_only : bool
            Optional flag to specify if wavefield will be used only for gradient computation [in]
        o_grad_type : int32
            Optional type of gradient [in]
        """
        pass


    def allocate(self, mesh, o_mem=None):
        """
        Allocate field arrays. If o_mem is provided, just update memory counter.

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh context [in]
        o_mem : Mem_size
            Optional Mem_size, if provided update memory estimation instead of allocating \
                arrays [in,out]
        """
        pass


    def init(self):
        """
        Initialize field arrays. Assume field has already been allocated.
        """
        pass


    def copy(self, field_in):
        """
        Copy one field to another one.

        Parameters
        ----------
        field_in : Field_template
            Field to be copied [in]
        """
        pass


    def free(self):
        """
        Free data allocated by the field
        """
        pass


    def get_n_gradient_field(self, grad_type):
        """
        Function returns the number of wavefields needed to compute model gradients

        Parameters
        ----------
        grad_type : int32
            Type of gradient computation (adjoint, true amplitude, etc.)

        Returns
        -------
        field_ac_cpml_get_n_gradient_field : int32
        """
        pass


    def get_gradient_field(self, grad_type, i, array):
        """
        Function returns the ith wavefield needed to compute model gradients.

        Parameters
        ----------
        grad_type : int32
            Type of gradient (adjoint, true amplitude, etc.)
        i : int32
            Index of field to return [in]
        array : Array_3d
            Array_3d containing pointer to ith wavefield
        """
        pass


    def get_gradient_field_size(self, mesh, grad_type, i):
        """
        Function returns the size of the ith wavefield needed to compute model \
            gradients.

        Parameters
        ----------
        mesh : Propagator_mesh
            Propagation mesh
        grad_type : int32
            Type of gradient (adjoint, true amplitude, etc.) [in]
        i : int32
            Index of field to return [in]

        Returns
        -------
        sizes : int array
            Sizes along 3 dimensions
        """
        pass


    def get_n_inner_field(self):
        """
        Function returns the number of inner wavefields needed for back-propagation

        Returns
        -------
        field_ac_cpml_get_n_inner_field : int32
        """
        pass


    def time_length(self):
        """
        Returns time length of field.
        Time length corresponds to the number of time steps between the youngest and the oldest wavefield.

        Returns
        -------
        Field_template_time_length : int32
            Time length of field
        """
        pass


    def reverse_time(self):
        """
        Reverse time order of wavefields.
        Called at the end of forward propagation when we want to start reverse time propagation.
        """
        pass


    def get_stats(self, paral):
        """
        Compute and output statistics (min/max values and presence of inf/NaN values)

        Parameters
        ----------
        paral : Parallelism
            Parallelism object containing MPI communicator

        Returns
        -------
        minfield : float32
        maxfield : float32
        flag_nan : bool
        flag_inf : bool

        Example:
        >>> field = Field_template()
        >>> minval, maxval, has_nan, has_inf = field.get_stats()
        """
        pass


    def __del__(self):
        """
        Destructor
        """
        pass
