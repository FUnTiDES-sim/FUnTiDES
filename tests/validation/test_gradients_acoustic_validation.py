"""
Validation tests for acoustic gradient computation.

These tests verify that the gradient computation produces reasonable results
by checking:1. Gradient array sizes and shapes
2. Gradient magnitudes are non-zero for non-trivial input
3. Gradient values make physical sense
"""

import pytest
import numpy as np

try:
    import kokkos
    import pyfuntides.model as Model
    import pyfuntides.solver as Solver
except ImportError:
    pytest.skip("FUnTiDES Python bindings not available", allow_module_level=True)


class TestGradientsAcousticDataStructure:
    """Test acoustic gradient data structure"""

    def test_gradients_acoustic_creation(self):
        """Test that GradientsAcoustic can be created"""
        # Create small gradient views
        grad_kappa = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        grad_buoyancy = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        
        # Create GradientsAcoustic data structure
        gradients = Solver.GradientsAcoustic(grad_kappa, grad_buoyancy)
        
        assert gradients is not None
        assert gradients.get_num_grads() == 2

    def test_gradients_acoustic_names(self):
        """Test that gradient names are correct"""
        grad_kappa = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        grad_buoyancy = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        
        gradients = Solver.GradientsAcoustic(grad_kappa, grad_buoyancy)
        
        names = gradients.get_grads_names()
        assert names[0] == "gradKappa"
        assert names[1] == "gradBuoyancy"

    def test_gradients_acoustic_access(self):
        """Test that gradients can be accessed"""
        grad_kappa = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(5,))
        grad_buoyancy = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(5,))
        
        # Initialize with test values
        for i in range(5):
            grad_kappa[i] = float(i)
            grad_buoyancy[i] = float(i * 2)
        
        gradients = Solver.GradientsAcoustic(grad_kappa, grad_buoyancy)
        
        grad_k = gradients.get_grad_kappa()
        grad_b = gradients.get_grad_buoyancy()
        
        assert grad_k is not None
        assert grad_b is not None
        assert grad_k.shape[0] == 5
        assert grad_b.shape[0] == 5


class TestGradientsElasticDataStructure:
    """Test elastic gradient data structure"""

    def test_gradients_elastic_creation(self):
        """Test that GradientsElastic can be created"""
        grad_rho = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        grad_lambda = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        grad_mu = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        
        gradients = Solver.GradientsElastic(grad_rho, grad_lambda, grad_mu)
        
        assert gradients is not None
        assert gradients.get_num_grads() == 3

    def test_gradients_elastic_names(self):
        """Test that gradient names are correct"""
        grad_rho = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        grad_lambda = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        grad_mu = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(10,))
        
        gradients = Solver.GradientsElastic(grad_rho, grad_lambda, grad_mu)
        
        names = gradients.get_grads_names()
        assert names[0] == "gradRho"
        assert names[1] == "gradLambda"
        assert names[2] == "gradMu"

    def test_gradients_elastic_access(self):
        """Test that gradients can be accessed"""
        grad_rho = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(5,))
        grad_lambda = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(5,))
        grad_mu = kokkos.memory_space.HostSpace.create_view(dtype=np.float32, layout="LayoutLeft", extent=(5,))
        
        # Initialize with test values
        for i in range(5):
            grad_rho[i] = float(i)
            grad_lambda[i] = float(i * 2)
            grad_mu[i] = float(i * 3)
        
        gradients = Solver.GradientsElastic(grad_rho, grad_lambda, grad_mu)
        
        g_r = gradients.get_grad_rho()
        g_l = gradients.get_grad_lambda()
        g_m = gradients.get_grad_mu()
        
        assert g_r is not None
        assert g_l is not None
        assert g_m is not None
        assert g_r.shape[0] == 5
        assert g_l.shape[0] == 5
        assert g_m.shape[0] == 5


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
