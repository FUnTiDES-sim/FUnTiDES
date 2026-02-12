# test_face_connectivity.py
import pyfuntides.model as Model
import pytest


# =============================================================================
# FaceConnectivity Bindings Tests
# =============================================================================

def test_face_connectivity_data_classes_exist():
    """Test that all FaceConnectivityUnstructData classes are bound"""
    assert hasattr(Model, 'FaceConnectivityUnstructData_f32_i32')
    assert hasattr(Model, 'FaceConnectivityUnstructData_f64_i32')
    assert hasattr(Model, 'FaceConnectivityUnstructData_f32_i64')
    assert hasattr(Model, 'FaceConnectivityUnstructData_f64_i64')


def test_face_connectivity_class_classes_exist():
    """Test that all FaceConnectivityUnstruct classes are bound"""
    assert hasattr(Model, 'FaceConnectivityUnstruct_f32_i32')
    assert hasattr(Model, 'FaceConnectivityUnstruct_f64_i32')
    assert hasattr(Model, 'FaceConnectivityUnstruct_f32_i64')
    assert hasattr(Model, 'FaceConnectivityUnstruct_f64_i64')


def test_face_connectivity_data_instantiation():
    """Test that FaceConnectivityUnstructData can be instantiated"""
    fc_data_f32_i32 = Model.FaceConnectivityUnstructData_f32_i32()
    fc_data_f64_i32 = Model.FaceConnectivityUnstructData_f64_i32()
    fc_data_f32_i64 = Model.FaceConnectivityUnstructData_f32_i64()
    fc_data_f64_i64 = Model.FaceConnectivityUnstructData_f64_i64()
    
    assert fc_data_f32_i32 is not None
    assert fc_data_f64_i32 is not None
    assert fc_data_f32_i64 is not None
    assert fc_data_f64_i64 is not None


def test_face_connectivity_data_members_accessible():
    """Test that FaceConnectivityUnstructData members are accessible"""
    fc_data = Model.FaceConnectivityUnstructData_f32_i32()
    
    # Should have all public members from the struct
    assert hasattr(fc_data, 'n_faces')
    assert hasattr(fc_data, 'ndofs_per_face')
    assert hasattr(fc_data, 'elem_to_faces')
    assert hasattr(fc_data, 'face_dofs')
    assert hasattr(fc_data, 'face_elem_owner')
    assert hasattr(fc_data, 'face_elem_neighbor')
    assert hasattr(fc_data, 'face_local_owner')
    assert hasattr(fc_data, 'face_local_neighbor')
    
    # Default values should be 0 or empty
    assert fc_data.n_faces == 0
    assert fc_data.ndofs_per_face == 0


def test_face_connectivity_data_members_writable():
    """Test that FaceConnectivityUnstructData members can be written"""
    fc_data = Model.FaceConnectivityUnstructData_f32_i32()
    
    # Should be able to set scalar members
    fc_data.n_faces = 100
    fc_data.ndofs_per_face = 9
    
    assert fc_data.n_faces == 100
    assert fc_data.ndofs_per_face == 9


def test_face_connectivity_class_instantiation_default():
    """Test that FaceConnectivityUnstruct can be instantiated with default constructor"""
    fc_f32_i32 = Model.FaceConnectivityUnstruct_f32_i32()
    fc_f64_i32 = Model.FaceConnectivityUnstruct_f64_i32()
    fc_f32_i64 = Model.FaceConnectivityUnstruct_f32_i64()
    fc_f64_i64 = Model.FaceConnectivityUnstruct_f64_i64()
    
    assert fc_f32_i32 is not None
    assert fc_f64_i32 is not None
    assert fc_f32_i64 is not None
    assert fc_f64_i64 is not None


def test_face_connectivity_class_instantiation_from_data():
    """Test that FaceConnectivityUnstruct can be instantiated from Data"""
    fc_data = Model.FaceConnectivityUnstructData_f32_i32()
    fc_data.n_faces = 50
    fc_data.ndofs_per_face = 9
    
    fc = Model.FaceConnectivityUnstruct_f32_i32(fc_data)
    assert fc is not None
    
    # If methods are bound, test them
    if hasattr(fc, 'get_number_of_faces'):
        assert fc.get_number_of_faces() == 50


def test_face_connectivity_class_methods_exist():
    """Test that FaceConnectivityUnstruct methods are bound"""
    fc = Model.FaceConnectivityUnstruct_f32_i32()
    
    # These methods should be bound (based on bindings_face_connectivity.h)
    expected_methods = [
        'get_number_of_faces',
        'get_global_face',
        'get_global_node_from_face',
        'is_boundary_face',
        'build'
    ]
    
    for method in expected_methods:
        assert hasattr(fc, method), f"Method {method} not found on FaceConnectivityUnstruct"


@pytest.mark.parametrize("suffix", ["f32_i32", "f64_i32", "f32_i64", "f64_i64"])
def test_face_connectivity_integration_with_model_unstruct_data(suffix):
    """Test that FaceConnectivityUnstructData can be assigned to ModelUnstructData"""
    # Get the appropriate classes
    fc_data_cls = getattr(Model, f'FaceConnectivityUnstructData_{suffix}')
    model_data_cls = getattr(Model, f'ModelUnstructData_{suffix}')
    
    # Create instances
    fc_data = fc_data_cls()
    fc_data.n_faces = 42
    fc_data.ndofs_per_face = 9
    
    model_data = model_data_cls()
    
    # Should be able to assign face_connectivity
    model_data.face_connectivity = fc_data
    
    # Should be able to read it back
    assert model_data.face_connectivity.n_faces == 42
    assert model_data.face_connectivity.ndofs_per_face == 9