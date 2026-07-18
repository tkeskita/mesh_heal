import sys
import os
from unittest.mock import MagicMock

# Remove the parent directory from sys.path to prevent importing the main addon
if os.path.dirname(os.path.dirname(__file__)) in sys.path:
    sys.path.remove(os.path.dirname(os.path.dirname(__file__)))

# Mock all Blender-related modules immediately to prevent import errors during test discovery  
modules_to_mock = ['bpy', 'bmesh', 'mathutils', 'numpy']

for module_name in modules_to_mock:
    if module_name not in sys.modules:
        mock_module = MagicMock()
        
        if module_name == 'bpy':
            # Comprehensive bpy mock
            mock_module.types = MagicMock()
            mock_module.props = MagicMock()
            mock_module.utils = MagicMock()
            mock_module.context = MagicMock()
            mock_module.data = MagicMock()
            mock_module.ops = MagicMock()
            
            # Mock property types
            mock_module.props.FloatProperty = MagicMock()
            mock_module.props.BoolProperty = MagicMock()
            mock_module.props.StringProperty = MagicMock()
            mock_module.props.PointerProperty = MagicMock()
            
            # Mock common types
            mock_module.types.PropertyGroup = type('PropertyGroup', (), {})
            mock_module.types.Operator = type('Operator', (), {})
            mock_module.types.Panel = type('Panel', (), {})
            mock_module.types.Scene = MagicMock()
            
        elif module_name == 'mathutils':
            # Mock mathutils
            mock_module.Vector = lambda *args: list(args) if args else [0, 0, 0]
            mock_module.Matrix = MagicMock()
            
        elif module_name == 'bmesh':
            # Mock bmesh
            mock_module.new = MagicMock()
            mock_module.from_mesh = MagicMock()
        
        elif module_name == 'numpy':
            # Use actual numpy if available, otherwise mock
            try:
                import numpy as np
                mock_module = np
            except ImportError:
                mock_module.array = MagicMock()
                mock_module.zeros = MagicMock()
        
        sys.modules[module_name] = mock_module

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock
from types import ModuleType


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    temp_path = tempfile.mkdtemp()
    yield Path(temp_path)
    shutil.rmtree(temp_path)


@pytest.fixture
def mock_bpy():
    """Mock Blender's bpy module for testing."""
    # Create a mock bpy module
    bpy_mock = MagicMock()
    
    # Mock common bpy attributes and methods
    bpy_mock.types = MagicMock()
    bpy_mock.props = MagicMock()
    bpy_mock.utils = MagicMock()
    bpy_mock.context = MagicMock()
    bpy_mock.data = MagicMock()
    bpy_mock.ops = MagicMock()
    
    # Mock property types
    bpy_mock.props.FloatProperty = MagicMock()
    bpy_mock.props.BoolProperty = MagicMock()
    bpy_mock.props.StringProperty = MagicMock()
    bpy_mock.props.PointerProperty = MagicMock()
    
    # Mock registration functions
    bpy_mock.utils.register_class = MagicMock()
    bpy_mock.utils.unregister_class = MagicMock()
    
    # Mock common types
    bpy_mock.types.PropertyGroup = type('PropertyGroup', (), {})
    bpy_mock.types.Operator = type('Operator', (), {})
    bpy_mock.types.Panel = type('Panel', (), {})
    bpy_mock.types.Scene = MagicMock()
    
    # Mock mesh operations
    bpy_mock.ops.mesh = MagicMock()
    bpy_mock.ops.object = MagicMock()
    
    return bpy_mock


@pytest.fixture
def mock_blender_environment(mock_bpy):
    """Set up a complete mock Blender environment."""
    # Add bpy to sys.modules
    sys.modules['bpy'] = mock_bpy
    
    # Mock mathutils if needed
    mathutils_mock = MagicMock()
    mathutils_mock.Vector = lambda *args: list(args) if args else [0, 0, 0]
    mathutils_mock.Matrix = MagicMock()
    sys.modules['mathutils'] = mathutils_mock
    
    # Mock bmesh if needed
    bmesh_mock = MagicMock()
    sys.modules['bmesh'] = bmesh_mock
    
    yield mock_bpy
    
    # Clean up
    modules_to_remove = ['bpy', 'mathutils', 'bmesh']
    for module in modules_to_remove:
        if module in sys.modules:
            del sys.modules[module]


@pytest.fixture
def sample_mesh_data():
    """Provide sample mesh data for testing."""
    return {
        'vertices': [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ],
        'faces': [
            [0, 1, 2, 3]
        ],
        'edges': [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0]
        ]
    }


@pytest.fixture
def mock_context():
    """Mock Blender context object."""
    context = MagicMock()
    context.active_object = MagicMock()
    context.active_object.type = 'MESH'
    context.mode = 'OBJECT'
    context.scene = MagicMock()
    context.scene.mesh_heal = MagicMock()
    return context


@pytest.fixture
def mock_mesh_object():
    """Mock Blender mesh object."""
    obj = MagicMock()
    obj.type = 'MESH'
    obj.data = MagicMock()
    obj.data.vertices = []
    obj.data.edges = []
    obj.data.polygons = []
    obj.name = 'TestMesh'
    return obj


@pytest.fixture
def mesh_heal_settings():
    """Mock mesh heal settings with default values."""
    settings = MagicMock()
    settings.vert_merge_distance = 0.001
    settings.sew_ratio_threshold = 0.3
    settings.max_abs_twist_angle = 3.0
    settings.max_abs_edge_overlap_angle = 1.0
    return settings


@pytest.fixture
def mock_operator():
    """Mock Blender operator base class."""
    operator = MagicMock()
    operator.report = MagicMock()
    operator.bl_idname = "mesh.test_operator"
    operator.bl_label = "Test Operator"
    return operator


@pytest.fixture
def cleanup_imports():
    """Clean up any module imports after test."""
    yield
    # Remove any test-specific modules that might have been imported
    modules_to_check = list(sys.modules.keys())
    for module_name in modules_to_check:
        if module_name.startswith('op_') or module_name.startswith('mesh_heal'):
            if module_name in sys.modules:
                del sys.modules[module_name]


class MockBlenderInfo:
    """Mock bl_info for testing addon metadata."""
    def __init__(self):
        self.name = "Test Addon"
        self.author = "Test Author"
        self.version = (1, 0, 0)
        self.blender = (2, 80, 0)


@pytest.fixture
def mock_bl_info():
    """Mock bl_info dictionary."""
    return MockBlenderInfo()


@pytest.fixture(autouse=True)
def setup_test_environment():
    """Automatically set up test environment for all tests."""
    # Mock bpy module before any tests run to prevent import errors
    if 'bpy' not in sys.modules:
        bpy_mock = MagicMock()
        bpy_mock.types = MagicMock()
        bpy_mock.props = MagicMock() 
        bpy_mock.utils = MagicMock()
        bpy_mock.context = MagicMock()
        bpy_mock.data = MagicMock()
        bpy_mock.ops = MagicMock()
        sys.modules['bpy'] = bpy_mock
    
    yield
    
    # Clean up any test-specific modules (but keep bpy mock for other tests)
    modules_to_remove = []
    for module_name in sys.modules:
        if module_name.startswith('op_') and 'tests' not in module_name:
            modules_to_remove.append(module_name)
    
    for module_name in modules_to_remove:
        if module_name in sys.modules:
            del sys.modules[module_name]