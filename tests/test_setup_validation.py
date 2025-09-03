"""
Validation tests to verify the testing infrastructure is working correctly.
These tests ensure that all testing components are properly configured.
"""

import pytest
import sys
from pathlib import Path


class TestInfrastructureValidation:
    """Test that the testing infrastructure is properly set up."""

    def test_pytest_runs(self):
        """Verify that pytest is working."""
        assert True

    def test_fixtures_available(self, temp_dir, mock_bpy):
        """Test that common fixtures are available and working."""
        # Test temp_dir fixture
        assert temp_dir.exists()
        assert temp_dir.is_dir()
        
        # Test mock_bpy fixture
        assert mock_bpy is not None
        assert hasattr(mock_bpy, 'types')
        assert hasattr(mock_bpy, 'props')
        assert hasattr(mock_bpy, 'utils')

    def test_mock_blender_environment(self, mock_blender_environment):
        """Test that Blender environment mocking works."""
        # Should be able to import mocked bpy
        import bpy
        assert bpy is not None
        assert hasattr(bpy, 'types')
        assert hasattr(bpy, 'props')
        assert hasattr(bpy, 'utils')

    def test_sample_fixtures(self, sample_mesh_data, mock_context, mesh_heal_settings):
        """Test that sample data fixtures work."""
        assert 'vertices' in sample_mesh_data
        assert 'faces' in sample_mesh_data
        assert len(sample_mesh_data['vertices']) == 4
        
        assert mock_context.active_object is not None
        assert mock_context.mode == 'OBJECT'
        
        assert mesh_heal_settings.vert_merge_distance == 0.001
        assert mesh_heal_settings.sew_ratio_threshold == 0.3

    @pytest.mark.unit
    def test_unit_marker(self):
        """Test that unit marker works."""
        assert True

    @pytest.mark.integration  
    def test_integration_marker(self):
        """Test that integration marker works."""
        assert True

    @pytest.mark.slow
    def test_slow_marker(self):
        """Test that slow marker works."""
        assert True

    def test_project_structure(self):
        """Test that the project structure is correct."""
        project_root = Path(__file__).parent.parent
        
        # Check that key files exist
        assert (project_root / 'pyproject.toml').exists()
        assert (project_root / 'README.md').exists()
        
        # Check addon structure
        addon_dir = project_root / 'mesh_heal_addon'
        assert addon_dir.exists()
        assert (addon_dir / '__init__.py').exists()
        
        # Check test structure
        tests_dir = project_root / 'tests'
        assert tests_dir.exists()
        assert (tests_dir / '__init__.py').exists()
        assert (tests_dir / 'conftest.py').exists()
        assert (tests_dir / 'unit').exists()
        assert (tests_dir / 'integration').exists()

    def test_coverage_works(self):
        """Test that coverage measurement can work."""
        # This is a simple test that will be measured by coverage
        def sample_function():
            return 42
        
        result = sample_function()
        assert result == 42

    def test_mocking_works(self, mock_operator):
        """Test that mocking utilities work."""
        assert mock_operator is not None
        assert hasattr(mock_operator, 'report')
        assert mock_operator.bl_idname == "mesh.test_operator"


class TestBlenderAddonStructure:
    """Tests specific to Blender addon structure."""
    
    def test_addon_files_exist(self):
        """Test that main addon files exist."""
        project_root = Path(__file__).parent.parent
        addon_dir = project_root / 'mesh_heal_addon'
        
        # Main addon files
        addon_files = [
            '__init__.py',
            'op_clean_mesh.py',
            'op_delete_overlap.py', 
            'op_fill_holes.py',
            'op_gen.py',
            'op_merge_overlapping_edges.py',
            'op_norms.py',
            'op_sew.py'
        ]
        
        for filename in addon_files:
            assert (addon_dir / filename).exists(), f"Missing addon file: {filename}"

    def test_bl_info_structure(self):
        """Test that bl_info has required structure."""
        project_root = Path(__file__).parent.parent
        addon_dir = project_root / 'mesh_heal_addon'
        init_file = addon_dir / '__init__.py'
        
        content = init_file.read_text()
        
        # Check for required bl_info keys
        required_keys = ['name', 'author', 'blender', 'location', 'description']
        for key in required_keys:
            assert f'"{key}":' in content, f"bl_info missing required key: {key}"


def test_standalone_function():
    """Test that standalone functions work."""
    assert 2 + 2 == 4