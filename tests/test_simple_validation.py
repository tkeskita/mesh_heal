"""
Simple validation test that doesn't depend on conftest.py
"""

def test_basic_functionality():
    """Test that pytest is working at all."""
    assert 2 + 2 == 4

def test_imports():
    """Test that basic Python imports work."""
    import sys
    import os
    import tempfile
    assert sys is not None
    assert os is not None
    assert tempfile is not None

def test_pytest_markers():
    """Test that pytest markers can be used."""
    import pytest
    
    @pytest.mark.unit
    def dummy_unit_test():
        pass
        
    @pytest.mark.integration
    def dummy_integration_test():
        pass
        
    assert True