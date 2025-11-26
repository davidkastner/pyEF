"""Unit tests for pyef.geometry module"""

import pytest
import os
import pandas as pd
from pathlib import Path
from pyef.geometry import Geometry


# Get the fixtures directory path
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestGeometry:
    """Tests for Geometry class"""

    @pytest.fixture
    def sample_xyz_file(self):
        """Provide path to sample XYZ file"""
        return str(FIXTURES_DIR / "sample_1" / "final_optim.xyz")

    @pytest.fixture
    def geometry_obj(self, sample_xyz_file):
        """Create a Geometry object for testing"""
        return Geometry(sample_xyz_file)

    def test_init(self, sample_xyz_file):
        """Test Geometry initialization"""
        geom = Geometry(sample_xyz_file)
        assert geom.xyzfile == sample_xyz_file
        assert hasattr(geom, 'periodic_table')
        assert hasattr(geom, 'amassdict')
        assert isinstance(geom.periodic_table, dict)

    def test_getGeomInfo(self, geometry_obj):
        """Test getGeomInfo method returns a DataFrame"""
        df = geometry_obj.getGeomInfo()

        assert isinstance(df, pd.DataFrame)
        assert 'Atom' in df.columns
        assert 'X' in df.columns
        assert 'Y' in df.columns
        assert 'Z' in df.columns
        assert len(df) > 0

    def test_getGeomInfo_structure(self, geometry_obj):
        """Test that getGeomInfo returns expected structure"""
        df = geometry_obj.getGeomInfo()

        # Check that we have the expected atoms (O, C, H from the sample)
        atoms = set(df['Atom'].unique())
        assert 'O' in atoms or 'C' in atoms or 'H' in atoms

        # Check that coordinates are numeric
        assert df['X'].dtype in [float, 'float64']
        assert df['Y'].dtype in [float, 'float64']
        assert df['Z'].dtype in [float, 'float64']

    def test_getBondedAtoms(self, geometry_obj):
        """Test getBondedAtoms method"""
        # Get bonded atoms for first atom (index 0)
        bonded = geometry_obj.getBondedAtoms(0)

        assert isinstance(bonded, list)
        # Check that bonded atoms are integers
        if len(bonded) > 0:
            assert all(isinstance(idx, (int, np.integer)) for idx in bonded)

    def test_periodic_table_completeness(self, geometry_obj):
        """Test that periodic table contains common elements"""
        pt = geometry_obj.periodic_table

        # Check for common elements
        assert pt['H'] == 1
        assert pt['C'] == 6
        assert pt['N'] == 7
        assert pt['O'] == 8
        assert pt['Fe'] == 26
        assert pt['Cu'] == 29

    def test_amassdict_structure(self, geometry_obj):
        """Test that amassdict has expected structure"""
        amass = geometry_obj.amassdict

        # Check common elements
        assert 'H' in amass
        assert 'C' in amass
        assert 'O' in amass

        # Check that values are tuples with expected structure
        h_data = amass['H']
        assert isinstance(h_data, tuple)
        assert len(h_data) == 4


import numpy as np


class TestGeometryEdgeCases:
    """Test edge cases and error handling"""

    def test_invalid_file_path(self):
        """Test initialization with invalid file path"""
        # This should not raise an error at init, but will fail when methods are called
        geom = Geometry("/nonexistent/file.xyz")
        assert geom.xyzfile == "/nonexistent/file.xyz"

    def test_getGeomInfo_with_invalid_file(self):
        """Test getGeomInfo with non-existent file"""
        geom = Geometry("/nonexistent/file.xyz")
        with pytest.raises(FileNotFoundError):
            geom.getGeomInfo()


class TestGeometryCalculations:
    """Test geometric calculations"""

    @pytest.fixture
    def simple_xyz(self, tmp_path):
        """Create a simple XYZ file for testing"""
        xyz_file = tmp_path / "simple.xyz"
        xyz_content = """2
Simple molecule
H 0.0 0.0 0.0
H 1.0 0.0 0.0
"""
        xyz_file.write_text(xyz_content)
        return str(xyz_file)

    def test_simple_molecule(self, simple_xyz):
        """Test with a simple two-atom molecule"""
        geom = Geometry(simple_xyz)
        df = geom.getGeomInfo()

        assert len(df) == 2
        assert all(df['Atom'] == 'H')
        assert df.iloc[0]['X'] == 0.0
        assert df.iloc[1]['X'] == 1.0
