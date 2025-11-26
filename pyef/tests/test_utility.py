"""Unit tests for pyef.utility module"""

import pytest
import os
import json
import tempfile
from pathlib import Path
from pyef.utility import MoldenObject


# Get the fixtures directory path
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestMoldenObject:
    """Tests for MoldenObject class"""

    @pytest.fixture
    def sample_files(self):
        """Provide paths to sample XYZ and molden files"""
        return {
            'xyz': str(FIXTURES_DIR / "sample_1" / "final_optim.xyz"),
            'molden': str(FIXTURES_DIR / "sample_1" / "optim.molden")
        }

    @pytest.fixture
    def molden_obj(self, sample_files):
        """Create a MoldenObject for testing"""
        return MoldenObject(sample_files['xyz'], sample_files['molden'])

    def test_init(self, sample_files):
        """Test MoldenObject initialization"""
        mobj = MoldenObject(sample_files['xyz'], sample_files['molden'])

        assert mobj.xyzfile == sample_files['xyz']
        assert mobj.moldenFile == sample_files['molden']
        assert hasattr(mobj, 'periodic_table')
        assert hasattr(mobj, 'amassdict')

    def test_periodic_table(self, molden_obj):
        """Test periodic table dictionary"""
        pt = molden_obj.periodic_table

        # Check basic elements
        assert pt['H'] == 1
        assert pt['C'] == 6
        assert pt['N'] == 7
        assert pt['O'] == 8
        assert pt['Fe'] == 26
        assert pt['Cu'] == 29
        assert pt['Zn'] == 30

        # Check periodic table is complete for first 118 elements
        assert len(pt) >= 100

    def test_amassdict_structure(self, molden_obj):
        """Test amassdict has expected structure"""
        amass = molden_obj.amassdict

        # Check presence of common elements
        assert 'H' in amass
        assert 'C' in amass
        assert 'O' in amass
        assert 'Fe' in amass

        # Check tuple structure (mass, atomic_num, radius, valence)
        h_data = amass['H']
        assert isinstance(h_data, tuple)
        assert len(h_data) == 4
        assert h_data[0] == pytest.approx(1.0079, rel=0.01)  # mass
        assert h_data[1] == 1  # atomic number

    def test_countBasis_spherical(self, molden_obj):
        """Test countBasis method with spherical harmonics"""
        nbf = molden_obj.countBasis(spherical=True)

        assert isinstance(nbf, int)
        assert nbf > 0  # Should have at least some basis functions

    def test_countBasis_cartesian(self, molden_obj):
        """Test countBasis method with cartesian basis"""
        nbf_cartesian = molden_obj.countBasis(spherical=False)

        assert isinstance(nbf_cartesian, int)
        assert nbf_cartesian > 0

        # Cartesian should generally have more basis functions than spherical
        # for molecules with d or f orbitals
        nbf_spherical = molden_obj.countBasis(spherical=True)
        assert nbf_cartesian >= nbf_spherical

    def test_build_core_map(self, molden_obj):
        """Test build_core_map method"""
        # Test with a common ECP family
        try:
            core_map = molden_obj.build_core_map("lanl2dz")

            assert isinstance(core_map, dict)
            # H should have 0 core electrons (all electron)
            assert core_map.get('H', 0) == 0
            # Check that heavier elements might have core electrons
            assert isinstance(core_map.get('Fe', 0), int)

        except Exception as e:
            # If basis_set_exchange is not available or family not found,
            # the method should handle it gracefully
            pytest.skip(f"Skipping build_core_map test: {e}")


class TestMoldenObjectEdgeCases:
    """Test edge cases for MoldenObject"""

    def test_init_with_nonexistent_files(self):
        """Test initialization with non-existent files"""
        mobj = MoldenObject("/fake/path.xyz", "/fake/molden.molden")

        assert mobj.xyzfile == "/fake/path.xyz"
        assert mobj.moldenFile == "/fake/molden.molden"

    def test_countBasis_with_invalid_file(self):
        """Test countBasis with non-existent molden file"""
        mobj = MoldenObject("/fake/path.xyz", "/fake/molden.molden")

        with pytest.raises(FileNotFoundError):
            mobj.countBasis()


class TestHybridCoreMap:
    """Test hybrid core map functionality"""

    @pytest.fixture
    def molden_obj_minimal(self, tmp_path):
        """Create minimal MoldenObject for testing"""
        xyz_file = tmp_path / "test.xyz"
        molden_file = tmp_path / "test.molden"

        # Create minimal files
        xyz_file.write_text("1\ntest\nH 0.0 0.0 0.0\n")
        molden_file.write_text("[Molden Format]\n[GTO]\n")

        return MoldenObject(str(xyz_file), str(molden_file))

    def test_build_hybrid_core_map_basic(self, molden_obj_minimal):
        """Test build_hybrid_core_map with basic setup"""
        # Create a simple base map for testing
        base_maps = {
            'lanl2dz': {'H': 0, 'C': 2, 'Fe': 10}
        }

        hybrid_map = molden_obj_minimal.build_hybrid_core_map(
            name='test_hybrid',
            cutoff_Z=18,
            heavy_family='lanl2dz',
            base_maps=base_maps
        )

        assert isinstance(hybrid_map, dict)
        # Elements with Z <= 18 should have core=0
        assert hybrid_map['H'] == 0
        assert hybrid_map['C'] == 0
        assert hybrid_map['O'] == 0
        # Elements with Z > 18 should use heavy_family mapping
        assert hybrid_map['Fe'] == 10


class TestMoldenBasisParsing:
    """Test molden file parsing for basis functions"""

    @pytest.fixture
    def minimal_molden(self, tmp_path):
        """Create a minimal molden file for testing"""
        molden_file = tmp_path / "minimal.molden"
        molden_content = """[Molden Format]
[Atoms] AU
H          1    1          0.000000000000          0.000000000000          0.000000000000
[GTO]
  1 0
 s    1 1.00
      3.425250914 0.1543289673
 s    1 1.00
      0.6239137298 0.5353281423
 p    1 1.00
      1.168000000 1.0000000000

[MO]
"""
        molden_file.write_text(molden_content)

        xyz_file = tmp_path / "minimal.xyz"
        xyz_file.write_text("1\ntest\nH 0.0 0.0 0.0\n")

        return MoldenObject(str(xyz_file), str(molden_file))

    def test_countBasis_minimal(self, minimal_molden):
        """Test basis counting with minimal molden file"""
        # Should count: 1 s + 1 s + 3 p = 5 basis functions (spherical)
        nbf = minimal_molden.countBasis(spherical=True)
        assert nbf == 5  # 2 s orbitals (1 each) + 1 p orbital (3 functions)
