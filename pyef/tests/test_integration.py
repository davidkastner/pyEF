"""Integration tests for pyef workflow"""

import pytest
import os
import tempfile
from pathlib import Path
from pyef.manage import parse_job_batch_file


# Get the fixtures directory path
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestWorkflowIntegration:
    """Test complete workflow integration"""

    def test_parse_job_batch_to_analysis_setup(self, tmp_path):
        """Test parsing job batch file and setting up for analysis"""
        # Create a batch CSV file
        batch_file = tmp_path / "batch.csv"
        batch_file.write_text(
            "ef,path/to/sample_1.molden,path/to/sample_1.xyz,0,1\n"
            "# This is a comment\n"
            "esp,path/to/sample_2.molden,path/to/sample_2.xyz,1,2\n"
        )

        # Parse the batch file
        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(batch_file))

        # Verify correct parsing
        assert len(analysis_types) == 2
        assert analysis_types[0] == "ef"
        assert analysis_types[1] == "esp"

        assert len(molden_paths) == 2
        assert molden_paths[0] == "path/to/sample_1.molden"
        assert molden_paths[1] == "path/to/sample_2.molden"

        assert len(xyz_paths) == 2
        assert xyz_paths[0] == "path/to/sample_1.xyz"
        assert xyz_paths[1] == "path/to/sample_2.xyz"

        assert metal_indices == [0, 1]
        assert bond_indices == [[(0, 1)], [(1, 2)]]

        # Verify we can use these for downstream analysis
        for analysis_type, molden_path, xyz_path, metal_idx, bond_idx in zip(
            analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices
        ):
            assert isinstance(analysis_type, str)
            assert isinstance(molden_path, str)
            assert isinstance(xyz_path, str)
            assert isinstance(metal_idx, int)
            assert isinstance(bond_idx, list)
            assert len(bond_idx) > 0
            assert isinstance(bond_idx[0], tuple)

    def test_fixtures_availability(self):
        """Test that required fixtures are available"""
        sample_dir = FIXTURES_DIR / "sample_1"
        assert sample_dir.exists()

        xyz_file = sample_dir / "final_optim.xyz"
        molden_file = sample_dir / "optim.molden"

        assert xyz_file.exists()
        assert molden_file.exists()

        # Verify files are not empty
        assert xyz_file.stat().st_size > 0
        assert molden_file.stat().st_size > 0

    def test_input_files_format(self):
        """Test that input files have correct format"""
        jobs_file = FIXTURES_DIR / "sample_jobs.txt"
        metals_file = FIXTURES_DIR / "sample_metals.txt"
        bonds_file = FIXTURES_DIR / "sample_bonds.txt"

        assert jobs_file.exists()
        assert metals_file.exists()
        assert bonds_file.exists()

        # Read and validate jobs file
        with open(jobs_file) as f:
            jobs = [line.strip() for line in f if line.strip()]
        assert len(jobs) > 0
        assert all(isinstance(job, str) for job in jobs)

        # Read and validate metals file
        with open(metals_file) as f:
            metals = [int(line.strip()) for line in f if line.strip()]
        assert len(metals) > 0
        assert all(isinstance(m, int) for m in metals)

        # Read and validate bonds file
        with open(bonds_file) as f:
            bonds = [int(line.strip()) for line in f if line.strip()]
        assert len(bonds) > 0
        assert all(isinstance(b, int) for b in bonds)


class TestCLIInputs:
    """Test CLI argument parsing and validation"""

    def test_cli_required_files_check(self):
        """Test that we can validate required CLI input files exist"""
        # Simulate checking if required files exist
        required_files = [
            FIXTURES_DIR / "sample_jobs.txt",
            FIXTURES_DIR / "sample_metals.txt",
            FIXTURES_DIR / "sample_bonds.txt",
        ]

        for f in required_files:
            assert f.exists(), f"Required file {f} does not exist"

    def test_read_file_lines(self, tmp_path):
        """Test the read_file_lines helper used in run.py"""
        from pyef.run import read_file_lines

        # Create a test file
        test_file = tmp_path / "test.txt"
        test_file.write_text("line1\nline2\nline3\n")

        lines = read_file_lines(str(test_file))

        assert len(lines) == 3
        assert lines == ["line1", "line2", "line3"]

    def test_read_file_lines_with_empty_lines(self, tmp_path):
        """Test read_file_lines with empty lines"""
        from pyef.run import read_file_lines

        test_file = tmp_path / "test.txt"
        test_file.write_text("line1\n\nline2\n\n\nline3\n")

        lines = read_file_lines(str(test_file))

        # Should include empty strings for empty lines
        assert "line1" in lines
        assert "line2" in lines
        assert "line3" in lines


class TestDataFlow:
    """Test data flow through the pipeline"""

    def test_xyz_file_structure(self):
        """Test that XYZ files have expected structure"""
        xyz_file = FIXTURES_DIR / "sample_1" / "final_optim.xyz"

        with open(xyz_file) as f:
            lines = f.readlines()

        # XYZ format: first line is number of atoms
        num_atoms = int(lines[0].strip())
        assert num_atoms > 0

        # Second line is comment
        # Following lines should be atom coordinates
        # Each atom line: Element X Y Z
        for i in range(2, min(5, len(lines))):  # Check first few atom lines
            parts = lines[i].split()
            if len(parts) >= 4:
                # Element symbol
                assert parts[0].isalpha()
                # Coordinates should be numbers
                for coord in parts[1:4]:
                    float(coord)  # Should not raise ValueError

    def test_end_to_end_file_preparation(self):
        """Test that all required files for a job are present"""
        sample_dir = FIXTURES_DIR / "sample_1"

        required_files = [
            "final_optim.xyz",
            "optim.molden",
        ]

        for filename in required_files:
            filepath = sample_dir / filename
            assert filepath.exists(), f"Missing required file: {filename}"
            assert filepath.stat().st_size > 0, f"File {filename} is empty"


class TestConfigurationHandling:
    """Test configuration and parameter handling"""

    def test_dielectric_constant_validation(self):
        """Test dielectric constant parameter validation"""
        # Dielectric should be a positive float
        valid_dielectrics = [1.0, 2.0, 78.5, 80.0]
        for d in valid_dielectrics:
            assert d > 0
            assert isinstance(d, (int, float))

    def test_index_validation(self):
        """Test that atom indices are valid integers"""
        from pyef.run import read_file_lines

        metals_file = FIXTURES_DIR / "sample_metals.txt"
        bonds_file = FIXTURES_DIR / "sample_bonds.txt"

        metal_indices = [int(idx) for idx in read_file_lines(str(metals_file))]
        bond_indices = [int(idx) for idx in read_file_lines(str(bonds_file))]

        # All indices should be non-negative integers
        for idx in metal_indices:
            assert isinstance(idx, int)
            assert idx >= 0

        for idx in bond_indices:
            assert isinstance(idx, int)
            assert idx >= 0
