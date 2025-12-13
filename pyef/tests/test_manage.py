"""Unit tests for pyef.manage module"""

import pytest
import tempfile
import os
from pyef.manage import parse_job_batch_file


class TestParseJobBatchFile:
    """Tests for parse_job_batch_file function"""

    def test_basic_csv_parsing(self, tmp_path):
        """Test basic CSV parsing without comments"""
        csv_file = tmp_path / "jobs.csv"
        csv_file.write_text(
            "ef,path/to/file1.molden,path/to/file1.xyz,0,1\n"
            "esp,path/to/file2.molden,path/to/file2.xyz,1,2\n"
            "estab,path/to/file3.molden,path/to/file3.xyz,2,3\n"
        )

        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(csv_file))

        assert analysis_types == ["ef", "esp", "estab"]
        assert molden_paths == ["path/to/file1.molden", "path/to/file2.molden", "path/to/file3.molden"]
        assert xyz_paths == ["path/to/file1.xyz", "path/to/file2.xyz", "path/to/file3.xyz"]
        assert metal_indices == [0, 1, 2]
        assert bond_indices == [[(0, 1)], [(1, 2)], [(2, 3)]]

    def test_csv_with_comments(self, tmp_path):
        """Test CSV parsing with Python-style comments"""
        csv_file = tmp_path / "jobs.csv"
        csv_file.write_text(
            "# This is a comment\n"
            "ef,path/to/file1.molden,path/to/file1.xyz,0,1\n"
            "esp,path/to/file2.molden,path/to/file2.xyz,1,2  # inline comment\n"
            "# Another comment\n"
            "estab,path/to/file3.molden,path/to/file3.xyz,2,3\n"
        )

        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(csv_file))

        assert analysis_types == ["ef", "esp", "estab"]
        assert molden_paths == ["path/to/file1.molden", "path/to/file2.molden", "path/to/file3.molden"]
        assert xyz_paths == ["path/to/file1.xyz", "path/to/file2.xyz", "path/to/file3.xyz"]
        assert metal_indices == [0, 1, 2]
        assert bond_indices == [[(0, 1)], [(1, 2)], [(2, 3)]]

    def test_csv_with_empty_lines(self, tmp_path):
        """Test CSV parsing with empty lines"""
        csv_file = tmp_path / "jobs.csv"
        csv_file.write_text(
            "ef,path/to/file1.molden,path/to/file1.xyz,0,1\n"
            "\n"
            "esp,path/to/file2.molden,path/to/file2.xyz,1,2\n"
            "\n\n"
            "estab,path/to/file3.molden,path/to/file3.xyz,2,3\n"
        )

        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(csv_file))

        assert analysis_types == ["ef", "esp", "estab"]
        assert molden_paths == ["path/to/file1.molden", "path/to/file2.molden", "path/to/file3.molden"]
        assert xyz_paths == ["path/to/file1.xyz", "path/to/file2.xyz", "path/to/file3.xyz"]
        assert metal_indices == [0, 1, 2]

    def test_csv_with_whitespace(self, tmp_path):
        """Test CSV parsing with extra whitespace"""
        csv_file = tmp_path / "jobs.csv"
        csv_file.write_text(
            "  ef  ,  path/to/file1.molden  ,  path/to/file1.xyz  ,  0  ,  1  \n"
            "esp,path/to/file2.molden,path/to/file2.xyz,1,2\n"
            "  estab  ,  path/to/file3.molden  ,  path/to/file3.xyz  ,  2  ,  3\n"
        )

        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(csv_file))

        assert analysis_types == ["ef", "esp", "estab"]
        assert molden_paths == ["path/to/file1.molden", "path/to/file2.molden", "path/to/file3.molden"]
        assert xyz_paths == ["path/to/file1.xyz", "path/to/file2.xyz", "path/to/file3.xyz"]
        assert metal_indices == [0, 1, 2]

    def test_empty_file(self, tmp_path):
        """Test parsing an empty file"""
        csv_file = tmp_path / "empty.csv"
        csv_file.write_text("")

        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(csv_file))

        assert analysis_types == []
        assert molden_paths == []
        assert xyz_paths == []
        assert metal_indices == []
        assert bond_indices == []

    def test_only_comments(self, tmp_path):
        """Test file with only comments"""
        csv_file = tmp_path / "comments.csv"
        csv_file.write_text(
            "# Comment 1\n"
            "# Comment 2\n"
            "\n"
            "# Comment 3\n"
        )

        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(csv_file))

        assert analysis_types == []
        assert molden_paths == []
        assert xyz_paths == []
        assert metal_indices == []
        assert bond_indices == []

    def test_integer_conversion(self, tmp_path):
        """Test that metal and bond indices are properly converted to integers"""
        csv_file = tmp_path / "jobs.csv"
        csv_file.write_text("ef,path/to/file.molden,path/to/file.xyz,10,20\n")

        analysis_types, molden_paths, xyz_paths, metal_indices, bond_indices = parse_job_batch_file(str(csv_file))

        assert isinstance(metal_indices[0], int)
        assert isinstance(bond_indices[0][0][0], int)
        assert isinstance(bond_indices[0][0][1], int)
        assert metal_indices[0] == 10
        assert bond_indices[0] == [(10, 20)]

    def test_file_not_found(self):
        """Test handling of non-existent file"""
        with pytest.raises(FileNotFoundError):
            parse_job_batch_file("/nonexistent/path/to/file.csv")
