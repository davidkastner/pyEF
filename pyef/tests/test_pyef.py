"""
Unit and regression test for the pyef package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pyef


def test_pyef_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "pyef" in sys.modules
