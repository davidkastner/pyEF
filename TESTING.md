# Testing Guide for pyEF

This guide explains how to set up your development environment and run tests for pyEF.

## Quick Start with Conda (Recommended)

The easiest way to install all dependencies including openbabel is using conda:

```bash
# Create conda environment from environment.yml
conda env create -f environment.yml

# Activate the environment
conda activate pyef

# Install pyef in development mode
pip install -e .

# Run tests
pytest pyef/tests/ -v
```

## Manual Installation

If you prefer to install dependencies manually:

### Using Conda

```bash
# Create a new conda environment
conda create -n pyef python=3.10

# Activate the environment
conda activate pyef

# Install openbabel and other dependencies
conda install -c conda-forge openbabel numpy pandas scipy pytest pytest-cov basis_set_exchange

# Install biopandas with pip
pip install biopandas

# Install pyef in development mode
pip install -e .
```

### Using Pip (Advanced)

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y libopenbabel-dev openbabel swig

# Install Python dependencies
pip install openbabel-wheel numpy pandas scipy biopandas basis_set_exchange pytest pytest-cov

# Install pyef in development mode
pip install -e .
```

## Running Tests

### Run All Tests

```bash
pytest pyef/tests/ -v
```

### Run Specific Test Files

```bash
# Run only manage tests (no dependencies required)
pytest pyef/tests/test_manage.py -v

# Run integration tests
pytest pyef/tests/test_integration.py -v

# Run geometry tests (requires openbabel)
pytest pyef/tests/test_geometry.py -v

# Run utility tests (requires openbabel)
pytest pyef/tests/test_utility.py -v
```

### Run with Coverage

```bash
# Note: pytest-cov must be installed
pytest pyef/tests/ -v --cov=pyef --cov-report=term-missing --cov-report=html
```

## Test Structure

```
pyef/tests/
├── __init__.py
├── test_manage.py         # CSV parsing and file management (8 tests)
├── test_geometry.py       # Geometry class and XYZ handling
├── test_utility.py        # MoldenObject and ECP handling
├── test_integration.py    # End-to-end workflow tests (10 tests)
└── fixtures/              # Test data
    ├── sample_1/
    │   ├── final_optim.xyz
    │   ├── optim.molden
    │   └── ...
    ├── sample_jobs.txt
    ├── sample_metals.txt
    └── sample_bonds.txt
```

## Test Results

Current test status:
- ✅ `test_manage.py`: 8/8 passing
- ✅ `test_integration.py`: 7/10 passing (3 require openbabel)
- ⚠️  `test_geometry.py`: Requires openbabel
- ⚠️  `test_utility.py`: Requires openbabel

## Continuous Integration

Tests run automatically on GitHub Actions for:
- Python 3.8, 3.9, 3.10, 3.11
- Every push and pull request to main/develop branches

See `.github/workflows/tests.yml` for CI configuration.

## Troubleshooting

### ImportError: No module named 'openbabel'

This means openbabel is not installed. Use conda to install it:

```bash
conda install -c conda-forge openbabel
```

### Tests fail with "FileNotFoundError"

Make sure you're running tests from the project root directory:

```bash
cd /path/to/pyEF
pytest pyef/tests/ -v
```

### pytest not found

Install pytest:

```bash
conda install pytest pytest-cov
# or
pip install pytest pytest-cov
```

## Adding New Tests

When adding new tests:

1. Create test files in `pyef/tests/` following the naming convention `test_*.py`
2. Add test fixtures to `pyef/tests/fixtures/` if needed
3. Run tests locally before committing
4. Ensure tests pass in CI before merging

## Test Coverage Goals

We aim for:
- 80%+ code coverage for core modules (manage, geometry, utility)
- 100% coverage for critical bug fixes
- All new features must include tests
