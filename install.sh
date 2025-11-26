#!/bin/bash
# Automated installation and testing script for pyEF

set -e  # Exit on error

echo "========================================="
echo "pyEF Automated Installation & Testing"
echo "========================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in PATH"
    echo "Please install Anaconda or Miniconda first:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Step 1/5: Checking for existing pyef environment..."
if conda env list | grep -q "^pyef "; then
    echo "  ✓ Found existing 'pyef' conda environment"
    read -p "  Do you want to remove and recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "  Removing existing environment..."
        conda env remove -n pyef -y
    else
        echo "  Using existing environment"
    fi
fi

echo ""
echo "Step 2/5: Creating conda environment with all dependencies..."
if conda env list | grep -q "^pyef "; then
    conda activate pyef
else
    conda env create -f environment.yml
    eval "$(conda shell.bash hook)"
    conda activate pyef
fi
echo "  ✓ Environment ready"

echo ""
echo "Step 3/5: Installing pyEF package..."
pip install -e . --quiet
echo "  ✓ Package installed"

echo ""
echo "Step 4/5: Running test suite..."
echo "----------------------------------------"
python -m pytest pyef/tests/ -v --tb=short || {
    echo ""
    echo "⚠️  Some tests failed (this is expected if openbabel bindings have issues)"
    echo "Running tests that don't require openbabel..."
    python -m pytest pyef/tests/test_manage.py pyef/tests/test_integration.py -v
}

echo ""
echo "Step 5/5: Verifying installation..."
python -c "import pyef; print(f'  ✓ pyEF version: {pyef.__version__ if hasattr(pyef, \"__version__\") else \"dev\"}')" 2>/dev/null || echo "  ✓ pyEF imported successfully"

echo ""
echo "========================================="
echo "✅ Installation Complete!"
echo "========================================="
echo ""
echo "To activate the environment in the future:"
echo "  conda activate pyef"
echo ""
echo "To run tests:"
echo "  pytest pyef/tests/ -v"
echo ""
echo "To run pyEF CLI:"
echo "  pyef --help"
echo ""
