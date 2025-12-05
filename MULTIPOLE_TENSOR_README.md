# Multipole Tensor Formalism for Electrostatic Stabilization

## Overview

This document describes the `getElectrostatic_stabilization()` function that implements the AMOEBA-style multipole tensor formalism for computing electrostatic stabilization energies.

## What's New?

### Key Addition: Direct Multipole-Multipole Interactions

The new tensor formalism computes **all** direct pairwise multipole interactions between substrate and environment atoms:

```
E_total = Σᵢ∈substrate Σⱼ∈environment M_i^T · T_ij · M_j
```

where:
- `M_i` = multipole vector for atom i: [q, μx, μy, μz, Qxx, Qxy, ...]
- `T_ij` = interaction tensor (derivatives of 1/r)
- `M_j` = multipole vector for atom j

### Included Interaction Terms

| Multipole Order | Terms Included |
|----------------|----------------|
| **Order 1** | • q×q (charge-charge) |
| **Order 2** | • q×q (charge-charge)<br>• **q×μ (charge-dipole)** ✨ NEW!<br>• **μ×μ (dipole-dipole)** ✨ NEW! |
| **Order 3** | All above PLUS:<br>• **q×Q (charge-quadrupole)** ✨ NEW!<br>• **μ×Q (dipole-quadrupole)** ✨ NEW!<br>• **Q×Q (quadrupole-quadrupole)** ✨ NEW! |

## Why Use This?

### Old Approach (Field-Based)
```python
# Environment creates fields, substrate responds
E_monopole = Σᵢ qᵢ × V(rᵢ)      # substrate charge × environment potential
E_dipole = -Σᵢ μᵢ · E(rᵢ)       # substrate dipole × environment field
```

**Limitations:**
- Only considers substrate multipoles responding to environment fields
- Ignores environment dipoles and quadrupoles
- Missing direct multipole-multipole coupling terms

### New Approach (Tensor-Based)
```python
# Direct pairwise multipole-multipole interactions
E_total = Σᵢ∈substrate Σⱼ∈environment M_i^T · T_ij · M_j
```

**Advantages:**
- ✅ Complete: includes ALL multipole-multipole interactions
- ✅ Symmetric: properly treats both substrate and environment multipoles
- ✅ Accurate: captures charge-dipole and dipole-dipole coupling
- ✅ Extensible: easy to add higher-order terms

## Usage

### Basic Example

```python
from pyef.analysis import Electrostatics
import numpy as np

# Setup (new simplified syntax)
estat = Electrostatics(
    ['AXEDEN/scr'],  # Complete folder path
    molden_filename='input.molden',
    xyzfilename='xyz.xyz',
    hasECP=True
)

# Define substrate atoms
substrate_idxs = np.arange(0, 11)  # Atoms 0-10

# Run calculation with dipole-dipole interactions using tensor formalism
df = estat.getElectrostatic_stabilization(
    multiwfn_path='/path/to/Multiwfn',
    multiwfn_module='module load multiwfn',
    atmrad_path='/path/to/atmrad',
    substrate_idxs=substrate_idxs,
    charge_type='Hirshfeld_I',
    multipole_order=2,  # Include dipole-dipole terms
    decompose_atomwise=True  # Get per-atom breakdown
)

# Results
print(f"Total energy: {df[0]['Total_Energy_kcal_mol'].iloc[0]:.4f} kcal/mol")
```

### Choosing Multipole Order

**Order 1** (monopole only):
- Use when: You want to reproduce the classic q×q interaction
- Fastest, but least accurate for polar systems

**Order 2** (+ dipole) **← RECOMMENDED**:
- Use when: You have polar molecules or want accurate intermolecular interactions
- Includes critical dipole-dipole interactions
- Good balance of accuracy and speed

**Order 3** (+ quadrupole):
- Use when: Highest accuracy is needed
- Important for systems with strong quadrupole moments (e.g., CO, benzene)
- Slower due to more tensor elements

## Function Signature

```python
def getElectrostatic_stabilization(
    self,
    multiwfn_path: str,
    multiwfn_module: str,
    atmrad_path: str,
    substrate_idxs: list,
    charge_type: str = 'Hirshfeld_I',
    name_dataStorage: str = 'estaticFile_tensor',
    env_idxs: list | None = None,
    decompose_atomwise: bool = False,
    visualize: bool | None = None,
    multipole_order: int = 2
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `multiwfn_path` | str | Path to Multiwfn executable |
| `multiwfn_module` | str | Module load command for Multiwfn |
| `atmrad_path` | str | Path to atmrad directory |
| `substrate_idxs` | list | List of substrate atom indices (per file) |
| `charge_type` | str | Charge partitioning scheme (default: 'Hirshfeld_I') |
| `name_dataStorage` | str | Output filename prefix |
| `env_idxs` | list\|None | Environment atom indices (None = all non-substrate) |
| `decompose_atomwise` | bool | Return per-atom breakdown (default: False) |
| `visualize` | bool\|None | Create PDB visualization (default: None) |
| `multipole_order` | int | 1=monopole, 2=+dipole, 3=+quadrupole (default: 2) |

### Returns

- If `decompose_atomwise=False`: DataFrame with total energies
- If `decompose_atomwise=True`: Tuple of (total_df, atomwise_df)

## Test Scripts

Two test scripts are provided:

### 1. `test_multipole_tensor.py`
Comprehensive comparison of all methods:
- Field-based monopole (original)
- Field-based dipole correction
- Tensor monopole
- Tensor order 2 (recommended)
- Tensor order 3 (quadrupoles)

Run with:
```bash
cd /path/to/your/calculations
python /home/gridsan/mmanetsch/pyEF/test_multipole_tensor.py
```

### 2. `example_tensor_usage.py`
Simple usage examples showing:
- Basic usage with multipole_order=2
- Comparing different orders
- Processing results

Run with:
```bash
cd /path/to/your/calculations
python /home/gridsan/mmanetsch/pyEF/example_tensor_usage.py
```

## Output Files

The function creates:

1. **`{name_dataStorage}.csv`**: Total energies
   - Columns: `Total_Energy_kcal_mol`, `Multipole_Order`, `Num_Substrate_Atoms`, etc.

2. **`{name_dataStorage}_atomwise.csv`** (if `decompose_atomwise=True`): Per-atom contributions
   - Shows which substrate atoms contribute most to stabilization
   - Includes multipole moments for each atom

3. **PDB files** (if `visualize=True`):
   - Contributions stored in B-factor column for visualization

## When to Use Each Method

| Scenario | Recommended Method | Order |
|----------|-------------------|-------|
| Quick approximation | `getElectrostatic_stabilization()` | 1 |
| Non-polar systems | `getElectrostatic_stabilization()` | 1 |
| **Polar molecules** | `getElectrostatic_stabilization()` | **2** |
| **Intermolecular interactions** | `getElectrostatic_stabilization()` | **2** |
| Highest accuracy | `getElectrostatic_stabilization()` | 3 |
| Systems with strong quadrupoles | `getElectrostatic_stabilization()` | 3 |

## Example Results

Typical output showing importance of higher-order terms:

```
RESULTS SUMMARY
================================================================================

Monopole-only energy:                -45.2341 kcal/mol
Full order-2 energy (tensor):        -52.8765 kcal/mol
Higher-order contribution (q×μ+μ×μ):  -7.6424 kcal/mol
Relative importance:                  16.89%

⚠️  WARNING: Higher-order terms are significant!
    The monopole-only approximation may not be accurate.
    Consider using the tensor formalism with multipole_order=2 or 3.
```

## Technical Details

### Interaction Tensor Elements

The interaction tensor `T_ij` contains spatial derivatives of the Coulomb potential:

**Order 1** (1×1):
```
T[0,0] = k/r  (monopole-monopole)
```

**Order 2** (4×4):
```
T[0,0]   = k/r                      (monopole-monopole)
T[0,1:4] = -k·r̂/r²                  (monopole-dipole)
T[1:4,0] = k·r̂/r²                   (dipole-monopole)
T[1:4,1:4] = k·(3r̂⊗r̂ - I)/r³        (dipole-dipole)
```

**Order 3** (10×10): Adds quadrupole terms with third and fourth derivatives

### Unit Conversions

- **Charges**: Dimensionless (in units of elementary charge)
- **Dipoles**: Bohr → meters (×5.29e-11)
- **Quadrupoles**: Bohr² → m² (×2.80e-21)
- **Energy**: Joules → kcal/mol (×6.02e23/4184)

## References

1. Ren, P. & Ponder, J. W. "Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation" *J. Phys. Chem. B* **2003**, 107, 5933-5947.
   - Original AMOEBA force field paper describing multipole tensor formalism

2. Stone, A. J. "The Theory of Intermolecular Forces" 2nd ed., Oxford University Press, 2013.
   - Theoretical background on multipole expansions

3. Multiwfn documentation: http://sobereva.com/multiwfn/
   - Source of atomic multipole moments

## Questions?

For issues or questions about this implementation, refer to:
- Original implementation: `/home/gridsan/mmanetsch/pyEF/pyef/analysis.py` (lines 3396-3829)
- Test scripts: `test_multipole_tensor.py`, `example_tensor_usage.py`
- AMOEBA paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC2918242/
