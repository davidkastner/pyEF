# pyEF Function Migration Guide

This guide shows how to migrate from deprecated functions to the new consolidated API.

## Summary of Changes

The following functions have been **removed** and replaced with unified methods:

### ESP Functions (3 removed → 1 unified)
- ❌ `getESPData()`
- ❌ `getESPMultipole()`
- ❌ `getESPDecay()`
- ✅ `getESP()` (unified replacement)

### E-field Functions (3 removed → 1 unified)
- ❌ `getEfield_acrossBond()`
- ❌ `getEfield_decomposable()`
- ❌ `getEFieldMultipole()`
- ✅ `getEfield()` (unified replacement)

### Coordination Shell Functions (2 removed → 1 unified)
- ❌ `esp_first_coord()`
- ❌ `esp_second_coord()`
- ✅ `esp_coord_shell()` (unified replacement)

---

## ESP Function Migrations

### 1. Replace `getESPData()` with `getESP()`

**Old Code:**
```python
df = electrostatics.getESPData(
    charge_types=['Hirshfeld', 'Becke', 'CHELPG'],
    ESPdata_filename='esp_output',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    dielectric=1
)
```

**New Code:**
```python
df = electrostatics.getESP(
    charge_types=['Hirshfeld', 'Becke', 'CHELPG'],
    ESPdata_filename='esp_output',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    use_multipole=False,              # NEW: Use monopole charges
    include_decay=False,               # NEW: No decay analysis
    include_coord_shells=False,        # NEW: No coordination shell analysis
    dielectric=1
)
```

**Key Changes:**
- Added `use_multipole=False` to specify monopole-based ESP calculation
- Added `include_decay=False` to disable distance decay analysis
- Added `include_coord_shells=False` to disable coordination shell analysis

---

### 2. Replace `getESPMultipole()` with `getESP()`

**Old Code:**
```python
df = electrostatics.getESPMultipole(
    charge_type='Hirshfeld_I',         # Note: singular, not plural
    ESPdata_filename='esp_multipole',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    dielectric=1
)
```

**New Code:**
```python
df = electrostatics.getESP(
    charge_types='Hirshfeld_I',        # Can be string or list
    ESPdata_filename='esp_multipole',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    use_multipole=True,                # NEW: Use multipole expansion
    include_decay=False,                # NEW: No decay analysis
    include_coord_shells=False,         # NEW: No coordination shell analysis
    dielectric=1
)
```

**Key Changes:**
- Changed parameter name from `charge_type` (singular) to `charge_types` (plural)
- Added `use_multipole=True` for multipole-based ESP calculation
- Added `include_decay=False` and `include_coord_shells=False`

---

### 3. Replace `getESPDecay()` with `getESP()`

**Old Code:**
```python
df = electrostatics.getESPDecay(
    charge_types=['Hirshfeld', 'Becke'],
    ESPdata_filename='esp_decay',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    dielectric=1
)
```

**New Code:**
```python
df = electrostatics.getESP(
    charge_types=['Hirshfeld', 'Becke'],
    ESPdata_filename='esp_decay',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    use_multipole=False,                # Monopole-based
    include_decay=True,                 # NEW: Include distance decay analysis
    include_coord_shells=True,          # NEW: Include coordination shell ESP
    dielectric=1
)
```

**Key Changes:**
- Added `use_multipole=False` for monopole-based calculation
- Added `include_decay=True` to enable distance-dependent ESP analysis
- Added `include_coord_shells=True` to calculate ESP at coordination shells

---

## E-field Function Migrations

### 4. Replace `getEfield_acrossBond()` with `getEfield()`

**Old Code:**
```python
df = electrostatics.getEfield_acrossBond(
    charge_type='Hirshfeld',
    Efielddata_filename='efield_bond',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    multipole_bool=False,
    input_bond_indices=[[0, 1], [0, 2]],
    dielectric=1
)
```

**New Code:**
```python
df = electrostatics.getEfield(
    charge_types='Hirshfeld',          # Parameter name changed to plural
    Efielddata_filename='efield_bond',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    multipole_bool=False,
    input_bond_indices=[[0, 1], [0, 2]],
    auto_find_bonds=False,              # NEW: Use manual bond indices
    decompose_atomwise=True,            # NEW: Include atom-wise decomposition
    visualize=None,                     # NEW: Optional visualization control
    dielectric=1
)
```

**Key Changes:**
- Changed `charge_type` → `charge_types`
- Added `auto_find_bonds=False` to use manual bond specification
- Added `decompose_atomwise=True` to enable per-atom E-field contributions
- Added `visualize=None` (uses config settings for visualization)

---

### 5. Replace `getEfield_decomposable()` with `getEfield()`

**Old Code:**
```python
df = electrostatics.getEfield_decomposable(
    charge_type='Hirshfeld_I',
    Efielddata_filename='efield_decomp',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    multipole_bool=True,
    input_bond_indices=[[0, 1]],
    dielectric=1
)
```

**New Code:**
```python
df = electrostatics.getEfield(
    charge_types='Hirshfeld_I',        # Parameter name changed
    Efielddata_filename='efield_decomp',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    multipole_bool=True,
    input_bond_indices=[[0, 1]],
    auto_find_bonds=False,              # NEW: Manual bond specification
    decompose_atomwise=True,            # NEW: Enable decomposition
    visualize=None,                     # NEW: Use config visualization settings
    dielectric=1
)
```

**Note:** This is nearly identical to `getEfield_acrossBond()` migration. Both had the same functionality.

---

### 6. Replace `getEFieldMultipole()` with `getEfield()`

**Old Code:**
```python
df = electrostatics.getEFieldMultipole(
    Efield_data_filename='efield_multipole',
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    input_bond_indices=[[0, 1], [2, 3]],
    excludeAtoms=[[5, 6], [7, 8]],
    polarization_scheme='Hirshfeld_I'
)
```

**New Code:**
```python
# If you need to exclude atoms, set them BEFORE calling getEfield:
electrostatics.config['excludeAtomfromEcalc'] = [[5, 6], [7, 8]]

df = electrostatics.getEfield(
    charge_types='Hirshfeld_I',        # Was 'polarization_scheme'
    Efielddata_filename='efield_multipole',  # Parameter name changed
    multiwfn_module='module load multiwfn',
    multiwfn_path='/path/to/Multiwfn',
    atmrad_path='/path/to/atmrad',
    multipole_bool=True,                # NEW: Enable multipole calculation
    input_bond_indices=[[0, 1], [2, 3]],
    auto_find_bonds=False,              # NEW: Manual bonds
    decompose_atomwise=False,           # NEW: No atom-wise decomposition (matching old behavior)
    visualize=None,                     # NEW: Use config settings
    dielectric=1
)
```

**Key Changes:**
- Renamed `polarization_scheme` → `charge_types`
- Renamed `Efield_data_filename` → `Efielddata_filename`
- `excludeAtoms` must now be set via `config['excludeAtomfromEcalc']` before calling
- Added `multipole_bool=True` explicitly
- Added `auto_find_bonds=False` and `decompose_atomwise=False`

---

## Coordination Shell Function Migrations

### 7. Replace `esp_first_coord()` with `esp_coord_shell()`

**Old Code:**
```python
esp_value = electrostatics.esp_first_coord(
    metal_idx=42,
    charge_file='charges.txt',
    path_to_xyz='structure.xyz'
)
```

**New Code:**
```python
esp_value = electrostatics.esp_coord_shell(
    metal_idx=42,
    charge_file='charges.txt',
    path_to_xyz='structure.xyz',
    n_shells=1                          # NEW: Specify 1st coordination shell
)
```

**Key Changes:**
- Added `n_shells=1` parameter to specify first coordination shell

---

### 8. Replace `esp_second_coord()` with `esp_coord_shell()`

**Old Code:**
```python
esp_value = electrostatics.esp_second_coord(
    metal_idx=42,
    charge_file='charges.txt',
    path_to_xyz='structure.xyz'
)
```

**New Code:**
```python
esp_value = electrostatics.esp_coord_shell(
    metal_idx=42,
    charge_file='charges.txt',
    path_to_xyz='structure.xyz',
    n_shells=2                          # NEW: Specify 1st + 2nd coordination shells
)
```

**Key Changes:**
- Added `n_shells=2` parameter to include first and second coordination shells

**Bonus - Third Shell (New Feature):**
```python
# You can now calculate third coordination shell (not possible before!)
esp_value = electrostatics.esp_coord_shell(
    metal_idx=42,
    charge_file='charges.txt',
    path_to_xyz='structure.xyz',
    n_shells=3                          # Would include 1st, 2nd, and 3rd shells
)
# Note: n_shells > 2 will raise NotImplementedError in current version
```

---

## Quick Reference Table

| Old Function | New Function | Key New Parameters |
|-------------|--------------|-------------------|
| `getESPData()` | `getESP()` | `use_multipole=False, include_decay=False, include_coord_shells=False` |
| `getESPMultipole()` | `getESP()` | `use_multipole=True, include_decay=False, include_coord_shells=False` |
| `getESPDecay()` | `getESP()` | `use_multipole=False, include_decay=True, include_coord_shells=True` |
| `getEfield_acrossBond()` | `getEfield()` | `auto_find_bonds=False, decompose_atomwise=True` |
| `getEfield_decomposable()` | `getEfield()` | `auto_find_bonds=False, decompose_atomwise=True` |
| `getEFieldMultipole()` | `getEfield()` | `multipole_bool=True, auto_find_bonds=False, decompose_atomwise=False` |
| `esp_first_coord()` | `esp_coord_shell()` | `n_shells=1` |
| `esp_second_coord()` | `esp_coord_shell()` | `n_shells=2` |

---

## Benefits of New Unified Functions

### `getESP()` Benefits:
- **Single interface** for monopole and multipole ESP calculations
- **Flexible analysis** - enable/disable decay and coordination shell analysis as needed
- **Consistent API** - same parameters regardless of calculation type
- **Better documentation** - comprehensive NumPy-style docstrings

### `getEfield()` Benefits:
- **Unified workflow** for all E-field calculation types
- **Flexible bond specification** - manual or automatic bond finding
- **Optional decomposition** - control atom-wise contribution analysis
- **Visualization control** - programmatic control over PDB generation
- **Cleaner code** - reduced duplication, easier maintenance

### `esp_coord_shell()` Benefits:
- **Extensible design** - easily add 3rd, 4th coordination shells
- **Single function** - no need to remember separate function names
- **Clear intent** - `n_shells` parameter makes calculation type explicit
- **Reduced code** - eliminates duplicate logic

---

## Common Migration Patterns

### Pattern 1: Simple Replacement
Most migrations just add new parameters with sensible defaults:
```python
# Old
getESPData(args...)

# New
getESP(args..., use_multipole=False, include_decay=False, include_coord_shells=False)
```

### Pattern 2: Parameter Renaming
Watch for singular → plural parameter name changes:
```python
# Old
charge_type='Hirshfeld'      # Singular
polarization_scheme='...'     # Different name

# New
charge_types='Hirshfeld'     # Plural (can also be list)
charge_types='...'            # Consistent naming
```

### Pattern 3: Config-based Settings
Some parameters moved to config object:
```python
# Old
getEFieldMultipole(..., excludeAtoms=[[1,2]])

# New
electrostatics.config['excludeAtomfromEcalc'] = [[1,2]]
getEfield(...)
```

---

## Need Help?

If you encounter issues during migration:

1. Check the function signature in `analysis.py` for the exact parameters
2. Look at the NumPy-style docstrings for parameter descriptions
3. The new unified functions accept **more** parameters than the old ones - extra flexibility!
4. All old parameter values should work - you just need to add new ones

## Location in Code

All consolidated functions are in:
- **File:** `/home/gridsan/mmanetsch/pyEF/pyef/analysis.py`
- **Class:** `Electrostatics`
- **Functions:**
  - `getESP()` - Line 1594
  - `getEfield()` - Line 1747
  - `esp_coord_shell()` - Line 1391
