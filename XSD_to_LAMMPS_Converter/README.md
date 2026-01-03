# XSD to LAMMPS Data File Converter

Converts Materials Studio XSD files to LAMMPS data file format.

## Features

- ✅ Parses Materials Studio XSD format (v20.1+)
- ✅ Extracts atoms with coordinates, charges, and element types
- ✅ Extracts bonds with proper atom connectivity
- ✅ **Differentiates bond types** based on:
  - Atom elements (C-C, C-H, C-O, O-H, C-N, etc.)
  - Bond order (Single, Double, Aromatic)
- ✅ Handles non-sequential atom IDs
- ✅ Generates LAMMPS-compatible data files
- ✅ Batch conversion of multiple files

## Requirements

- Python 3.6+
- No external dependencies (uses only standard library)

## Usage

### Single File Conversion

```python
from xsd_to_lammps import XSDtoLAMMPS

converter = XSDtoLAMMPS('molecule.xsd')
converter.parse_xsd()
converter.write_lammps_data('molecule.data')
```

### Batch Conversion

Convert all XSD files in current directory:

```bash
python convert_all_xsd.py
```

### Command Line (One-liner)

```bash
python -c "from xsd_to_lammps import XSDtoLAMMPS; c = XSDtoLAMMPS('molecule.xsd'); c.parse_xsd(); c.write_lammps_data('molecule.data')"
```

## Output Format

The converter generates LAMMPS data files with:

```
LAMMPS data file - converted from XSD

8 atoms
7 bonds
0 angles
0 dihedrals
0 impropers

2 atom types
2 bond types
0 angle types
0 dihedral types

-16.841545 5.485442 xlo xhi
-11.587840 10.245384 ylo yhi
-11.038180 10.760274 zlo zhi

Masses

1 12.011  # C
2 1.008  # H

# Bond Types Reference:
#   Type 1: C-C (Single)
#   Type 2: C-H (Single)

Atoms

1 1 1 0.000000 -6.45019656 -0.54394656 0.00002697
2 1 1 0.000000 -4.91038533 -0.52056858 0.00006581
...

Bonds

1 1 1 2
2 2 1 3
...
```

## Bond Type Examples

| Molecule | Bond Types Generated |
|----------|---------------------|
| H₂O | H-O (Single) |
| CH₄ | C-H (Single) |
| C₂H₆ | C-C (Single), C-H (Single) |
| C₂H₄ | C-C (Double), C-H (Single) |
| Graphene Oxide | C-C (Aromatic), C-C (Single), C-H (Single), C-O (Single), C-O (Double), H-O (Single), C-N (Single), H-N (Single) |

## Examples

The `examples/` directory contains sample XSD files:

- `Single_H2O.xsd` - Water molecule (3 atoms, 2 bonds)
- `Single_CH4.xsd` - Methane (5 atoms, 4 bonds)
- `Single_C2H6.xsd` - Ethane (8 atoms, 7 bonds)
- `Single_C2H4.xsd` - Ethene with C=C double bond (6 atoms, 5 bonds)

## Verification with LAMMPS

Test the generated data file:

```bash
# Create test input file
echo "units real
atom_style full
boundary p p p
read_data molecule.data
print 'LAMMPS reading successful'" > test.lmp

# Run LAMMPS
lmp -i test.lmp
```

## XSD File Structure

The converter handles Materials Studio XSD format:

```xml
<Molecule ID="2" NumChildren="5" Name="Water">
    <Atom3d ID="3" XYZ="-2.536,0.674,0.025" 
            Components="O" Charge="-0.82" 
            ForcefieldType="o2*" Connections="4,5"/>
    <Atom3d ID="6" XYZ="-2.792,1.327,0.677" 
            Components="H" Charge="0.41"/>
    <Bond ID="4" Connects="3,6"/>
    <Atom3d ID="7" XYZ="-1.579,0.686,0.036" 
            Components="H" Charge="0.41"/>
    <Bond ID="5" Connects="3,7" Type="Double"/>
</Molecule>
```

Key points:
- `Atom3d` elements contain coordinates, element type, charge
- `Bond` elements define connectivity via `Connects` attribute
- `Type` attribute specifies bond order (Single, Double, Aromatic)
- Atom IDs are non-sequential and mapped to sequential LAMMPS indices

## Limitations

- **Angles/Dihedrals**: Not explicitly stored in XSD files; LAMMPS auto-generates from bond topology
- **Force field parameters**: Not included; add manually or use LAMMPS `pair_coeff`, `bond_coeff` commands

## License

MIT License - Feel free to use and modify.

## Author

Created for converting Materials Studio structures to LAMMPS format.
