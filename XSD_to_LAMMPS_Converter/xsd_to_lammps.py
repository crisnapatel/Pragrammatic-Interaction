"""
XSD to LAMMPS Data File Converter
Converts Materials Studio XSD files to LAMMPS data format.
"""

import xml.etree.ElementTree as ET
import numpy as np
from collections import defaultdict, OrderedDict
import re


class XSDtoLAMMPS:
    def __init__(self, xsd_file):
        self.xsd_file = xsd_file
        self.atoms = {}
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.atom_types = {}
        self.bond_types = {}
        self.angle_types = {}
        self.dihedral_types = {}
        self.atom_type_counter = 0
        self.bond_type_counter = 0
        self.angle_type_counter = 0
        self.dihedral_type_counter = 0
        self.box_bounds = {'xlo': 0, 'xhi': 0, 'ylo': 0, 'yhi': 0, 'zlo': 0, 'zhi': 0}
        self.xsd_id_to_lammps_idx = {}  # Map XSD atom IDs to LAMMPS sequential indices
        
    def parse_xsd(self):
        """Parse XSD file and extract atomic information"""
        tree = ET.parse(self.xsd_file)
        root = tree.getroot()
        
        # Find all Atom3d elements
        ns = {'xsd': 'http://www.accelrys.com/xsd'}  # Default namespace if present
        
        # Try with and without namespace
        atoms_list = root.findall('.//Atom3d')
        if not atoms_list:
            atoms_list = root.findall('.//{http://www.accelrys.com/xsd}Atom3d')
        
        print(f"Found {len(atoms_list)} atoms")
        
        # Extract atoms
        for idx, atom_elem in enumerate(atoms_list, 1):
            atom_id = atom_elem.get('ID')
            # Store mapping from XSD ID to LAMMPS index
            self.xsd_id_to_lammps_idx[int(atom_id)] = idx
            atom_name = atom_elem.get('Name')
            forcefield_type = atom_elem.get('ForcefieldType', 'unknown')
            charge = float(atom_elem.get('Charge', 0.0))
            xyz_str = atom_elem.get('XYZ', '0,0,0')
            
            try:
                x, y, z = map(float, xyz_str.split(','))
            except:
                x, y, z = 0.0, 0.0, 0.0
            
            # Get element from atom components
            components_elem = atom_elem.find('.//Components') or atom_elem.find('./{http://www.accelrys.com/xsd}Components')
            element = 'X'
            if components_elem is not None:
                element = components_elem.text if components_elem.text else 'X'
            else:
                components_attr = atom_elem.get('Components', 'X')
                element = components_attr
            
            self.atoms[idx] = {
                'id': atom_id,
                'name': atom_name,
                'element': element,
                'forcefield_type': forcefield_type,
                'charge': charge,
                'x': x,
                'y': y,
                'z': z,
                'atom_type': None
            }
            
            # Track unique atom types
            # Use forcefield_type if available, otherwise use element
            if forcefield_type == 'unknown' or not forcefield_type:
                type_key = element  # Use element symbol as type key
            else:
                type_key = forcefield_type
            
            if type_key not in self.atom_types:
                self.atom_type_counter += 1
                self.atom_types[type_key] = self.atom_type_counter
            
            self.atoms[idx]['atom_type'] = self.atom_types[type_key]
            self.atoms[idx]['type_key'] = type_key  # Store for reference
        
        # Extract bonds from Bond elements
        self._extract_bonds_from_bond_elements(root)
        
        # Calculate box bounds
        self._calculate_box_bounds()
        
        return len(self.atoms)
    
    def _extract_bonds_from_bond_elements(self, root):
        """
        Extract bond information from explicit Bond elements.
        Bond elements have format: <Bond ID="4" Connects="3,6" Type="Double"/>
        where Connects specifies the atom IDs being connected.
        
        Bond types are determined by:
        1. Element types of both atoms (e.g., C-C, C-H, O-H)
        2. Bond order (Single, Double, Triple, Aromatic)
        """
        # Find all Bond elements
        bonds_list = root.findall('.//Bond')
        if not bonds_list:
            bonds_list = root.findall('.//{http://www.accelrys.com/xsd}Bond')
        
        bond_id = 1
        processed_bonds = set()
        
        print(f"Found {len(bonds_list)} bonds in Bond elements")
        
        for bond_elem in bonds_list:
            connects = bond_elem.get('Connects', '')
            bond_order = bond_elem.get('Type', 'Single')  # Default to Single if not specified
            
            if connects:
                try:
                    # Parse Connects attribute: "atom_id1,atom_id2"
                    parts = [int(p.strip()) for p in connects.split(',')]
                    if len(parts) == 2:
                        atom1_xsd_id, atom2_xsd_id = parts
                        
                        # Map XSD atom IDs to LAMMPS sequential indices
                        atom1_idx = self.xsd_id_to_lammps_idx.get(atom1_xsd_id)
                        atom2_idx = self.xsd_id_to_lammps_idx.get(atom2_xsd_id)
                        
                        # Skip if atoms not found in mapping
                        if atom1_idx is None or atom2_idx is None:
                            print(f"Warning: Bond references unknown atom IDs {atom1_xsd_id},{atom2_xsd_id}")
                            continue
                        
                        # Create bond tuple (ensure consistent ordering)
                        bond_tuple = tuple(sorted([atom1_idx, atom2_idx]))
                        
                        if bond_tuple not in processed_bonds:
                            processed_bonds.add(bond_tuple)
                            
                            # Get element types for both atoms
                            elem1 = self.atoms[atom1_idx]['element']
                            elem2 = self.atoms[atom2_idx]['element']
                            
                            # Create bond type key based on elements and bond order
                            # Sort element names for consistency (C-H same as H-C)
                            elem_pair = tuple(sorted([elem1, elem2]))
                            bond_type_key = (elem_pair[0], elem_pair[1], bond_order)
                            
                            # Assign bond type number
                            if bond_type_key not in self.bond_types:
                                self.bond_type_counter += 1
                                self.bond_types[bond_type_key] = self.bond_type_counter
                            
                            self.bonds.append({
                                'id': bond_id,
                                'type': self.bond_types[bond_type_key],
                                'atom1': bond_tuple[0],
                                'atom2': bond_tuple[1],
                                'order': bond_order,
                                'elements': bond_type_key
                            })
                            bond_id += 1
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse bond: {connects}, Error: {e}")
                    continue
        
        # Print bond type summary
        if self.bond_types:
            print(f"Bond types found:")
            for bond_key, bond_num in sorted(self.bond_types.items(), key=lambda x: x[1]):
                elem1, elem2, order = bond_key
                print(f"  Type {bond_num}: {elem1}-{elem2} ({order})")
    
    def _calculate_box_bounds(self):
        """Calculate simulation box bounds from atomic coordinates"""
        if not self.atoms:
            return
        
        x_coords = [atom['x'] for atom in self.atoms.values()]
        y_coords = [atom['y'] for atom in self.atoms.values()]
        z_coords = [atom['z'] for atom in self.atoms.values()]
        
        # Add 10 Angstrom padding
        padding = 10.0
        
        self.box_bounds['xlo'] = min(x_coords) - padding
        self.box_bounds['xhi'] = max(x_coords) + padding
        self.box_bounds['ylo'] = min(y_coords) - padding
        self.box_bounds['yhi'] = max(y_coords) + padding
        self.box_bounds['zlo'] = min(z_coords) - padding
        self.box_bounds['zhi'] = max(z_coords) + padding
    
    def write_lammps_data(self, output_file):
        """Write LAMMPS data file"""
        with open(output_file, 'w') as f:
            # Header
            f.write("LAMMPS data file - converted from XSD\n")
            f.write("\n")
            
            # Counts
            f.write(f"{len(self.atoms)} atoms\n")
            f.write(f"{len(self.bonds)} bonds\n")
            f.write(f"{len(self.angles)} angles\n")
            f.write(f"{len(self.dihedrals)} dihedrals\n")
            f.write("0 impropers\n")
            f.write("\n")
            
            # Types
            f.write(f"{len(self.atom_types)} atom types\n")
            f.write(f"{len(self.bond_types)} bond types\n")
            f.write(f"{len(self.angle_types)} angle types\n")
            f.write(f"{len(self.dihedral_types)} dihedral types\n")
            f.write("\n")
            
            # Box bounds
            f.write(f"{self.box_bounds['xlo']:.6f} {self.box_bounds['xhi']:.6f} xlo xhi\n")
            f.write(f"{self.box_bounds['ylo']:.6f} {self.box_bounds['yhi']:.6f} ylo yhi\n")
            f.write(f"{self.box_bounds['zlo']:.6f} {self.box_bounds['zhi']:.6f} zlo zhi\n")
            f.write("\n")
            
            # Masses section with atom type comments
            f.write("Masses\n")
            f.write("\n")
            
            atomic_masses = {
                'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
                'S': 32.065, 'P': 30.974, 'Cl': 35.45, 'F': 18.998,
                'Br': 79.904, 'I': 126.904, 'Na': 22.990, 'K': 39.098,
                'Ca': 40.078, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085
            }
            
            for ff_type in sorted(self.atom_types.keys(), key=lambda x: self.atom_types[x]):
                atom_type_id = self.atom_types[ff_type]
                # Extract element from forcefield type or use directly if it's an element
                element = re.sub(r'[^a-zA-Z]', '', ff_type)
                if not element:
                    element = 'X'
                # Get first letter capitalized for mass lookup
                elem_key = element[0].upper() if len(element) == 1 else element[:2].capitalize()
                if elem_key not in atomic_masses:
                    elem_key = element[0].upper()
                mass = atomic_masses.get(elem_key, 12.0)
                f.write(f"{atom_type_id} {mass:.3f}  # {ff_type}\n")
            
            f.write("\n")
            
            # Bond Coeffs section (commented, for user to fill in)
            if self.bond_types:
                f.write("# Bond Types Reference:\n")
                for bond_key, bond_num in sorted(self.bond_types.items(), key=lambda x: x[1]):
                    elem1, elem2, order = bond_key
                    f.write(f"#   Type {bond_num}: {elem1}-{elem2} ({order})\n")
                f.write("\n")
            
            # Atoms section
            f.write("Atoms\n")
            f.write("\n")
            
            for atom_idx in sorted(self.atoms.keys()):
                atom = self.atoms[atom_idx]
                f.write(f"{atom_idx} 1 {atom['atom_type']} {atom['charge']:.6f} "
                       f"{atom['x']:.8f} {atom['y']:.8f} {atom['z']:.8f}\n")
            
            f.write("\n")
            
            # Bonds section
            if self.bonds:
                f.write("Bonds\n")
                f.write("\n")
                for bond in self.bonds:
                    f.write(f"{bond['id']} {bond['type']} {bond['atom1']} {bond['atom2']}\n")
                f.write("\n")
        
        print(f"LAMMPS data file written to {output_file}")
