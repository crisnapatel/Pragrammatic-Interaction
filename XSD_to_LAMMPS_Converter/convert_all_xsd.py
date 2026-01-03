#!/usr/bin/env python3
"""
XSD to LAMMPS Data File Converter - Batch Script
Converts all XSD files in current directory to LAMMPS data format.
"""

from xsd_to_lammps import XSDtoLAMMPS
import glob
import os

def main():
    # Find all XSD files
    xsd_files = sorted(glob.glob('*.xsd'))
    
    if not xsd_files:
        print("No XSD files found in current directory")
        return
    
    print("=" * 70)
    print("XSD to LAMMPS Converter")
    print("=" * 70)
    print(f"\nFound {len(xsd_files)} XSD file(s):\n")
    
    for f in xsd_files:
        print(f"  - {f}")
    
    results = {}
    
    for xsd_file in xsd_files:
        print("\n" + "=" * 70)
        print(f"Converting: {xsd_file}")
        print("=" * 70)
        
        try:
            converter = XSDtoLAMMPS(xsd_file)
            num_atoms = converter.parse_xsd()
            
            # Generate output filename
            output_file = xsd_file.replace('.xsd', '.data')
            
            converter.write_lammps_data(output_file)
            
            print(f"\n✓ SUCCESS")
            print(f"  {num_atoms} atoms")
            print(f"  {len(converter.bonds)} bonds")
            print(f"  {len(converter.atom_types)} atom types")
            print(f"  Output: {output_file}")
            
            results[xsd_file] = True
            
        except Exception as e:
            print(f"\n✗ FAILED")
            print(f"  Error: {str(e)}")
            import traceback
            traceback.print_exc()
            results[xsd_file] = False
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    for xsd_file, success in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"  {status}: {xsd_file}")
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    print(f"\nCompleted: {passed}/{total} conversions successful\n")

if __name__ == '__main__':
    main()
