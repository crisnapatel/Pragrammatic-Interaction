# Pragrammatic Interaction

A collection of Python, Perl, and Bash scripts for processing molecular dynamics simulation data and trajectories. Each objective has its own directory with related scripts and utilities.

## Project Structure

### [`lammps_log_parsing/`](lammps_log_parsing/)
Tools for parsing and analyzing LAMMPS simulation log files.

- **`lammps_log_parser.ipynb`** - Jupyter notebook that extracts thermodynamic properties from LAMMPS log files
  - Parses thermo data (temperature, pressure, stress components, etc.)
  - Returns data as pandas DataFrame for analysis
  - Handles multiple log files and property filtering

### [`trajectory_export/`](trajectory_export/)
Scripts for converting Materials Studio trajectory files to LAMMPS-compatible formats.

- **`export_trajectory_chunked.pl`** - Perl script that converts Materials Studio (2020.1v) forcite trajectory files to LAMMPS dump format
  - Exports trajectories for **all atoms**
  - Output compatible with OVITO for visualization
  
- **`export_sets_trajectory.pl`** - Perl script that exports trajectories for selected atom sets
  - Exports only atoms from predefined `Sets` (see Materials Studio documentation)
  - Output compatible with OVITO for visualization

### [`XSD_to_LAMMPS_Converter/`](XSD_to_LAMMPS_Converter/)
Tools for converting Materials Studio XSD files to LAMMPS data file format.

- **`xsd_to_lammps.py`** - Main converter class for parsing XSD files and generating LAMMPS data files
  - Extracts atoms with coordinates, charges, and element types
  - Extracts bonds with proper atom connectivity
  - Differentiates bond types based on atom elements and bond order
  - Handles non-sequential atom IDs
  
- **`convert_all_xsd.py`** - Batch conversion utility
  - Converts all XSD files in a directory to LAMMPS format
  - Automated processing for high-throughput workflows

## Requirements

- Python 3.x with pandas
- Perl 5.x
- LAMMPS (for generating log files)
- OVITO (optional, for trajectory visualization)

## Usage

Refer to the README files in each subdirectory for specific usage instructions and examples.
