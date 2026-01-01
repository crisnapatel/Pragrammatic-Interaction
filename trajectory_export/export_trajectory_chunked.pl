#!perl
# =============================================================================
# Export MS Trajectory to LAMMPS Dump Format (Chunked)
# =============================================================================
# Exports NVT1.xtd trajectory to multiple LAMMPS dump files for faster I/O.
# Each file contains FRAMES_PER_FILE frames.
#
# Output: trajectory_chunks/chunk_XXXX.lammpstrj
# =============================================================================

use strict;
use warnings;
use MaterialsScript qw(:all);

# Configuration
my $trajectory_file = "NVT1.xtd";
my $output_dir = "D:/MS_Projects/Krishna/aPS30_AND_aPS40_AND_aPS50_GRA_TOL_Files/Documents/Production/aPS30-3_GRA_TOL/NVT1 Forcite Dynamics";
my $chunk_dir = "$output_dir/trajectory_chunks";

# Frame sampling
my $frame_start = 1;      # Start frame
my $frame_step = 1;       # Export every Nth frame
my $FRAMES_PER_FILE = 50; # Frames per chunk file

# Create chunk directory
mkdir($chunk_dir) unless -d $chunk_dir;

# Open trajectory document
my $doc = $Documents{$trajectory_file};
die "Cannot open trajectory file: $trajectory_file\n" unless $doc;

my $trajectory = $doc->Trajectory;
my $numFrames = $trajectory->NumFrames;
print "Trajectory: $trajectory_file\n";
print "Total frames: $numFrames\n";
print "Frames per chunk: $FRAMES_PER_FILE\n";
print "Output directory: $chunk_dir\n";

# Get cell parameters
my $symmetryDef = $doc->SymmetryDefinition;
my $cellA = $symmetryDef->LengthA;
my $cellB = $symmetryDef->LengthB;
my $cellC = $symmetryDef->LengthC;
print "Cell dimensions: $cellA x $cellB x $cellC Angstrom\n";

# Get all atoms
my @all_atoms = @{$doc->UnitCell->Atoms};
my $num_atoms = scalar(@all_atoms);
print "Total atoms: $num_atoms\n";

# Build element to type mapping
my %element_type_map;
my $type_counter = 1;
foreach my $atom (@all_atoms) {
    my $elem = $atom->ElementSymbol;
    if (!exists $element_type_map{$elem}) {
        $element_type_map{$elem} = $type_counter++;
    }
}
print "Element types: ";
foreach my $elem (sort keys %element_type_map) {
    print "$elem=$element_type_map{$elem} ";
}
print "\n";

# =============================================================================
# Main loop: Export frames in chunks
# =============================================================================
print "Processing frames...\n";

my $chunk_num = 0;
my $fh;
my $frames_in_current_chunk = 0;
my $chunk_start_frame = $frame_start;

for (my $frame = $frame_start; $frame <= $numFrames; $frame += $frame_step) {
    
    # Start new chunk file if needed
    if ($frames_in_current_chunk == 0) {
        my $chunk_file = sprintf("%s/chunk_%04d.lammpstrj", $chunk_dir, $chunk_num);
        open($fh, '>', $chunk_file) or die "Cannot open $chunk_file: $!";
        $chunk_start_frame = $frame;
    }
    
    $trajectory->CurrentFrame = $frame;
    
    # Write LAMMPS dump frame header
    print $fh "ITEM: TIMESTEP\n";
    print $fh "$frame\n";
    
    print $fh "ITEM: NUMBER OF ATOMS\n";
    print $fh "$num_atoms\n";
    
    print $fh "ITEM: BOX BOUNDS pp pp pp\n";
    print $fh "0.0 $cellA\n";
    print $fh "0.0 $cellB\n";
    print $fh "0.0 $cellC\n";
    
    print $fh "ITEM: ATOMS id type element mass q x y z\n";
    
    # Write atom data
    my $atom_id = 1;
    foreach my $atom (@all_atoms) {
        my $elem = $atom->ElementSymbol;
        my $type = $element_type_map{$elem};
        my $mass = $atom->Mass;
        my $charge = $atom->Charge || 0;
        my $x = $atom->X;
        my $y = $atom->Y;
        my $z = $atom->Z;
        
        printf $fh "%d %d %s %.4f %.6f %.6f %.6f %.6f\n",
            $atom_id, $type, $elem, $mass, $charge, $x, $y, $z;
        
        $atom_id++;
    }
    
    $frames_in_current_chunk++;
    
    # Close chunk if full
    if ($frames_in_current_chunk >= $FRAMES_PER_FILE) {
        close($fh);
        print "  Chunk $chunk_num: frames $chunk_start_frame-$frame\n";
        $chunk_num++;
        $frames_in_current_chunk = 0;
    }
}

# Close final chunk if not already closed
if ($frames_in_current_chunk > 0) {
    close($fh);
    print "  Chunk $chunk_num: frames $chunk_start_frame-$numFrames ($frames_in_current_chunk frames)\n";
    $chunk_num++;
}

print "\n=== Export Complete ===\n";
print "Total chunks: $chunk_num\n";
print "Chunk directory: $chunk_dir\n";
print "\nNext: Run merge_chunks.py or directly load chunks in Python/OVITO\n";
