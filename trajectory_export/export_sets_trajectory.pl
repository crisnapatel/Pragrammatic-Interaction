#!/usr/bin/perl
# export_sets_trajectory.pl
# Export trajectory coordinates for selected sets (PS1, PS2, PS3, GRA) in LAMMPS dump format
# This excludes toluene solvent for efficiency

use strict;
use warnings;
use MaterialsScript qw(:all);

# Configuration
my $xtd_file = "NVT1.xtd";
my $output_dir = "D:/MS_Projects/Krishna/aPS30_AND_aPS40_AND_aPS50_GRA_TOL_Files/Documents/Production/aPS30-3_GRA_TOL/NVT1 Forcite Dynamics";
my $output_file = "$output_dir/trajectory_sets.lammpstrj";
my $start_frame = 1;
my $end_frame = 10000000;  # Set to actual last frame
my $frame_step = 10;    # Export every 10th frame

# Sets to export (no toluene)
my @set_names = ("PS1", "PS2", "PS3", "GRA");

# Open trajectory
my $doc = $Documents{$xtd_file};
die "Cannot open $xtd_file\n" unless $doc;

my $trajectory = $doc->Trajectory;
die "No trajectory found\n" unless $trajectory;

my $num_frames = $trajectory->NumFrames;
print "Trajectory: $xtd_file\n";
print "Total frames: $num_frames\n";

# Adjust end_frame if needed
$end_frame = $num_frames if $end_frame > $num_frames;

# Get cell dimensions (assuming constant)
$trajectory->CurrentFrame = 1;
my $lattice = $doc->Lattice3D;
my ($lx, $ly, $lz) = ($lattice->LengthA, $lattice->LengthB, $lattice->LengthC);
print "Cell dimensions: $lx x $ly x $lz Angstrom\n";

# Collect atoms from all sets and assign mol_id
my @atoms_info = ();
my $mol_id = 0;

foreach my $set_name (@set_names) {
    $mol_id++;
    my $set = $doc->UnitCell->Sets($set_name);
    if ($set) {
        my $atoms = $set->Atoms;
        my $count = $atoms->Count;
        print "Set $set_name: $count atoms (mol_id=$mol_id)\n";
        
        for (my $i = 0; $i < $count; $i++) {
            my $atom = $atoms->Item($i);
            push @atoms_info, {
                atom => $atom,
                atom_id => $i + 1,
                mol_id => $mol_id,
                set_name => $set_name,
                element => $atom->ElementSymbol,
            };
        }
    } else {
        print "Warning: Set $set_name not found\n";
    }
}

my $total_atoms = scalar(@atoms_info);
print "Total atoms to export: $total_atoms\n";

# Open output file
open(my $fh, ">", $output_file) or die "Cannot open $output_file: $!\n";

# Process frames
my $frames_processed = 0;
for (my $frame = $start_frame; $frame <= $end_frame; $frame += $frame_step) {
    $trajectory->CurrentFrame = $frame;
    
    # Write LAMMPS dump header
    print $fh "ITEM: TIMESTEP\n";
    print $fh "$frame\n";
    print $fh "ITEM: NUMBER OF ATOMS\n";
    print $fh "$total_atoms\n";
    print $fh "ITEM: BOX BOUNDS pp pp pp\n";
    print $fh "0.0 $lx\n";
    print $fh "0.0 $ly\n";
    print $fh "0.0 $lz\n";
    print $fh "ITEM: ATOMS id mol type element x y z\n";
    
    # Write atom data
    my $global_id = 0;
    foreach my $info (@atoms_info) {
        $global_id++;
        my $atom = $info->{atom};
        my $x = $atom->X;
        my $y = $atom->Y;
        my $z = $atom->Z;
        
        printf $fh "%d %d %s %s %.6f %.6f %.6f\n",
            $global_id,
            $info->{mol_id},
            $info->{set_name},
            $info->{element},
            $x, $y, $z;
    }
    
    $frames_processed++;
    if ($frames_processed % 100 == 0) {
        print "Processed $frames_processed frames (current: $frame)\n";
    }
}

close($fh);

print "\n=== Export Complete ===\n";
print "Frames processed: $frames_processed\n";
print "Output file: $output_file\n";
print "Atoms per frame: $total_atoms\n";
