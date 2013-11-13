#!/usr/bin/env python

"""Convert atomic protein coordinates in PDB format into corresponding
Debye spheres for small angle scattering calculations"""

import argparse
import sys
from bisect import bisect_right
import numpy as np
import sluv2
from sct_seq import *

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description= 'Convert atomic pdb to sphere model\n')

    parser.add_argument('-i','--input_filename', nargs='?', type=str,
        dest='input_filename', help = 'Path to the input PDB file',
        required=True)

    # Modes for output:
    # info - info on sphere model to screen only
    # archive - info and sphere model to file
    # model - model to file (no info)
    # project - model to file and info to screen
    parser.add_argument('-m','--mode', choices = ['info', 'archive', 'model', 'project'],
        default = 'project', help = 'Type of output to be produced')

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
        dest='output_filename', default=None,
        help = 'Path to the output file')

    parser.add_argument('-c','--cutoff', nargs='?', type=int,
        default=4, help = 'Cutoff number of atoms to create sphere in grid box')

    parser.add_argument('-b','--box_side', nargs='?', type=float,
        default=5.0, help = 'Side length of grid boxes')

    args = parser.parse_args()

    if args.mode != 'info':
        if args.output_filename == None:
            print "For mode {0:s} an output file must be specified!".format(args.mode)
            sys.exit(1)

    return args

def read_pdb_atom_data (filename):
    """Read PDB file and return residue ferquencies and atom coordinates"""

    # Initialize a dictionary of all recognized residues (frequency = 0)
    res_freq = residue_freq_dict()

    # Track the last residues seen
    last_res_no = 0

    atom_coords = []

    # Read in lines of the pdb using sct_seq function
    # ATOM and HETATM records are interpretted - others ignored
    with open(filename) as f:
        for line in f:
            data = pdb_res_line_parse(line)

            if len(data) != 0:
                atom_coords.append(data['coords'])
                # If residue number has changed increment res_id tally
                if data['res_no'] != last_res_no:
                    res_freq[data['res_id']] += 1
                    last_res_no = data['res_no']

    return res_freq, atom_coords

def create_grid(atom_coords, x_axis, y_axis, z_axis):
    """Count atoms in grid boxes defined by passed axes intervals"""

    # Initialize grid
    # dimension -1 as looking at intervals rather than axes points
    grid = np.zeros([len(x_axis) - 1, len(y_axis) - 1, len(z_axis) - 1])

    # Use a binary search to find which interval each atom coordinate fits in
    # then increment count in appropriate grid box
    for atom_coord in atom_coords:
        x = bisect_right(x_axis, atom_coord[0]) - 1
        y = bisect_right(y_axis, atom_coord[1]) - 1
        z = bisect_right(z_axis, atom_coord[2]) - 1
        grid[x, y, z] += 1

    return grid

def define_axes(coords, box_side):
    """Create intervals of length box_side along x, y and z axes which contain coords"""

    # Create limits of the axes which contain all coords
    coord_mins = np.floor(np.amin(coords, 0))
    coord_maxs = np.ceil(np.amax(coords, 0))

    # Create float valued arrays of intervals
    x_axis = np.arange(coord_mins[0], coord_maxs[0] + box_side, box_side)
    y_axis = np.arange(coord_mins[1], coord_maxs[1] + box_side, box_side)
    z_axis = np.arange(coord_mins[2], coord_maxs[2] + box_side, box_side)

    return x_axis, y_axis, z_axis

def grid_to_spheres(grid, cutoff, x_axis, y_axis, z_axis):

    # Create a grid containing the number of atoms in each cube defined
    # by the axes created above
    grid = create_grid(atom_coords, x_axis, y_axis, z_axis)

    # Find the locations in the grid containing more than cutoff atoms
    # A sphere will be placed at the centre of each of these cubes
    full = np.where(grid > cutoff)

    sphere_coords = []

    for ndx in xrange(0, no_spheres):
        x = x_axis[full[0][ndx]] + radius
        y = y_axis[full[1][ndx]] + radius
        z = z_axis[full[2][ndx]] + radius
        sphere_coords.append([x,y,z])

    return sphere_coords

def write_spheres(coord_list, radius, out):

    for coords in coord_list:
        write_sphere_line(coords[0], coords[1], coords[2], radius, out)

def write_sphere_line(x, y, z, radius, out):

    out.write("{0:10.2f}{1:10.2f}{2:10.2f}{3:10.2f}\n".format(x, y, z, radius))

def create_sphere_model(atom_coords, cutoff, box_side):

    # Conversion to numpy array speeds up operations later
    atom_coords = np.array(atom_coords)

    # Get axes for the grid used to create sphere models
    # The grid we create contains all atoms and is divided into cubes with
    # dimension box_side
    x_axis, y_axis, z_axis = define_axes(atom_coords, box_side)

    # Create a grid containing the number of atoms in each cube defined
    # by the axes created above
    grid = create_grid(atom_coords, x_axis, y_axis, z_axis)

    # Convert grid to spheres
    # A sphere centre of a box is created if > cutoff atoms are within it
    sphere_coords = grid_to_spheres(grid, cutoff, x_axis, y_axis, z_axis)

    return sphere_coords

def main():

    # Process command line inputs
    args = parse_arguments()

    # Read in the residues frequencies (to calculate target volume) and
    # atomic coordinates from input PDB
    res_freq, atom_coords = read_pdb_atom_data(args.input_filename)

    # Bin atomic coordinates to give sphere coordinates
    sphere_coords = create_sphere_model(atom_coords, args.cutoff, args.box_side)

    # Set the radius for each sphere
    radius = box_side / 2.0

    # Output the sphere coordinates and radii
    if args.mode != 'info':

        output = open(args.output_filename,'w')
        write_spheres(sphere_coords, radius, output)

    # Depending on output mode chosen print out information on the models
    if args.mode != 'model':

        if args.mode == 'info':
            output = sys.stdout

        no_atoms = len(atom_coords)
        volume = sluv2.sum_volume(sluv2.all_residues, res_freq, 'chothia1975')
        no_res = sluv2.sum_res_no(sluv2.all_residues, res_freq)
        no_spheres = len(sphere_coords)
        volume_spheres = no_spheres * box_side**3
        no_x_box = len(x_axis)
        no_y_box = len(y_axis)
        no_z_box = len(z_axis)

        output.write("pdb_to_sphere: version 0.5 - 05 November 2013\n")
        output.write("No. Of Atoms: {0:d}\n".format(no_atoms))
        output.write("Total Of Amino Acid Residues {0:d}\n".format(no_res))
        output.write("One Side Of The Box: {0:f}\n".format(box_side))
        output.write("Cutoff For Creating A Ball: {0:d}\n".format(cutoff))
        output.write("No Of Balls: {0:d}\n".format(no_spheres))
        output.write("Volume Of Cubes: {0:f}\n".format(volume_spheres))
        output.write("Volume of Protein: {0:f}\n".format(volume))
        output.write("Max And Min X: {0:f} {1:f}\n".format(x_axis[0],x_axis[no_x_box - 1]))
        output.write("Max And Min Y: {0:f} {1:f}\n".format(y_axis[0],y_axis[no_y_box - 1]))
        output.write("Max And Min Z: {0:f} {1:f}\n".format(z_axis[0],z_axis[no_z_box - 1]))
        output.write("No. grid points in X, Y, Z: {0:d} {1:d} {2:d}\n".format(no_x_box, no_y_box, no_z_box))

    if args.mode != 'info':
        output.close()

if __name__ == "__main__":
    main()
