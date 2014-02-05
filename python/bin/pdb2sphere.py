#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert atomic protein coordinates in PDB format into corresponding Debye
spheres for small angle scattering calculations
"""

# Copyright 2014 University College London

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import sys

import sct

def parse_arguments():
    """
    Parse command line arguments and ensure correct combinations present
    """

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
    parser.add_argument('-m','--mode',
                        choices = ['info', 'archive', 'model', 'project'],
                        default = 'project',
                        help = 'Type of output to be produced')

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
                        dest='output_filename', default=None,
                        help = 'Path to the output file')

    parser.add_argument('-c','--cutoff', nargs='?', type=int,
                        default=4,
                        help = 'Cutoff number of atoms to create sphere in grid box')

    parser.add_argument('-b','--box_side', nargs='?', type=float,
                        default=5.0,
                        help = 'Side length of grid boxes')

    args = parser.parse_args()

    if args.mode != 'info':
        if args.output_filename == None:
            print "For mode {0:s} an output file must be specified!".format(args.mode)
            sys.exit(1)

    return args

def main():

    # Process command line inputs
    args = parse_arguments()

    # Read in the residues frequencies (to calculate target volume) and
    # atomic coordinates from input PDB
    res_freq, atom_coords = sct.pdb.read_pdb_atom_data(args.input_filename)

    # Bin atomic coordinates to give sphere coordinates
    sphere_coords, x_axis, y_axis, z_axis = sct.sphere.create_sphere_model(atom_coords,
                                                                           args.cutoff,
                                                                           args.box_side)

    # Set the radius for each sphere
    radius = args.box_side / 2.0

    # Output the sphere coordinates and radii
    if args.mode != 'info':

        output = open(args.output_filename,'w')
        sct.sphere.write_spheres(sphere_coords, radius, output)

    # Depending on output mode chosen print out information on the models
    if args.mode != 'model':

        if args.mode == 'info':
            output = sys.stdout

        no_atoms = len(atom_coords)
        volume = sct.seq.sum_volume(sct.seq.all_residues, res_freq, 'perkins1986b')
        no_res = sct.seq.sum_res_no(sct.seq.all_residues, res_freq)
        no_spheres = len(sphere_coords)
        volume_spheres = no_spheres * args.box_side**3
        no_x_box = len(x_axis)
        no_y_box = len(y_axis)
        no_z_box = len(z_axis)

        output.write("pdb_to_sphere: version 0.5 - 05 November 2013\n")
        output.write("No. Of Atoms: {0:d}\n".format(no_atoms))
        output.write("Total Of Amino Acid Residues {0:d}\n".format(no_res))
        output.write("One Side Of The Box: {0:f}\n".format(args.box_side))
        output.write("Cutoff For Creating A Ball: {0:d}\n".format(args.cutoff))
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
