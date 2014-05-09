#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Take an input PDB and optimize the box size and cutoff used to create sphere
models. The optimization is run to recreate the volume predicted from sluv2.
The box size determines the size of the grid used to create the sphere models
from a PDB. The cutoff is used in the creation of the hydration layer.

These parameters can then be put into the input file (yaml) for the
analyse_pdb_models.py script to run on a series of models.
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

import sct


def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(description = 'Optimize the box side used to create dry sphere models and cutoff used to creat a hydration layer that reproduce expected volumes from an input PDB.\n')

    parser.add_argument('-i','--input_pdb', nargs='?', type=str,
        help = 'Path to the input PDB file', required=True)

    parser.add_argument('-s','--input_seq', nargs='?', type=str,
        help = 'Path to input sequence file if different from PDB', 
        default = None)
        
    parser.add_argument('-c','--cutoff', nargs='?', type=int,
                        default=4,
                        help = 'Cutoff number of atoms to create sphere in grid box')
                        
    parser.add_argument('-f','--input_format', choices = ['fas','yml'],
        help = 'Input file format (fasta or sluv yaml)', default = None)

    parser.add_argument('-o','--output_file', nargs='?', type=str,
        default='optimized_params.txt', help = 'Path to output file')

    parser.add_argument('-p','--parameter_file', nargs='?', type=str,
        help = 'Path to a file containing input parameters', default = None)

    parser.add_argument('-m','--min_max', nargs=2, type=float,
        default = [2.0, 11.0], help = 'Minimum and maximum values of box side')

    parser.add_argument('-w','--wcut_min_max', nargs=2, type=int,
        default = [2, 20], help = 'Minimum and maximum values of water cutoff')

    return parser.parse_args()

def get_box_opt_input(pdb_filename, seq_filename, seq_type):
    """
    Get the target volume and atomic coordinates for box_side optimization. If 
    specified use a sequence other than that of the input PDB providing the 
    coordinates.
    @type  pdb_filename: string
    @param pdb_filename: Path to PDB file
    @type  seq_filename: string
    @param seq_filename: Path to sequence file (either fasta or YAML)
    @type  file_type:    string
    @param file_type:    Is the input a fasta ('fas') or YAML ('yml') file
    @rtype:              float, list
    @return:             1. Target volume from sequence

                         2. A list containing lists of x, y & z coordinates 
                         (3 * floats)     
    """

    # Read in the residues frequencies (to calculate target volume) and
    # atomic coordinates from input PDB
    res_freq, atom_coords = sct.pdb.read_pdb_atom_data(pdb_filename)

    # If an additional sequence file was provided then use this as the sequence
    # to optimize 
    if seq_type != None:
        res_freq = sct.seq.seq_file_to_freq(seq_filename, seq_type)

    dry_volume = sct.seq.sum_volume(sct.seq.all_residues, 
                                    res_freq, 'perkins1986a')
    # Calculate the expected volume of the hydration layer
    # Hydration is estimated to be 0.3 g water per 1 g protein
    # Bound water volume (< bulk volume) is given as a sluv parameter    
    wet_volume = sct.seq.calc_hydration_volume(res_freq) + dry_volume 
    
    return dry_volume, wet_volume, atom_coords

def main():
    
    args = parse_arguments()

    if args.parameter_file != None:

        print "WARNING: A SCT parameter file was specified, so the modelling parameters from the command line flags will be ignored!"        

        # Read in parameters and check we have those we need
        param, err = sct.param.read_parameter_file(args.parameter_file)
        
        if err != None:
            sct.param.output_error(err, args.parameter_file)        
        
        sct.param.check_parameters(param, ['curve','sphere','rg','rxs1','rfac'])

        if err != None:
            sct.param.output_error(err, args.parameter_file)

        cutoff = param['sphere']['cutoff']
    else:
        cutoff = args.cutoff
    
    # Get target volumes for the dry and hydrated protein based on the sequence 
    # from input PDB or separate sequence file if provided and coordinates from 
    # the PDB
    dry_volume, wet_volume, atom_coords = get_box_opt_input(args.input_pdb, 
                                                            args.input_seq, 
                                                            args.input_format)

    # Optimize box side for sphere models
    box_side, err = sct.sphere.optimize_box_side(cutoff, atom_coords, 
                                                 dry_volume,
                                                 args.min_max[0],
                                                 args.min_max[1],
                                                 0.01)

    dry_spheres, x_axis, y_axis, z_axis = sct.sphere.create_sphere_model(
                                                                    atom_coords,
                                                                    cutoff,
                                                                    box_side)
    model_dry_vol = len(dry_spheres) * box_side**3

    # Optimize the water cutoff used to filter out excessive/overlapping 
    # hydration spheres
    # Use 26 additional spheres in the trial hydration
    # 10 - 14 is a reasonable window over which to refine the cutoff 
    # for proteins, this didn't work for linear carbohydrates
    hydr_cutoff, err = sct.sphere.optimize_watercut(box_side,
                                                    dry_spheres,
                                                    27,
                                                    wet_volume,
                                                    args.wcut_min_max[0],
                                                    args.wcut_min_max[1],
                                                    0.01)
                                                    
    wet_spheres = sct.sphere.hydrate_sphere_model(dry_spheres,
                                            27,
                                            box_side,
                                            hydr_cutoff,
                                            xaxis = x_axis,
                                            yaxis = y_axis,
                                            zaxis = z_axis)
    model_wet_vol = len(wet_spheres) * box_side**3

    out_file = open(args.output_file,'w')
    out_file.write("Sphere model box side: {0:7.4f}\n".format(box_side))
    out_file.write("Target dry volume: {0:7.4f}\n".format(dry_volume))
    out_file.write("Model dry volume: {0:7.4f}\n".format(model_dry_vol))
    for ii in range(args.wcut_min_max[0], args.wcut_min_max[1] + 1):
        wet_spheres = sct.sphere.hydrate_sphere_model(dry_spheres,
                                            27,
                                            box_side,
                                            ii,
                                            xaxis = x_axis,
                                            yaxis = y_axis,
                                            zaxis = z_axis)
        check_wet_vol = len(wet_spheres) * box_side**3        
        out_file.write("Hydrated volume, water cutoff = {0:d}: {1:7.4f}\n".format(ii,check_wet_vol))
        
    out_file.write("Suggested Hydration cutoff [Experimental]: {0:d}\n".format(hydr_cutoff))
    out_file.write("Target hydrated volume: {0:7.4f}\n".format(wet_volume))
    out_file.write("Model hydrated volume: {0:7.4f}\n".format(model_wet_vol))


if __name__ == "__main__":
    main()
