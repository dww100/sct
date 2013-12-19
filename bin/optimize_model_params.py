#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Take an input PDB and optimize the box size and cutoff used to create sphere
models. The optimization is run to recreate the volume predicted from sluv2.
The box size determines the size of the grid used to create the sphere models
from a PDB. The cutoff is used in the creation of the hydration layer.

These parameters can then be put into the input file (yaml) for the
analyse_pdb_models.py script to run on a series of models.
"""

import argparse
import pdb2sphere as p2s
import hydrate_spheres as hydrate
import get_boxside
import sluv2
import yaml

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(description= 'Optimize the box side used to create dry sphere models and cutoff used to creat a hydration layer that reproduce expected volumes from an input PDB.\n')

    parser.add_argument('-i','--input_pdb', nargs='?', type=str,
        help = 'Path to the input PDB file', required=True)

    parser.add_argument('-o','--output_file', nargs='?', type=str,
        default='.', help = 'Path to output file')

    parser.add_argument('-p','--parameter_file', nargs='?', type=str,
        help = 'Path to a file containing input parameters', required=True)

    parser.add_argument('-m','--min_max', nargs=2, type=float,
        default = [2.0, 11.0], help = 'Path to a file containing input parameters')

    return parser.parse_args()

args = parse_arguments()

# Read in parameters
param_file = file(args.parameter_file)
param = yaml.load(param_file)

# Read in initial PDB
res_freq, atom_coords = p2s.read_pdb_atom_data(args.input_pdb)

# Get volume from sluv for PDB
volume = sluv2.sum_volume(sluv2.all_residues, res_freq, 'perkins1986b')

# Optimize box side for sphere models
cutoff = param['sphere']['cutoff']
box_side, err = get_boxside.optimize_side(cutoff, atom_coords, volume,
                                          args.min_max[0],
                                          args.min_max[1],
                                          0.01)

dry_spheres, x_axis, y_axis, z_axis = p2s.create_sphere_model(atom_coords,
                                                              cutoff,
                                                              box_side)

# Calculate the expected volume of the hydration layer
# Hydration is estimated to be 0.3 g water per 1 g protein
# Bound water volume (< bulk volume) is given as a sluv parameter
water_vol = sluv2.calc_model_wet_volume(res_freq)

wet_volume = volume + water_vol

hydr_cutoff, err = hydrate.optimize_cut(box_side,
                                        dry_spheres,
                                        26,
                                        wet_volume,
                                        10,
                                        14,
                                        0.01)

out_file = open(args.output_file,'w')
out_file.write("Target dry volume: {0:7.4f}\n".format(volume))
out_file.write("Target hydrated volume: {0:7.4f}\n".format(wet_volume))
out_file.write("Sphere model box side: {0:7.4f}\n".format(box_side))
out_file.write("Hydration cutoff: {0:d}\n".format(hydr_cutoff))
