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
import yaml

import sct


def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(description= 'Optimize the box side used to create dry sphere models and cutoff used to creat a hydration layer that reproduce expected volumes from an input PDB.\n')

    parser.add_argument('-i','--input_pdb', nargs='?', type=str,
        help = 'Path to the input PDB file', required=True)

    parser.add_argument('-o','--output_file', nargs='?', type=str,
        default='optimized_params.txt', help = 'Path to output file')

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
res_freq, atom_coords = sct.pdb.read_pdb_atom_data(args.input_pdb)

# Get volume from sluv for PDB
volume = sct.seq.sum_volume(sct.seq.all_residues, res_freq, 'perkins1986b')

# Optimize box side for sphere models
cutoff = param['sphere']['cutoff']
box_side, err = sct.sphere.optimize_box_side(cutoff, atom_coords, volume,
                                          args.min_max[0],
                                          args.min_max[1],
                                          0.01)

dry_spheres, x_axis, y_axis, z_axis = sct.sphere.create_sphere_model(atom_coords,
                                                                     cutoff,
                                                                     box_side)

# Calculate the expected volume of the hydration layer
# Hydration is estimated to be 0.3 g water per 1 g protein
# Bound water volume (< bulk volume) is given as a sluv parameter
water_vol = sct.seq.calc_model_wet_volume(res_freq)

wet_volume = volume + water_vol

# Optimize the water cutoff used to filter out excessive/overlapping hydration
# spheres
# Use 26 additional spheres in teh trial hydration
# 10 - 14 is a reasonable window over which to refine the cutoff
hydr_cutoff, err = sct.sphere.optimize_cut(box_side,
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
