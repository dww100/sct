#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compares pre-calculated theorerical scattering curves to
experimental x-ray and neutron scattering curves.
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
import os
import glob
import yaml

import sct

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description = 'Compute R factor from a comparison of several scattering curves.')
    parser.add_argument(
        '-c','--calc', nargs='+', type=str, dest='calc_paths',
        help = 'Path to the input calculated curve or directory containing mutiple curves (with scn, scx extensions)', required=True)
    parser.add_argument(
        '-e','--expt', nargs='+', type=str, dest='expt_curves',
        help = 'List of input experimental curves', required=True)
    parser.add_argument('-p','--parameter_file', nargs='?', type=str,
        help = 'Path to a file containing input parameters', required=True)
    parser.add_argument(
        '-o','--outfile', nargs='?', type=str, dest='out_file',
        help = 'Path to the output file', required=True)
    parser.add_argument('-u','--unit', choices = ['nm', 'a'],
                        default = 'a', help = 'Unit for Q in input experimental data')

    args = parser.parse_args()
    return args

args = parse_arguments()

# Read in parameters
param_file = file(args.parameter_file)
param = yaml.load(param_file)

# Open file for output data
output = open(args.out_file,'w')

# Get list of scattering curves in the input directory

curve_filters = []
calc_curves = []

for path in args.calc_paths:
    if os.path.isdir(path):
        curve_filters.append(os.path.join(path,'*.scn'))
        curve_filters.append(os.path.join(path,'*.scx'))
    else:
        calc_curves.append(path)

for curve_filter in curve_filters:
    calc_curves.extend(glob.glob(curve_filter))

if len(calc_curves) < 1:
    print "No calculated curve files found to analyze"
    sys.exit(1)

if args.unit == 'nm':
    qmin = param['rfac']['qmin'] * 10
    qmax = param['rfac']['qmax'] * 10
else:
    qmin = param['rfac']['qmin']
    qmax = param['rfac']['qmax']

output.write("Expt\tCalc\tI0\tRfac\n")
for expt_curve in args.expt_curves:

    # Load experimental curve and if necessary convert Q to angstrom
    expt_data = sct.curve.load_scatter_curve(expt_curve, qmin, qmax)

    if args.unit == 'nm':
        expt_data[:,0] = expt_data[:,0] / 10.0

    for calc_curve in calc_curves:

        # Load the theoretically calculated scattering curve
        calc_data = sct.curve.load_scatter_curve(calc_curve,
                                                 param['rfac']['qmin'],
                                                 param['rfac']['qmax'])

        # Calculate the R factor comparing the calculated and experimental curves
        rfactor, scale = sct.curve.calculate_rfactor(expt_data,
                                                     calc_data,
                                                     param['rfac']['qmin'],
                                                     param['rfac']['qmax'])

        output.write("{0:s}\t{1:s}\t{2:7.4f}\t{3:7.4f}\n".format(expt_curve,
                                                                  calc_curve,
                                                                  1.0/scale,
                                                                  rfactor))
