#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate a theoretical scattering curve from a sphere model input.

Note: This code ignores command line optons silently if a parameter
file is passed.
"""

# Copyright 2014 University College London

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse

import sct


def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description='Calculate SAS curve from sphere model\n')

    parser.add_argument(
        '-i',
        '--input_filename',
        nargs='?',
        type=str,
        dest='input_filename',
        help='Path to the input sphere or PDB file',
        required=True)

    parser.add_argument(
        '-f',
        '--input_format',
        choices=[
            'pdb',
            'sph'],
        help='Input file format (PDB or SCT sphere)',
        default='pdb')

    parser.add_argument('-o', '--output_filename', nargs='?', type=str,
                        dest='output_filename', default=None,
                        help='Path to the output file')

    parser.add_argument(
        '-p',
        '--parameter_file',
        nargs='?',
        type=str,
        help='Path to a file containing input parameters',
        default=None)

    parser.add_argument('-q', '--q_max', nargs='?', type=float,
                        default=0.16, help='Maximum q value in output curve')

    parser.add_argument('-x', '--xray', action='store_true',
                        help='Flag to specify this is an x-ray curve')

    parser.add_argument('-r', '--radius', nargs='?', type=float,
                        default=3.77, help='Sphere radius')

    parser.add_argument(
        '-b',
        '--n_bins',
        nargs='?',
        type=int,
        default=400,
        help='No. bins to use in histogram of sphere separation')

    parser.add_argument('-n', '--n_points', nargs='?', type=int,
                        default=100, help='No. points in output curve')

    parser.add_argument('-s', '--smear', action='store_true',
                        help='Apply smearing to neutron curve')

    parser.add_argument(
        '-e',
        '--spread',
        nargs='?',
        type=float,
        default=0.1,
        help='Wavelength spread used to calculate smearing')

    parser.add_argument(
        '-w',
        '--wavelength',
        nargs='?',
        type=float,
        default=6.0,
        help='Wavelength used to calculate smearing')

    parser.add_argument(
        '-d',
        '--divergence',
        nargs='?',
        type=float,
        default=0.016,
        help='Beam divergence used to calculate smearing')

    return parser.parse_args()


def main():

    args = parse_arguments()

    if args.input_format == 'pdb':
        res_freq, coords = sct.pdb.read_pdb_atom_data(args.input_filename)
    else:
        coords, radius = sct.sphere.read_mono_spheres(args.input_filename)

    if args.parameter_file is None:
        curve = sct.sphere.spheres_to_sas_curve(coords, args.radius,
                                                args.q_max, args.n_points,
                                                rbins=args.n_bins)

        if args.smear and not args.xray:
            q_delta = args.q_max / args.n_points
            sct.curve.smear_sas_curve(curve, q_delta, args.wavelength,
                                      args.spread, args.divergence)

    else:

        print "WARNING: A SCT parameter file was specified, so the modelling parameters from the command line flags will be ignored!"

        # Read in parameters
        param = sct.param.parse_parameter_file(
            args.parameter_file, [
                'curve', 'sphere'])

        curve = sct.sphere.spheres_to_sas_curve(
            coords,
            param['sphere']['boxside'],
            param['curve']['qmax'],
            param['curve']['npoints'],
            rbins=param['curve']['radbins'])

        if param['curve']['smear'] and not args.xray:
            sct.curve.smear_sas_curve(
                curve,
                param['curve']['q_delta'],
                param['curve']['wavelength'],
                param['curve']['spread'],
                param['curve']['divergence'])

    sct.curve.output_sas_curve(curve, args.output_filename)

if __name__ == "__main__":
    main()
