#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute the R factor of a comparison between two scattering curves.
This will usually be one calculated and one experimental data set.

R factor is used by analogy with crystallography where:
R = sum (abs(F_expt - F_calc)) / sum (abs(F_expt))
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

    parser = argparse.ArgumentParser(
        description = 'Compute R factor from a comparison of two '
                      'scattering curves.')
    parser.add_argument(
        '-c','--calc', nargs='?', type=str, dest='calc_curve',
        help = 'Path to the input calculated curve', required=True)
    parser.add_argument(
        '-e','--expt', nargs='?', type=str, dest='expt_curve',
        help = 'Path to the input experimental curve', required=True)
    parser.add_argument(
        '-q','--qrange', nargs=2, type=float,
        help = 'Minimum and maximum Q values used in curve fitting',
        required=True)
    parser.add_argument(
        '-o','--outfile', nargs='?', type=str, dest='out_file',
        help = 'Path to the output file', required=False)
    parser.add_argument(
        '--header', action='store_true',
        help = 'Output a header (no data analysis performed)')

    args = parser.parse_args()
    return args

def main():
    # Interpret command line arguments
    args = parse_arguments()

    if (not args.header) :

        q_min = args.qrange[0]
        q_max = args.qrange[1]

        # Load the data from the two input files
        calc_data = sct.curve.load_scatter_curve(args.calc_curve, q_min, q_max)
        expt_data = sct.curve.load_scatter_curve(args.expt_curve, q_min, q_max)

        rfactor, scale = sct.curve.calculate_rfactor(expt_data,
                                                     calc_data,
                                                     q_min,
                                                     q_max)

        output_data = '{0:s}\t{1:s}\t{2:0.5f}\t{3:0.5f}\t{4:0.5f}\t{5:0.5f}\n'.format(
                                        args.calc_curve,
                                        args.expt_curve,
                                        q_min,
                                        q_max,
                                        1.0/scale,
                                        rfactor)

        if args.out_file == None:
            print output_data
        else:
            with open(args.out_file, "a") as fle:
                fle.write(output_data)
    else:

        header = 'Calculated file\tExperimental file\tQ min\tQ max\tScaling factor\tR factor\n'

        if args.out_file == None:
            print header
        else:
            with open(args.out_file, "a") as fle:
                fle.write(header)

if __name__ == "__main__":
    main()
