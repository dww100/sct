#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate Rg from an input sphere model.
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
        description='Calculate radius of gyration from a SCT sphere model file\n')

    parser.add_argument(
        '-i',
        '--input_filename',
        nargs='?',
        type=str,
        dest='input_filename',
        help='Path to the input sphere model',
        required=True)

    parser.add_argument('-o', '--output_filename', nargs='?', type=str,
                        default=None, help='Path to the output file')

    return parser.parse_args()


def main():

    args = parse_arguments()

    # Read in SCT formatted sphere model file
    # -returns coordinates and sphere radius
    coords, radius = sct.sphere.read_mono_spheres(args.input_filename)

    # Calculate the radius of gyration for the sphere distribution
    r_gyr = sct.sphere.sphere_model_rg(coords, radius)

    if args.output_filename is not None:
        output = open(args.output_filename, 'w')
    else:
        output = sys.stdout

    output.write("Radius of gyration: {0:7.4f}\n".format(r_gyr))

if __name__ == "__main__":
    main()
