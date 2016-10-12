#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert file containing sphere co-ordinates and radii into a PDB.
Each sphere is set to be a C1 atom.
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

from __future__ import absolute_import
import argparse

import sct


def parse_arguments():
    """
    Parse command line arguments and ensure correct combinations present
    """

    parser = argparse.ArgumentParser(
        description='Convert atomic pdb to sphere model\n')

    parser.add_argument('-i', '--input_filename', nargs='?', type=str,
                        dest='input_filename', help='Path to the input file',
                        required=True)

    parser.add_argument('-o', '--output_filename', nargs='?', type=str,
                        dest='output_filename', default=None,
                        help='Path to the output file', required=True)

    return parser.parse_args()


def main():

    args = parse_arguments()

    out_file = open(args.output_filename, 'w')

    # We will count the number of spheres read
    count = 1

    # Read input from spheres file and output each sphere with a different
    # atom and residue number
    with open(args.input_filename) as f:
        for line in f:
            sphere = sct.sphere.parse_sphere_line(line)
            # As in SJP original codes use SER residue type and C1 for atom type
            # Store the sphere radius at the temperature factor (beta)
            pdb_line = sct.pdb.create_pdb_atom(
                count,
                'SER',
                count,
                'C1',
                sphere['coords'],
                beta=sphere['radius'])
            out_file.write(pdb_line)
            count += 1

    out_file.close()

if __name__ == "__main__":
    main()
