#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Automatic optimization of the box side used to create sphere models from
atomistic PDBs in SCT tools.
The target is the theoretical volume of the PDB structure calculated from sluv2
Takes a PDB, a range of box sizes to try, a cutoff for the number of atoms
per box and a percentage tolerance for the optimization as inputs.

Within the Perkins lab this replaces the do_cubeside script.
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
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description= 'Optimize the box size used in creation of sphere models\n'
        )

    parser.add_argument('-i','--input_pdb', nargs='?', type=str,
        dest='input_filename', help = 'Path to the input PDB file',
        required=True)

    parser.add_argument('-s','--input_seq', nargs='?', type=str,
        help = 'Path to input sequence file if different from PDB', 
        default = None)    

    parser.add_argument('-t','--input_type', choices = ['fas','yml'],
        help = 'Input file format (fasta or sluv yaml)', default = None)

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
        dest='output_filename', default=None,
        help = 'Path to the output file')

    parser.add_argument('-b','--box_range', nargs=2, type=float,
        default=[1.1,10.0], help = 'Minimum and maximum box sides to try')

    parser.add_argument('-t','--tolerance', nargs='?', type=float,
        default=0.01, help = 'Tolerance for box size determination')

    parser.add_argument('-c','--cutoff', nargs='?', type=int,
        default=4, help = 'Atom number cutoff for sphere models')

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

    volume = sct.seq.sum_volume(sct.seq.all_residues, res_freq, 'perkins1986b')    
    
    return volume, atom_coords
            

def main ():

    args = parse_arguments()

    target_volume, atom_coords = get_box_opt_input(args.input_pdb, 
                                                   args.input_seq, 
                                                   args.input_type)

    best_side, dev = sct.sphere.optimize_box_side(args.cutoff,
                                                  atom_coords, 
                                                  target_volume, 
                                                  args.box_range[0], 
                                                  args.box_range[1],
                                                  args.tolerance)

    if args.output_filename == None:
        output = sys.stdout
    else:
        output = open(args.output_filename,'w')

    output.write("# Optimized box_side for {0:s}, using cutoff {1:d}\n".format(
                 args.input_filename, args.cutoff))
    output.write("# Deviation from target volume was {0:f}\n".format(dev))
    output.write("box_side: {0:f}\n".format(best_side))

    output.close()

if __name__ == "__main__":
    main()
