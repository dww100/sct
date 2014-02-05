#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Take sequence from a PDB or FASTA file and return a yaml file
containing the frequencies of all residue types.
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
        description= 'Convert fasta or pdb file to yaml residue frequency file\n')
    parser.add_argument('-i','--input_filename', nargs='?', type=str,
        dest='input_filename', help = 'Path to the input fasta/pdb file',
        required=True)
    parser.add_argument('-p', '--pdb', action='store_true', default=False,
        help = 'Flag to indicate input is a PDB rather than fasta file')
    parser.add_argument('-o','--output_filename', nargs='?', type=str,
        dest='output_filename', help = 'Path to the output file')

    return parser.parse_args()

args = parse_arguments()

# Read input
# Input should either be PDB or fasta format - chosen by input flag
if args.pdb:
    res_freq, coords = sct.pdb.pdb_res_freq(args.input_filename)
else:
    res_freq = sct.seq.fasta_res_freq(args.input_filename)

if args.output_file != None:
    sys.stdout = open(args.output_file,'w')

# Output pdb frequencies in a yaml style format that can be read by SLUV2
for res_name in sct.seq.amino_acids:
    print res_name + ': ' + str(res_freq[res_name]) + '\n'
for res_name in sct.seq.monosaccharides:
    print res_name + ': ' + str(res_freq[res_name]) + '\n'
