#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create SASSIE style weights file from SCT output files
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

from __future__ import absolute_import
from __future__ import print_function
import sys, os
import argparse

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description='Create SASSIE style weights file from SCT output files\n')

    parser.add_argument('-i', '--infile', nargs='?', type=str,
                        help='Path to the input SCT data file', required=True)

    parser.add_argument(
        '-c',
        '--column',
        nargs='?',
        type=int,
        default=1,
        help='Number of column in the data file to be filtered on (columns numbered from 0)')

    parser.add_argument(
        '-m',
        '--min',
        nargs='?',
        type=float,
        default=0,
        help='Minimum acceptable value (default = 0)')

    parser.add_argument(
        '-x',
        '--max',
        nargs='?',
        type=float,
        default=100000,
        help='Maximum acceptable value (default = 100000)')

    parser.add_argument('-o', '--outfile', nargs='?', type=str,
                        help='Filename for the output weights file', required=True)

    return parser.parse_args()

args = parse_arguments()

infile = args.infile
filename, fileext = os.path.splitext(infile)
if fileext != '.sum':
    print("WARNING: Input does not use the standard SCT extension")

filter_col = args.column
outfile = args.outfile

filter_min = args.min
filter_max = args.max

if not filter_min and not filter_max:
    print("ERROR: No filter range selected\n")
    sys.exit()    

output = open(outfile, 'w')

title = "#Weight file for {0:s} column {1:d}, min-max filter: {2:f}-{3:f}\n".format(
            filename, filter_col, filter_min, filter_max)

output.write(title)

output.write("# Frame, Value, Weight\n")

header_count = 0
stucture_count = 0

try:

    for line in open(infile, 'r'):
        header_count += 1
        if header_count == 4:
            col_head = line.split()[filter_col]
            if col_head not in ['Rg_model','Rg_curve','Rfactor']:
                print("ERROR: Selected column unsuitable for filtering: " + col_head) 
                print("ERROR: Column heading must be: Rg_model,Rg_curve or Rfactor")
                print("ERROR: Column heading must be: Rg_model,Rg_curve or Rfactor\n")
                sys.exit()
        if header_count > 4:
            
            cols = line.split()
            
            stucture_count += 1
            
            col_value = float(cols[filter_col])     
            
            if  filter_min <= col_value <= filter_max:
                        
                output.write(
                    "{0:d}\t{1:f}\t{2:f}\n".format(stucture_count, col_value, 1.0))
            else:
                output.write(
                    "{0:d}\t{1:f}\t{2:f}\n".format(stucture_count, col_value, 0.0))

    output.close()
    
except:
    print("ERROR: Unable to process file\n")
    
