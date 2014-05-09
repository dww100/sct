#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Add hydration shell of spheres around an existing sphere model read from a SCT
format sphere file.
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
    """
    Parse command line arguments and ensure correct combinations present
    """

    parser = argparse.ArgumentParser(
        description= 'Add hydration shell of spheres around an existing sphere model\n')

    parser.add_argument('-i','--input_filename', nargs='?', type=str,
                        dest='input_filename', required=True,
                        help = 'Path to the input file')

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
                        dest='output_filename', default=None,
                        help = 'Path to the output file', required=True)

    parser.add_argument('-p','--parameter_file', nargs='?', type=str,
        help = 'Path to a file containing input parameters', default=None)

    parser.add_argument('-n','--hydration_no', nargs='?', type=int,
                        default=26,
                        help = 'No. spheres to add as a hydration shell (1-26)')

    parser.add_argument('-c','--cutoff', nargs='?', type=int,
                        default=5.0,
                        help = 'Cutoff used to reduce the number of water spheres')

    return parser.parse_args()

def main ():

    args = parse_arguments()

    # Select number of solvent positions to occupy before filtering
    # Position 0 = original sphere position,
    # Positions 1 to 26 positions on cube centred on original sphere

    if args.parameter_file == None:    
        hydration_no = args.hydration_no + 1
        cutoff = args.cutoff
    else:
        # Read in parameters
        print "WARNING: A SCT parameter file was specified, so the modelling parameters from the command line flags will be ignored!"        
        
        # Read in parameters
        param, err = sct.param.read_parameter_file(args.parameter_file)
        
        if err != None:
            sct.param.output_error(err, args.parameter_file)
            
        err = sct.param.check_parameters(param, ['hydrate'])
        
        if err != None:
            sct.param.output_error(err, args.parameter_file)
            
        hydration_no = param['hydrate']['positions']
        cutoff = param['hydrate']['cutoff']

    # Read dry sphere model from file
    dry_spheres, radius = sct.sphere.read_mono_spheres(args.input_filename)

    # Create hydrated model
    wet_spheres = sct.sphere.hydrate_sphere_model(dry_spheres, hydration_no, radius, cutoff)

    out = open(args.output_filename, 'w')
    sct.sphere.write_spheres(wet_spheres, radius, out)
    out.close()

if __name__ == "__main__":
    main()
