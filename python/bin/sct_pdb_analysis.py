#!/usr0/usr3/people/davidw/linux_apps/anaconda/bin/python
# -*- coding: utf-8 -*-
"""Compares theorerical scattering curves generated from PDB structure files to
experimental x-ray and neutron scattering curves.

PDBs are converted into sphere models and the Debye equation used to compute
the theoretical curves. The original sphere models are surrounded with a
hydration layer of spheres before the creation of curves to be compared to x-ray
data.

Within the Perkins lab this replaces the do_curve script
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

import sct

def parse_arguments():
    """
    Parse command line arguments and ensure correct combinations present
    """

    parser = argparse.ArgumentParser(
        description= 'Compare theoretical curves generated from PDB files to experimental SAS curves\n')

    parser.add_argument('-i','--input_path', nargs='?', type=str,
        help = 'Path to the input PDB files', required=True)

    parser.add_argument('-o','--output_path', nargs='?', type=str,
        default='.', help = 'Path in which to save output files')

    parser.add_argument('-p','--parameter_file', nargs='?', type=str,
        help = 'Path to a file containing input parameters', required=True)

    parser.add_argument('-x','--xray', nargs='+', type=str,
        help = 'Paths to files containing experimental x-ray scattering curve', default = None)

    parser.add_argument('-n','--neutron', nargs='+', type=str,
        help = 'Paths to files containing experimental neutron scattering curve', default = None)

    parser.add_argument('-t','--title', nargs='?', type=str,
        help = 'Title to use for summary output file', default = 'sct_output')

    parser.add_argument('-xu','--xray_unit', choices = ['nm', 'a'],
                        default = 'a', help = 'Unit for Q in input x-ray data')

    parser.add_argument('-nu','--neutron_unit', choices = ['nm', 'a'],
                        default = 'a', help = 'Unit for Q in input neutron data')

    parser.add_argument('-ou','--output_unit', choices = ['nm', 'a'],
                        default = 'a', help = 'Unit for Q in output data')
                        
    parser.add_argument('--chi2', action='store_true', default=False)
    
    args = parser.parse_args()

    if (args.neutron == None) and (args.xray == None):
        print "At least one experimental curve is required for comparison (xray, neutron or both).\n"
        sys.exit(1)

    return args

def main():

    args = parse_arguments()
    
    # Read in parameters and check we have those we need for the workflow

    needed = ['curve','sphere','rg','rxs1','rfac']
    
    if args.xray is not None:
        needed.append('hydrate')
    
    param = sct.param.parse_parameter_file(args.parameter_file, needed)
    
    # Create output directory and open file for summary output
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    # Read in experimental curves and calculate Rg and Rxs
    # Setup output directories for theoretical curves and sphere models
    # Output summary analysis of the experimental data curves
    neut_data, xray_data, out_paths = sct.tasks.process_expt_data(args.neutron, 
                                                        args.neutron_unit, 
                                                        args.xray, 
                                                        args.xray_unit, 
                                                        args.output_path, 
                                                        args.title, 
                                                        param)
    
    # Create the file for model output
    summary_name = os.path.join(args.output_path, args.title + '.sum')
    summary_data = open(summary_name,'w')
    
    # Print the header for the summary data from sphere model analysis (Rxs2 added after Rxs1 if a range is supplied in param):
    # Path to input PDBs
    # Neutron                                                              X-ray                                       
    # Model Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * neutron curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves
    sct.tasks.write_summary_header(args.input_path, args.neutron, args.xray, param, summary_data, args.chi2)
    
    # Get list of PDBs in the input directory
    pdb_filter = os.path.join(args.input_path, '*.pdb')
    pdb_files = glob.glob(pdb_filter)
    
    if len(pdb_files) < 1:
        print "Error: No PDB files found to analyze"
        sys.exit(1)
    
    # Loop over input PDBs
    for pdb in pdb_files:
    
        # Create sphere models, compute scattering curves and compare to 
        # experimental curves
        # Dry models are compared to neutron data, wet to xray data.
        dry_data, wet_data = sct.tasks.perform_sas_analysis_pdb(pdb, 
                                                      neut_data, 
                                                      xray_data, 
                                                      param, 
                                                      out_paths)
    
        pdb_basename = os.path.basename(pdb)
        pdb_id = os.path.splitext(pdb_basename)[0]
    
        # Format the modelling output data for printing
        neut_summ = sct.tasks.sas_model_summary_output(dry_data, param)
        xray_summ = sct.tasks.sas_model_summary_output(wet_data, param)
    
        # Output all summary data to file
        summary_data.write('{0:s}\t{1:s}{2:s}\n'.format(pdb_id,
                                                        neut_summ,
                                                        xray_summ))
    
    summary_data.close()

if __name__ == "__main__":
    main()