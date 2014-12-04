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
import re
import glob
import yaml

import sct


def parse_arguments():
    """
    Parse command line arguments and ensure correct combinations present
    """

    parser = argparse.ArgumentParser(
        description='Compare theoretical curves generated from PDB files to experimental SAS curves\n')

    parser.add_argument('-i', '--input_path', nargs='?', type=str,
                        help='Path to the input PDB files', required=True)

    parser.add_argument('-o', '--output_path', nargs='?', type=str,
                        default='.', help='Path in which to save output files')

    parser.add_argument(
        '-p',
        '--parameter_file',
        nargs='?',
        type=str,
        help='Path to a file containing input parameters',
        required=True)

    parser.add_argument(
        '-x',
        '--xray',
        nargs='+',
        type=str,
        help='Paths to files containing experimental x-ray scattering curve',
        default=None)

    parser.add_argument(
        '-n',
        '--neutron',
        nargs='+',
        type=str,
        help='Paths to files containing experimental neutron scattering curve',
        default=None)

    parser.add_argument(
        '-t',
        '--title',
        nargs='?',
        type=str,
        help='Title to use for summary output file',
        default='sct_output')

    parser.add_argument(
        '-a',
        '--add_res',
        nargs='?',
        type=str,
        default=None,
        help='Path to YAML file containing mass and volume for residues not originally used by sluv/SCT')

    parser.add_argument('-xu', '--xray_unit', choices=['nm', 'a'],
                        default='a', help='Unit for Q in input x-ray data')

    parser.add_argument('-nu', '--neutron_unit', choices=['nm', 'a'],
                        default='a', help='Unit for Q in input neutron data')

    parser.add_argument('-ou', '--output_unit', choices=['nm', 'a'],
                        default='a', help='Unit for Q in output data')

    parser.add_argument(
        '--chi2',
        action='store_true',
        default=False,
        help='Select comparison metric to be Chi squared not R factor')

    args = parser.parse_args()

    if (args.neutron is None) and (args.xray is None):
        print "At least one experimental curve is required for comparison (xray, neutron or both).\n"
        sys.exit(1)

    return args

# Human sort taken from from:
# http://nedbatchelder.com/blog/200712/human_sorting.html
# Provided with no formal license but site indicates that you can do with
# it as you wish


def tryint(s):

    if s.isdigit():
        return int(s)
    else:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryint(c) for c in re.split('([0-9]+)', s.lower())]


def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """

    return sorted(l, key=alphanum_key)


def main():

    print "Running modern SCT workflow"
    print "---------------------------\n"

    args = parse_arguments()

    if args.add_res:
        res_file = file(args.add_res, 'r')
        add_res = yaml.load(res_file)
        for res, data in add_res.iteritems():
            sct.seq.all_residues.append(res)
            sct.seq.res_vols['perkins1986a']['residue'][res] = data['vol']
            sct.seq.params['mass'][res] = data['mass']

    # Read in parameters and check we have those we need for the workflow

    needed = ['curve', 'sphere', 'rg', 'rxs1', 'rfac']

    if args.xray is not None:
        needed.append('hydrate')

    param = sct.param.parse_parameter_file(args.parameter_file, needed)

    # Create output directory and open file for summary output
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    print "> Processing experimental data"
    # Read in experimental curves and calculate Rg and Rxs
    # Setup output directories for theoretical curves and sphere models
    # Output summary analysis of the experimental data curves
    neut_data, xray_data, out_paths = sct.tasks.process_expt_data(
        args.neutron, args.neutron_unit, args.xray, args.xray_unit, args.output_path, args.title, param)

    # Create the file for model output
    summary_name = os.path.join(args.output_path, args.title + '.sum')
    summary_data = open(summary_name, 'w')

    # Print the header for the summary data from sphere model analysis (Rxs2 added after Rxs1 if a range is supplied in param):
    # Path to input PDBs
    # Neutron                                                              X-ray
    # Model Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * neutron
    # curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves
    sct.tasks.write_summary_header(
        args.input_path,
        args.neutron,
        args.xray,
        param,
        summary_data,
        args.chi2)

    # Get list of PDBs in the input directory
    pdb_filter = os.path.join(args.input_path, '*.pdb')
    pdb_files = glob.glob(pdb_filter)

    if len(pdb_files) < 1:
        print "Error: No PDB files found to analyze"
        sys.exit(1)

    # Sort files so that they are in human expected alpha numerical order
    # This means that XXXX2.pdb will sort before XXXX100.pdb
    pdb_files = sort_nicely(pdb_files)

    print "> Analyzing input PDBs"

    # Loop over input PDBs
    for count, pdb in enumerate(pdb_files):

        if (count % 20 == 0):
            print "\tProcessing PDB number " + str(count)

        try:
            # Create sphere models, compute scattering curves and compare to
            # experimental curves
            # Dry models are compared to neutron data, wet to xray data.
            dry_data, wet_data = sct.tasks.perform_sas_analysis_pdb(pdb,
                                                                    neut_data,
                                                                    xray_data,
                                                                    param,
                                                                    out_paths)

        except IOError as e:
            print "Error loading PDB file name %s: %s" % (args.input_filename, e)
            continue

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

    print "Done"

if __name__ == "__main__":
    main()
