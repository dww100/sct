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
import shutil

import sct


def parse_arguments():
    """
    Parse command line arguments and ensure correct combinations present
    """

    parser = argparse.ArgumentParser(
        description='Compare theoretical curves generated from PDB files to experimental SAS curves\n')

    parser.add_argument('-i', '--input_path', nargs='?', type=str,
                        help='Path to the input files (should be level above xray/neuton directories)', required=True)

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
        default=[])

    parser.add_argument(
        '-n',
        '--neutron',
        nargs='+',
        type=str,
        help='Paths to files containing experimental neutron scattering curve',
        default=[])

    parser.add_argument(
        '-t',
        '--title',
        nargs='?',
        type=str,
        help='Title to use for summary output file',
        default='sct_output')

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

    if (args.neutron and args.xray):
        model_type = 'both'
    elif args.neutron:
        model_type = 'neutron'
    else:
        model_type = 'xray'

    data = {}

    if model_type in ['neutron','both']:
        
        neutron_dir = os.path.join(args.input_path, 'neutron')
        
        neutron_curve_filter = os.path.join(neutron_dir, 'curves','*.scn')
        neut_curve_files = glob.glob(neutron_curve_filter)
        neut_curve_files = sort_nicely(neut_curve_files)        
        
        neutron_model_filter = os.path.join(neutron_dir, 'models','*.pdb')
        neut_model_files = glob.glob(neutron_model_filter)
        neut_model_files = sort_nicely(neut_model_files)
        
        if len(neut_curve_files) != len(neut_model_files):
            err = 'Neutron curve and model numbers do not match!'
            print(err)
            sys.ext(1)
        else:
            data['neut'] = zip(neut_curve_files, neut_model_files)
            
    if model_type in ['xray','both']:
        xray_dir = os.path.join(args.input_path, 'xray')
        
        xray_curve_filter = os.path.join(xray_dir, 'curves','*.scx')
        xray_curve_files = glob.glob(xray_curve_filter)
        xray_curve_files = sort_nicely(xray_curve_files)        
        
        xray_model_filter = os.path.join(xray_dir, 'models','*.pdb')
        xray_model_files = glob.glob(xray_model_filter)
        xray_model_files = sort_nicely(xray_model_files)
        
        if len(xray_curve_files) != len(xray_model_files):
            err = 'X-ray curve and model numbers do not match!'
            print(err)
            sys.ext(1)
        else:
            data['xray'] = zip(xray_curve_files, xray_model_files)    
    
    if model_type in 'neutron':
        nf = len(data['neut'])
    else:
        nf = len(data['xray'])

    # Read in parameters and check we have those we need for the workflow

    needed = ['curve', 'sphere', 'rg', 'rxs1', 'rfac']

    if model_type in ['xray', 'both']:
        needed.append('hydrate')

    param = sct.param.parse_parameter_file(args.parameter_file, needed)
    box_side3 = param['sphere']['boxside3']

    # Create output directory and open file for summary output
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    print "> Processing experimental data"

    title = args.title

    neut_expt = sct.curve.read_scatter_curves(
        args.neutron,
        args.neutron_unit,
        param)
    xray_expt = sct.curve.read_scatter_curves(args.xray, args.xray_unit, param)
    sct.tasks.output_expt_summary(
        neut_expt,
        xray_expt,
        args.output_path,
        title)

    print "> Processing calculated data"

    # Create the file for model output
    summary_name = os.path.join(args.output_path, title + '.sum')
    summary_data = open(summary_name, 'w')

    # Print the header for the summary data from sphere model analysis (Rxs2 added after Rxs1 if a range is supplied in param):
    # Path to input data
    #       Neutron                                                              X-ray
    # Model Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * neutron
    # curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves
    sct.tasks.write_summary_header(
        args.input_path,
        args.neutron,
        args.xray,
        param,
        summary_data,
        args.chi2)

    out_text = 'There are %d models to process\n' % (nf)
    print out_text

    for i in xrange(0,nf):
    
        if model_type in ['neutron','both']:
            neut_curve = data['neut'][i][0]
            neut_model = data['neut'][i][1]
        if model_type in ['xray','both']:
            xray_curve = data['xray'][i][0]
            xray_model = data['xray'][i][1]

        xray_theor = {}
        neut_theor = {}

        if model_type in ['neutron', 'both']:

            dry_spheres, radius = sct.pdb.read_sphere_pdb(neut_model)
            # Get Rg from sphere model
            neut_theor['model_rg'] = sct.sphere.sphere_model_rg(
                dry_spheres,
                radius)
            # Get volume of sphere model
            neut_theor['volume'] = box_side3 * len(dry_spheres)

            theor_curve = sct.curve.load_scatter_curve(
                neut_curve,
                param['rfac']['qmin'],
                param['rfac']['qmax'])

            # Calculate curve metrics - Rg, Rxs? and Rfactor/Chi^2 comparing
            # theoretical and experimental curves
            neut_theor.update(
                sct.tasks.analyse_theor_curve(
                    theor_curve,
                    neut_expt,
                    radius,
                    param,
                    chi2=args.chi2))

        if model_type in ['xray', 'both']:

            wet_spheres, radius = sct.pdb.read_sphere_pdb(xray_model)
            # Get Rg from sphere model
            xray_theor['model_rg'] = sct.sphere.sphere_model_rg(
                wet_spheres,
                radius)
            # Get volume of sphere model
            xray_theor['volume'] = box_side3 * len(wet_spheres)

            theor_curve = sct.curve.load_scatter_curve(
                xray_curve,
                param['rfac']['qmin'],
                param['rfac']['qmax'])

            # Calculate curve metrics - Rg, Rxs? and Rfactor/Chi^2 comparing
            # theoretical and experimental curves
            xray_theor.update(
                sct.tasks.analyse_theor_curve(
                    theor_curve,
                    xray_expt,
                    radius,
                    param,
                    chi2=args.chi2))

        file_basename = os.path.basename(data.itervalues().next()[i][0])
        file_id = os.path.splitext(file_basename)[0]

        # Format the modelling output data for printing
        neut_summ = sct.tasks.sas_model_summary_output(neut_theor, param)
        xray_summ = sct.tasks.sas_model_summary_output(xray_theor, param)

        # Output all summary data to file
        summary_data.write('{0:s}\t{1:s}{2:s}\n'.format(file_id,
                                                        neut_summ,
                                                        xray_summ))

        num_proc = i + 1
        fraction_done = (float(num_proc) / float(nf))
        progress_string = 'Processed ' + \
            str(num_proc) + ' of ' + str(nf) + ' : ' + str(fraction_done * 100.0) + ' % done\n'
        if (num_proc % 20) == 0:
            print(progress_string)


if __name__ == "__main__":
    main()
