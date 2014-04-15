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
import yaml

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

    args = parser.parse_args()

    if (args.neutron == None) and (args.xray == None):
        print "At least one experimental curve is required for comparison (xray, neutron or both).\n"
        sys.exit(1)

    return args

def read_curves(curve_files, units, param):
    """
    Read in a list of scattering curve files and return a list of dictionaries
    which contain the file name ('file'), an array of Q and I(Q) vales ('data'),
    the radius of gyration ('rg') and the cross-section ('rxs1').

    @type  curve_files:  list
    @param curve_files:  List of files containing scattering curves Q, I(Q)
                         pairs
    @type  units:        string
    @param units:        String containing either 'nm' or 'a' to indicate the
                         length units used in the scattering curve
    @type  param:        dictionary
    @param param:        Dictionary containing parameters to use when creating
                         models and analysing curves.
    @return:             List of dictionaries containing the following key/value
                         pairs:
                         - data: numpy array of Q, I(Q)
                         - rg: radius of gyration calculated from input curve
                         - rxs1: cross-section calculated from input curve
                         - file: input file name
    """

    curves = []

    for curve_file in curve_files:

        curve = {}

        curve['file'] = curve_file
        # Read in the scattering curve
        # Modeling is performed in angstroms so convert files in nm to a
        if units == 'nm':
            curve['data'] = sct.curve.load_scatter_curve(curve_file,
                                                param['rfac']['qmin'] * 10.0,
                                                param['rfac']['qmax'] * 10.0)
            curve['data'][:,0] = curve['data'][:,0] / 10.0
        else:
            curve['data'] = sct.curve.load_scatter_curve(curve_file,
                                                param['rfac']['qmin'],
                                                param['rfac']['qmax'])


        curve['rg'], curve['rxs'] = sct.curve.get_curve_descriptors(curve['data'],
                                                  param['rg']['fitmin'],
                                                  param['rg']['fitmax'],
                                                  param['rxs1']['fitmin'],
                                                  param['rxs1']['fitmax'])

        curves.append(curve)

    return curves

def create_data_dir(basename, expt_type, data_type):
    """
    Set the output to a path basename/expt_type/data_type, if this does not
    exist create it.

    @type  basename:   string
    @param basename:   Basename for the output directory
    @type  expt_type:  string
    @param expt_type:  Experiment type
    @type  data_type:  string
    @param data_type:  Type of data being stored

    """

    full_path = os.path.join(basename, expt_type, data_type)
    if not os.path.exists(full_path):
        os.makedirs(full_path)

    return full_path

def analyse_sphere_model(model, expt_curves, sphere_radius, param, neutron=False):
    """
    Take a sphere model and calculate the theoretical scattering curve. The
    Rg and Rxs1 are also calculated from the curve and returned alongside the
    rfac comparing the theoretical curve to an experimental one.

    @type  model:          list
    @param model:          List of lists containing x, y & z coordinates
                           (3 * floats) for each sphere in a sphere model.
    @type  expt_curves:    list
    @param expt_curves:    List of two dimensional numpy arrays containing
                           scattered vector magnitudes, q, and intensities, I,
                           from experiment.
    @type  sphere_radius:  float
    @param sphere_radius:  Sphere radius
    @type  param:          dictionary
    @param param:          Dictionary containing parameters to use when creating
                           models and analysing curves.
    @param neutron:        Flag set if aiming to create a theoretical
                           neutron curve
    @type  neutron:        boolean
    @rtype:                dictionary
    @return:                Dictionary containing the following key/value pairs:
                           - model_rg: Radius of gyration calculated directly
                           from sphere model.
                           - curve_rg: Radius of gyration calculated from the
                           theoretical scattering curve derived from the sphere
                           model.
                           - curve_rxs1: Cross-section calculated from the
                           theoretical scattering curve derived from the sphere
                           model.
                           - rfac: R factor comparing experimental and
                           theoretical scattering curves.
    """

    result = {}

    # Calculate Rg from sphere model
    result['model_rg'] = sct.sphere.sphere_model_rg(model, sphere_radius)
    # Create theoretical curve from model
    result['curve'] = sct.sphere.spheres_to_sas_curve(model,
                                                sphere_radius,
                                                param['curve']['qmax'],
                                                param['curve']['npoints'],
                                                rbins = param['curve']['radbins'])

    # Rg and Rxs from theoretical curve
    result['curve_rg'], result['curve_rxs'] = sct.curve.get_curve_descriptors(result['curve'],
                                                            param['rg']['fitmin'],
                                                            param['rg']['fitmax'],
                                                            param['rxs1']['fitmin'],
                                                            param['rxs1']['fitmax'])

    # Neutron curves are usually smeared with parameters for the instrument used
    if (neutron and param['curve']['smear']):

        sct.curve.smear_sas_curve(result['i'], result['q'],
                                    param['curve']['q_delta'],
                                    param['curve']['wavelength'],
                                    param['curve']['spread'],
                                    param['curve']['divergence'])

    # Calculate Rfactors for theoretical vs experimental curves
    result['rfac'] = []
    for expt_curve in expt_curves:
        result['rfac'].append(sct.curve.calculate_rfactor(expt_curve['data'],
                                                          result['curve'],
                                                          param['rfac']['qmin'],
                                                          param['rfac']['qmax']))

    return result

args = parse_arguments()

# Read in parameters
param_file = file(args.parameter_file)
param = yaml.load(param_file)
param['curve']['q_delta'] = param['curve']['qmax'] / param['curve']['npoints']

# Create output directory and open file for summary output
if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

# Read in experimental curves and calculate Rg and Rxs
# Setup output directories for theoretical curves and sphere models
if args.neutron is not None:

    neut_data = read_curves(args.neutron, args.neutron_unit, param)

    scn_path = create_data_dir(args.output_path, 'neutron','curves')
    dry_model_path = create_data_dir(args.output_path, 'neutron','models')

    no_neut = len(neut_data)
    neut_rfact_head = ['\t'.join(['\t']*no_neut),
                       '\t\t'.join(args.neutron),
                       '\t'.join(['Rfactor\tscale']*no_neut)]

else:
    neut_data = []
    neut_rfact_head = ['\t\t','\t\t','Rfactor\tscale']

if args.xray is not None:

    xray_data = read_curves(args.xray, args.xray_unit, param)

    scx_path = create_data_dir(args.output_path, 'xray','curves')
    wet_model_path = create_data_dir(args.output_path, 'xray','models')

    no_xray = len(xray_data)
    xray_rfact_head = ['\t'.join(['\t']*no_xray),
                       '\t\t'.join(args.xray),
                       '\t'.join(['Rfactor\tscale']*no_xray)]

else:
    xray_data = []
    xray_rfact_head = ['\t\t','\t\t','Rfactor\tscale']

# Output summary analysis of the experimental data curves
expt_name = os.path.join(args.output_path, args.title + '_expt.sum')
expt_data = open(expt_name,'w')
expt_data.write("Filename\tRg\tRxs1\n")

for curve in neut_data + xray_data:
    expt_data.write("{0:s}\t{1:7.4f}\t{2:7.4f}\n".format(curve['file'],
                                                         curve['rg'],
                                                         curve['rxs']))
expt_data.close()

# Create the file for model output
summary_name = os.path.join(args.output_path, args.title + '.sum')
summary_data = open(summary_name,'w')

# Print header - note the inputs and put in column headings
summary_data.write("Input PDB path: {0:s}\n".format(args.input_path))
summary_data.write("\tNeutron\t\t\t" + neut_rfact_head[0] + "\tX-ray\n")
summary_data.write("\t\t\t\t\t" + neut_rfact_head[1] + "\t\t\t\t" + xray_rfact_head[1] + "\n")

basic_head = "Rg_model\tRg_curve\tRxs1_curve\tVolume\t"
col_head = "Model\t" + basic_head + neut_rfact_head[2] + '\t' + basic_head + xray_rfact_head[2] + "\n"
summary_data.write(col_head)

# Get list of PDBs in the input directory
pdb_filter = os.path.join(args.input_path, '*.pdb')
pdb_files = glob.glob(pdb_filter)

if len(pdb_files) < 1:
    print "No PDB files found to analyze"
    sys.exit(1)

# Parameters for turning PDB into spheres
cutoff = param['sphere']['cutoff']
box_side = param['sphere']['boxside']
box_side3 = box_side**3

# Model spheres just fit into grid box (no overlap) therefore:
radius = box_side / 2.0

# Loop over input PDBs
for pdb in pdb_files:

    pdb_basename = os.path.basename(pdb)
    pdb_id = os.path.splitext(pdb_basename)[0]

    # Read PDB
    res_freq, atom_coords = sct.pdb.read_pdb_atom_data(pdb)

    if len(atom_coords) > 0:

        # PDB to dry sphere model
        dry_spheres, x_axis, y_axis, z_axis = sct.sphere.create_sphere_model(atom_coords,
                                                                             cutoff,
                                                                             box_side)

        # If neutron data provided compare with curve computed from dry sphere model
        if args.neutron is not None:

            neut_theor = analyse_sphere_model(dry_spheres, neut_data, radius, param, True)

            # Write curve to file
            curve_file = os.path.join(scn_path, pdb_id + '.scn')
            output = open(curve_file,'w')
            sct.curve.output_sas_curve(neut_theor['curve'], output)
            output.close()

            # Write model to file
            model_file = os.path.join(dry_model_path,  pdb_id + '.pdb')
            sct.pdb.write_sphere_pdb(dry_spheres, radius, model_file)

            volume = box_side3 * len(dry_spheres)

            # Format results for output to file
            neut_summ = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\t".format(neut_theor['model_rg'],
                                                                     neut_theor['curve_rg'],
                                                                     neut_theor['curve_rxs'],
                                                                     volume)
            for dataset in neut_theor['rfac']:
                neut_summ += "{0:7.4f}\t{1:7.4f}\t".format(dataset[0],dataset[1])

        else:
            neut_summ = "NA\tNA\tNA\tNA\tNA\tNA"

        # If x-ray data provided compare with curve computed from wet sphere model
        if args.xray is not None:
            # Hydrate model - > wet sphere model
            wet_spheres = sct.sphere.hydrate_sphere_model(dry_spheres,
                                                          param['hydrate']['positions'],
                                                          box_side,
                                                          param['hydrate']['cutoff'],
                                                          xaxis = x_axis,
                                                          yaxis = y_axis,
                                                          zaxis = z_axis)

            xray_theor = analyse_sphere_model(wet_spheres, xray_data, radius, param)

            # Write curve to file
            curve_file = os.path.join(scx_path, pdb_id + '.scx')
            output = open(curve_file,'w')
            sct.curve.output_sas_curve(xray_theor['curve'], output)
            output.close()

            # Write model to file
            model_file = os.path.join(wet_model_path,  pdb_id + '.pdb')
            sct.pdb.write_sphere_pdb(wet_spheres, radius, model_file)

            volume = box_side3 * len(wet_spheres)

            # Format results for output to file
            xray_summ = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\t".format(xray_theor['model_rg'],
                                                                           xray_theor['curve_rg'],
                                                                           xray_theor['curve_rxs'],
                                                                           volume)

            for dataset in xray_theor['rfac']:
                xray_summ += "{0:7.4f}\t{1:7.4f}\t".format(dataset[0],dataset[1])

        else:
            xray_summ = "NA\tNA\tNA\tNA\tNA\tNA"

        # Output all summary data to file
        summary_data.write('{0:s}\t{1:s}\t{2:s}\n'.format(pdb_id,
                                                          neut_summ,
                                                          xray_summ))

summary_data.close()
