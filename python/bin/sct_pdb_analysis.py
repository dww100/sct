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
    @return:               Dictionary containing the following key/value pairs:
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

        sct.curve.smear_sas_curve(result['curve'],
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

def sas_model_summary_output(theor_data):
    """
    Format summary data from modelling for output. Columns are Rg from model, 
    Rg from curve, Rxs1 from curve, volume of sphere model.
    
    @type theor_data:   dictionary, dictionary
    @param theor_data:  Dictionaries containing the following key/value pairs 
                        for the dry model/neutron and wet model/xray 
                        comparisons:
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
                        - volume: Sphere model volume  
    @rtype:             string
    @return:            String containing the input for the sphere model Rg, 
                        curve Rg, curve Rxs1 and volume
    """
    
    if len(theor_data) > 0:
        summ = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\t".format(theor_data['model_rg'],
                                                                 theor_data['curve_rg'],
                                                                 theor_data['curve_rxs'],
                                                                 theor_data['volume'])        
        
        for dataset in theor_data['rfac']:
            summ += "{0:7.4f}\t{1:7.4f}\t".format(dataset[0],dataset[1])
    else:
        summ = "NA\tNA\tNA\tNA\tNA\tNA"
        
    return summ

def create_rfac_header_text(data_files):
    """
    Create the spacing for needed to provide two columns (Rfactor and scale) 
    for each input experimental curve in the summary output header. The header 
    has the following format (in the final version it is tab separated) except 
    without the line numbering:
    0. Path to input PDBs
    1. Neutron                                                                   X-ray                                       
    2. Model Rg_model Rg_curve Rxs1_curve Volume Rfactor scale) * neutron curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves
    This function provides the spacing for one set of curves (either x-ray or 
    neutron).
    
    @type data_files:   list
    @param data_files:  List of input files names containing experimental 
                        scattering curves
    @rtype:             list
    @return:            List of strings containing the spacing and column 
                        headings required for the three lines described above.
    """
    
    if data_files is not None:
    
        no_files = len(data_files)
        rfact_head = ['\t'.join(['\t']*no_files),
                      '\t\t'.join(data_files),
                      '\t'.join(['Rfactor\tscale']*no_files)]
    else:
        rfact_head = ['\t\t','\t\t','Rfactor\tscale']
                       
    return rfact_head

def  write_summary_header(in_pdb, in_neut, in_xray, output):
    """
    Print the header for the summary data from sphere model analysis. 
    When properties of the sphere model are references dry models are 
    associated with neutrons and hydrated (wet) ones with x-rays.
    The header has the following format (in the final version it is tab 
    separated) except without the line numbering:
    0. Path to input PDBs
    1. Neutron                                                                    X-ray                                       
    2. Model Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * neutron curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves
    
    Each experimental curve input requires two columns - Rfactor and scale
    Where the latter is the scaling factor used to match the theoretical and 
    experimental curves in the Rfactor calculation.
    
    @type theor_data:   dictionary, dictionary
    @param theor_data:  Dictionaries containing the following key/value pairs
    """
    # Print header - note the inputs and print as first header row
    output.write("Input PDB path: {0:s}\n".format(in_pdb))

    # Create appropriate spacing for headings of the Rfactor related columns 
    # for both neutron and xray input. 
    neut_rfact_head = create_rfac_header_text(in_neut)
    xray_rfact_head = create_rfac_header_text(in_xray)

    # Header 0
    summary_data.write("\tNeutron\t\t\t" + neut_rfact_head[0] + "\tX-ray\n")
    # Header 1
    summary_data.write("\t\t\t\t\t" + neut_rfact_head[1] + "\t\t\t\t" + xray_rfact_head[1] + "\n")
    # Header 2
    basic_head = "Rg_model\tRg_curve\tRxs1_curve\tVolume\t"
    col_head = "Model\t" + basic_head + neut_rfact_head[2] + '\t' + basic_head + xray_rfact_head[2] + "\n"
    output.write(col_head)

def perform_sas_analysis_pdb(pdb, neut_data, xray_data, param, radius, box_side3):
    """
    Create sphere models from PDBs and use to generate theoretical scattering
    curves and compare these to experimental inputs
    
    @type  pdb:        string
    @param pdb:        Filename of the input PDB
    @type  neut_data:  list
    @param neut_data:  List of numpy arrays with two columns (Q and I) 
                       containing neutron experimental data
    @type  xray_data:  list
    @param xray_data:  List of numpy arrays with two columns (Q and I)
                       containing xray experimental data
    @type  param:      dictionary
    @param param:      Dictionary containing parameters to use when creating
                       models and analysing curves.    
    @type  radius:     float
    @param radius:     Radius of spheres to be used in sphere models
    @type  box_side3:  float
    @param box_side3:  Cubed box length for the grid used to make sphere 
                       models (sphere diameter).
    @rtype:            dictionary, dictionary
    @return:           Dictionaries containing the following key/value pairs 
                       for the dry model/neutron and wet model/xray 
                       comparisons:
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
                       - volume: Sphere model volume
    """
    
    cutoff = param['sphere']['cutoff']
    box_side = param['sphere']['boxside']    
    
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
        if len(neut_data) != 0:

            neut_theor = analyse_sphere_model(dry_spheres, neut_data, radius, param, True)

            # Write curve to file
            curve_file = os.path.join(scn_path, pdb_id + '.scn')
            sct.curve.output_sas_curve(neut_theor['curve'], curve_file)

            # Write model to file
            model_file = os.path.join(dry_model_path,  pdb_id + '.pdb')
            sct.pdb.write_sphere_pdb(dry_spheres, radius, model_file)

            neut_theor['volume'] = box_side3 * len(dry_spheres)

        else:
            neut_theor = {}

        # If x-ray data provided compare with curve computed from wet sphere model
        if len(xray_data) != 0:
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
            sct.curve.output_sas_curve(xray_theor['curve'], curve_file)

            # Write model to file
            model_file = os.path.join(wet_model_path,  pdb_id + '.pdb')
            sct.pdb.write_sphere_pdb(wet_spheres, radius, model_file)

            xray_theor['volume'] = box_side3 * len(wet_spheres)

        else:
            xray_theor = {}
            
        return neut_theor, xray_theor

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

else:
    neut_data = []

if args.xray is not None:

    xray_data = read_curves(args.xray, args.xray_unit, param)

    scx_path = create_data_dir(args.output_path, 'xray','curves')
    wet_model_path = create_data_dir(args.output_path, 'xray','models')

else:
    xray_data = []

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

# Print the header for the summary data from sphere model analysis:
# Path to input PDBs
# Neutron                                                              X-ray                                       
# Model Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * neutron curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves
write_summary_header(args.input_path, args.neutron, args.xray, summary_data)

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

    # Create sphere models, compute scattering curves and compare to 
    # experimental curves
    # Dry models are compared to neutron data, wet to xray data.
    dry_data, wet_data = perform_sas_analysis_pdb(pdb, neut_data, xray_data, 
                                                  param, radius, box_side3)

    pdb_basename = os.path.basename(pdb)
    pdb_id = os.path.splitext(pdb_basename)[0]

    # Format the modelling output data for printing
    neut_summ = sas_model_summary_output(dry_data)
    xray_summ = sas_model_summary_output(wet_data)

    # Output all summary data to file
    summary_data.write('{0:s}\t{1:s}\t{2:s}\n'.format(pdb_id,
                                                      neut_summ,
                                                      xray_summ))

summary_data.close()
