#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Compares theorerical scattering curves generated from PDB structure files to
experimental x-ray and neutron scattering curves.

PDBs are converted into sphere models and the Debye equation used to compute
the theoretical curves. The original sphere models are surrounded with a
hydration layer of spheres before the creation of curves to be compared to x-ray
data.

Within the Perkins lab this replaces the do_curve script
"""

import argparse
import sys
import os
import glob
import pdb2sphere as p2s
import sphere2pdb as s2pdb
import hydrate_spheres as hydrate
import sas_curve_analysis
import calculate_curve
import calculate_rfactor as rfactor
import sphere_rg
import yaml
import numpy as np

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description= 'Compare theoretical curves generated from PDB files to experimental SAS curves\n')

    parser.add_argument('-i','--input_path', nargs='?', type=str,
        help = 'Path to the input PDB files', required=True)

    parser.add_argument('-o','--output_path', nargs='?', type=str,
        default='.', help = 'Path in which to save output files')

    parser.add_argument('-p','--parameter_file', nargs='?', type=str,
        help = 'Path to a file containing input parameters', required=True)

    parser.add_argument('-x','--xray', nargs='?', type=str,
        help = 'Path to a file containing experimental x-ray scattering curve', default = None)

    parser.add_argument('-n','--neutron', nargs='?', type=str,
        help = 'Path to a file containing experimental neutron scattering curve', default = None)

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

def format_qi_data(q,i):
    """Convert separate lists/arrays of q and i values into a single 2d array"""
    data = []

    for element in range(0, len(q)):
        data.append([q[element], i[element]])

    return np.array(data)

def analyse_model(model, expt_curve, sphere_radius, param):
    """Take a sphere model and calculate the theoretical scattering curve. The
    Rg and Rxs1 are also calculated from the curve and returned alongside the
    rfac comparing the theoretical curve to an experimental one."""
    result = {}

    # Calculate Rg from sphere model
    result['model_rg'] = sphere_rg.sphere_model_rg(model, sphere_radius)
    # Create theoretical curve from model
    result['q'], result['i'] = calculate_curve.spheres_to_sas_curve(model,
                                                sphere_radius,
                                                param['curve']['qmax'],
                                                param['curve']['npoints'],
                                                rbins = param['curve']['radbins'])

    theor_curve = format_qi_data(result['q'], result['i'])
    # Rg and Rxs from theoretical curve
    result['curve_rg'], result['curve_rxs'] = get_curve_descriptors(theor_curve,
                                                            param['rg']['fitmin'],
                                                            param['rg']['fitmax'],
                                                            param['rxs1']['fitmin'],
                                                            param['rxs1']['fitmax'])

    # Calculate Rfactor for theoretical vs experimental curves
    result['rfac'] = rfactor.calculate_rfactor(expt_curve,
                                               theor_curve,
                                               param['rfac']['qmin'],
                                               param['rfac']['qmax'])

    return result

def get_curve_descriptors(data, rg_min, rg_max, rxs_min, rxs_max):
    """Calculate the Rg and Rxs1 from the input curve"""

    # Create mask to select range of Q values for Rg fitting
    rg_mask = (data[:,0] > rg_min) & (data[:,0] < rg_max)
    # Create mask to select range of Q values for Rxs fitting
    rxs_mask = (data[:,0] > rxs_min) & (data[:,0] < rxs_max)

    # Fitting is performed on:
    # Q^2 vs ln(I) for Rg
    # Q^2 vs ln(I*Q) for Rxs1
    x = data[:,0]**2
    y_rg = np.log(data[:,1])
    y_rxs = np.log(data[:,1] * data[:,0])

    rg_result = sas_curve_analysis.guinier_fit(x[rg_mask], y_rg[rg_mask], 'rg')
    rxs_result = sas_curve_analysis.guinier_fit(x[rxs_mask], y_rxs[rxs_mask], 'rxs1')

    return rg_result['r'], rxs_result['r']

args = parse_arguments()

# Read in parameters
param_file = file(args.parameter_file)
param = yaml.load(param_file)

# Create output directory and open file for summary output
if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

summary_name = os.path.join(args.output_path, args.title + '.sum')
summary_data = open(summary_name,'w')
summary_data.write("Neutron\t\t\t\t\tX-ray\n")
summary_data.write("Rg_model\tRg_curve\tRxs1_curve\tRfactor\tVolume\tRg_model\tRg_curve\tRxs1_curve\tRfactor\tVolume\n")

# Read in experimental curves and calculate Rg and Rxs
# Setup output directories for theoretical curves and sphere models

expt_name = os.path.join(args.output_path, args.title + '_expt.sum')
expt_data = open(expt_name,'w')

if args.neutron is not None:

    neut_data = rfactor.load_scatter_curve(args.neutron,
                                           param['rfac']['qmin'],
                                           param['rfac']['qmax'])

    if args.nu == 'nm':
        neut_data[:,0] = neut_data[:,0] / 10.0

    neut_rg, neut_rxs = get_curve_descriptors(neut_data,
                                              param['rg']['fitmin'],
                                              param['rg']['fitmax'],
                                              param['rxs1']['fitmin'],
                                              param['rxs1']['fitmax'])

    neut_expt_summ = "{0:7.4f}\t{1:7.4f}".format(neut_rg, neut_rxs)

    scn_path = os.path.join(args.output_path, 'neutron','curves')
    if not os.path.exists(scn_path):
        os.makedirs(scn_path)

    dry_model_path = os.path.join(args.output_path, 'neutron','models')
    if not os.path.exists(dry_model_path):
        os.makedirs(dry_model_path)

else:
    neut_expt_summ = "NA\tNA"


if args.xray is not None:

    xray_data = rfactor.load_scatter_curve(args.xray,
                                           param['rfac']['qmin'],
                                           param['rfac']['qmax'])

    if args.xu == 'nm':
        xray_data[:,0] = xray_data[:,0] / 10.0

    xray_rg, xray_rxs = get_curve_descriptors(xray_data,
                                              param['rg']['fitmin'],
                                              param['rg']['fitmax'],
                                              param['rxs1']['fitmin'],
                                              param['rxs1']['fitmax'])

    xray_expt_summ = "{0:7.4f}\t{1:7.4f}".format(xray_rg, xray_rxs)

    scx_path = os.path.join(args.output_path, 'xray','curves')
    if not os.path.exists(scx_path):
        os.makedirs(scx_path)

    wet_model_path = os.path.join(args.output_path, 'xray','models')
    if not os.path.exists(wet_model_path):
        os.makedirs(wet_model_path)

else:
    xray_expt_summ = "NA\tNA"

expt_data.write("Neutron\t\tX-ray\n")
expt_data.write("Rg\tRxs1\tRg\tRxs1\n")
expt_data.write("\t".join([neut_expt_summ, xray_expt_summ])+"\n")
expt_data.close()

# Get list of PDBs in the input directory
pdb_filter = os.path.join(args.input_path, '*.pdb')
pdb_files = glob.glob(pdb_filter)

if len(pdb_files) < 1:
    print "No PDB files found to analyze"
    sys.exit(1)

# Read in initial PDB
res_freq, atom_coords = p2s.read_pdb_atom_data(pdb_files[0])

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
    res_freq, atom_coords = p2s.read_pdb_atom_data(pdb)

    # PDB to dry sphere model
    dry_spheres, x_axis, y_axis, z_axis = p2s.create_sphere_model(atom_coords,
                                                                  cutoff,
                                                                  box_side)

    # If neutron data provided compare with curve computed from dry sphere model
    if args.neutron is not None:

        neut_theor = analyse_model(dry_spheres, neut_data, radius, param)

        # Write curve to file
        curve_file = os.path.join(scn_path, pdb_id + '.scn')
        output = open(curve_file,'w')
        calculate_curve.output_sas_curve(neut_theor['q'],neut_theor['i'], output)
        output.close()

        # Write model to file
        model_file = os.path.join(dry_model_path,  pdb_id + '.pdb')
        s2pdb.write_sphere_pdb(dry_spheres, radius, model_file)

        volume = box_side3 * len(dry_spheres)

        # Format results for output to file
        neut_summ = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\t{4:7.4f}".format(neut_theor['model_rg'],
                                                                 neut_theor['curve_rg'],
                                                                 neut_theor['curve_rxs'],
                                                                 neut_theor['rfac'][0],
                                                                 volume)
    else:
        neut_summ = "NA\tNA\tNA\tNA\tNA"

    # If x-ray data provided compare with curve computed from wet sphere model
    if args.xray is not None:
        # Hydrate model - > wet sphere model
        wet_spheres = hydrate.hydrate_model(dry_spheres,
                                            param['hydrate']['positions'],
                                            box_side,
                                            param['hydrate']['cutoff'],
                                            xaxis = x_axis,
                                            yaxis = y_axis,
                                            zaxis = z_axis)

        xray_theor = analyse_model(wet_spheres, xray_data, radius, param)

        # Write curve to file
        curve_file = os.path.join(scx_path, pdb_id + '.scx')
        output = open(curve_file,'w')
        calculate_curve.output_sas_curve(xray_theor['q'],xray_theor['i'], output)
        output.close()

        # Write model to file
        model_file = os.path.join(wet_model_path,  pdb_id + '.pdb')
        s2pdb.write_sphere_pdb(wet_spheres, radius, model_file)

        volume = box_side3 * len(wet_spheres)

        # Format results for output to file
        xray_summ = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\t{4:7.4f}".format(xray_theor['model_rg'],
                                                                 xray_theor['curve_rg'],
                                                                 xray_theor['curve_rxs'],
                                                                 xray_theor['rfac'][0],
                                                                 volume)
    else:
        xray_summ = "NA\tNA\tNA\tNA\tNA"

    # Output all summary data to file
    summary_data.write('{0:s}\t{1:s}\t{2:s}\n'.format(pdb_id,
                                                    neut_summ,
                                                    xray_summ))

summary_data.close()
