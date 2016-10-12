# -*- coding: utf-8 -*-
"""
Created on Mon May 12 12:13:28 2014

@author: dave
"""

from __future__ import absolute_import
import os
import sct.curve as curve
import sct.sphere as sphere
import sct.pdb as pdb
import sct.seq as seq
import numpy as np


def get_box_opt_input(pdb_filename, *args):
    """
    Get the target volume and atomic coordinates for box_side optimization. If
    specified use a sequence other than that of the input PDB providing the
    coordinates.
    @type  pdb_filename: string
    @param pdb_filename: Path to PDB file
    @type  *args:        list
    @param *args:        Allow entry of sequence files path and format.
                         Format is either fasta ('fas') or YAML ('yml')
    @rtype:              float, float, list
    @return:             1. Target unhydrated volume from sequence

                         2. Target hydrated volume from sequence

                         3. A list containing lists of x, y & z coordinates
                         (3 * floats)
    """

    try:
        seq_filename = args[0]
        seq_type = args[1]
    except:
        seq_filename = None
        seq_type = None

    # Read in the residues frequencies (to calculate target volume) and
    # atomic coordinates from input PDB
    res_freq, atom_coords = pdb.read_pdb_atom_data(pdb_filename)

    # If an additional sequence file was provided then use this as the sequence
    # to optimize
    if seq_type is not None:
        res_freq = seq.seq_file_to_freq(seq_filename, seq_type)

    dry_volume = seq.sum_volume(seq.all_residues,
                                res_freq, 'perkins1986a')
    # Calculate the expected volume of the hydration layer
    # Hydration is estimated to be 0.3 g water per 1 g protein
    # Bound water volume (< bulk volume) is given as a sluv parameter
    wet_volume = seq.calc_hydration_volume(res_freq) + dry_volume

    return dry_volume, wet_volume, atom_coords


def read_expt_data(
        neutron_files,
        neutron_unit,
        xray_files,
        xray_unit,
        output_path,
        param):
    """
    Read in lists of neutron and xray experimental data files. Create
    directories to hold theoretical models and curves used to compare against
    the experimental data.

    Returns the experimental data (converted in to Angstrom if required and
    including Rg and Rxs1 values).

    @type  neutron_files:  list
    @param neutron_files:  List of files containing Q, I(Q) data from neutron
                           experiments.
    @type  neutron_unit:   string
    @param neutron_unit:   Length unit used in neutron data ('a' or 'nm')
    @type  xray_files:     list
    @param xray_files:     List of files containing Q, I(Q) data from xray
                           experiments.
    @type  xray_unit:      string
    @param xray_unit:      Length unit used in xray data ('a' or 'nm')
    @type  output_path:    string
    @param output_path:    Path in which to create output directories
    @type  param:          dictionary
    @param param:          Dictionary containing standard SCT parameters
    @rtype:                list, list, dictionary
    @return:               Two lists containing dictionaries of containing the
                           experimental data ('data'), Rg ('rg'), Rxs1 ('rxs')
                           and filename ('file'). One dictionary containing the
                           paths created for the neutron models and data
                           ('dry_model' and 'scn') and the same for xrays
                           ('wet_model', 'scx').
    """
    out_paths = {}

    # Read in experimental curves and calculate Rg and Rxs
    # Setup output directories for theoretical curves and sphere models
    if neutron_files is not None:

        neut_data = curve.read_scatter_curves(
            neutron_files,
            neutron_unit,
            param)

        out_paths['scn'] = create_data_dir(output_path, 'neutron', 'curves')
        out_paths['dry_model'] = create_data_dir(
            output_path,
            'neutron',
            'models')

    else:
        neut_data = []

    if xray_files is not None:

        xray_data = curve.read_scatter_curves(xray_files, xray_unit, param)

        out_paths['scx'] = create_data_dir(output_path, 'xray', 'curves')
        out_paths['wet_model'] = create_data_dir(output_path, 'xray', 'models')

    else:
        xray_data = []

    return neut_data, xray_data, out_paths


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
    @rtype:            string
    @return:           Full path to the output directory basename/expt_type/data_type
    """

    full_path = os.path.join(basename, expt_type, data_type)
    if not os.path.exists(full_path):
        os.makedirs(full_path)

    return full_path


def analyse_theor_curve(
        theor_curve,
        expt_curves,
        param,
        chi2=False):
    """
    Take a theoretical curve and calculate the Rg and Rxs1 and the R factor or Chi^2
    values for comparisons to input experimental curves.

    @type  model:          list
    @param model:          List of lists containing x, y & z coordinates
                           (3 * floats) for each sphere in a sphere model.
    @type  expt_curves:    list
    @param expt_curves:    List of two dimensional numpy arrays containing
                           scattered vector magnitudes, q, and intensities, I,
                           from experiment.
    @type  param:          dictionary
    @param param:          Dictionary containing parameters to use when creating
                           models and analysing curves.
    @type  neutron:        boolean
    @param neutron:        Flag set if aiming to create a theoretical
                           neutron curve
    @type  chi2:           boolean
    @param chi2:           Flag set if using Chi^2 to compare curves

    @rtype:                dictionary
    @return:               Dictionary containing the following key/value pairs:
                           - model_rg: Radius of gyration calculated directly
                           from sphere model.
                           - curve_rg: Radius of gyration calculated from the
                           theoretical scattering curve derived from the sphere
                           model.
                           - curve_rxs1: Cross-section calculated from the
                           theoretical scattering curve derived from the sphere
                           model. If an rxs2 range is provided in param then
                           a 'curve_rxs2' will also be returned.
                           - rfac: R factor comparing experimental and
                           theoretical scattering curves.
    """
    result = {}

    # Rg and Rxs from theoretical curve
    if 'rxs2' in param:
        result.update(curve.get_curve_descriptors(theor_curve,
                                                  param['rg']['fitmin'],
                                                  param['rg']['fitmax'],
                                                  param['rxs1']['fitmin'],
                                                  param['rxs1']['fitmax'],
                                                  param['rxs2']['fitmin'],
                                                  param['rxs2']['fitmax']))
    else:
        result.update(curve.get_curve_descriptors(theor_curve,
                                                  param['rg']['fitmin'],
                                                  param['rg']['fitmax'],
                                                  param['rxs1']['fitmin'],
                                                  param['rxs1']['fitmax']))

    # Calculate comparison metric for theoretical vs experimental curves
    result['rfac'] = []

    for expt_curve in expt_curves:
        result['rfac'].append(curve.compare_curves(expt_curve['data'],
                                                   theor_curve,
                                                   param['rfac']['qmin'],
                                                   param['rfac']['qmax'],
                                                   chi2))

    return result


def analyse_sphere_model(
        model,
        expt_curves,
        sphere_radius,
        param,
        chi2=False,
        neutron=False):
    """
    Take a sphere model and calculate the theoretical scattering curve. The
    Rg and Rxs1 are also calculated from the curve and returned alongside the
    R factor (or Chi^2) comparing the theoretical curve to one or more
    experimental curves.

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
    @type  neutron:        boolean
    @param neutron:        Flag set if aiming to create a theoretical
                           neutron curve
    @type  chi2:           boolean
    @param chi2:           Flag set if using Chi^2 to compare curves

    @rtype:                dictionary
    @return:               Dictionary containing the following key/value pairs:
                           - model_rg: Radius of gyration calculated directly
                           from sphere model.
                           - curve_rg: Radius of gyration calculated from the
                           theoretical scattering curve derived from the sphere
                           model.
                           - curve_rxs1: Cross-section calculated from the
                           theoretical scattering curve derived from the sphere
                           model. If an rxs2 range is provided in param then
                           a 'curve_rxs2' will also be returned.
                           - rfac: R factor comparing experimental and
                           theoretical scattering curves.
    """

    result = {}

    # Calculate Rg from sphere model
    result['model_rg'] = sphere.sphere_model_rg(model, sphere_radius)
    # Create theoretical curve from model
    result['curve'] = sphere.spheres_to_sas_curve(
        model,
        sphere_radius,
        param['curve']['qmax'],
        param['curve']['npoints'],
        rbins=param['curve']['radbins'])

    # Neutron curves are usually smeared with parameters for the instrument
    # used
    if (neutron and param['curve']['smear']):

        curve.smear_sas_curve(result['curve'],
                              param['curve']['q_delta'],
                              param['curve']['wavelength'],
                              param['curve']['spread'],
                              param['curve']['divergence'])

    # Get curve derived metrics:
    #   Rg and Rxs from theoretical curve
    #   Calculate comparison metric for theoretical vs experimental curves
    result.update(
        analyse_theor_curve(
            result['curve'],
            expt_curves,
            param,
            chi2=chi2))

    return result


def sas_model_summary_output(theor_data, param):
    """
    Format summary data from modelling for output. Columns are Rg from model,
    Rg from curve, Rxs1 from curve and  volume of sphere model. An Rxs2 from
    curve column will also be produced if a 'rxs2' fit range is specified in
    param.

    @type theor_data:   dictionary, dictionary
    @param theor_data:  Dictionaries containing the following key/value pairs
                        for the dry model/neutron and wet model/xray
                        comparisons:

                        model_rg: Radius of gyration calculated directly from sphere model.

                        curve_rg: Radius of gyration calculated from the
                        theoretical scattering curve derived from the sphere
                        model.

                        curve_rxs1: Cross-section calculated from the
                        theoretical scattering curve derived from the sphere
                        model. May also contain curve_rxs2.

                        rfac: R factor comparing experimental and
                        theoretical scattering curves.

                        volume: Sphere model volume
    @type  param:       dictionary
    @param param:       Dictionary containing parameters to use when creating
                        models and analysing curves.
    @rtype:             string
    @return:            String containing the input for the sphere model Rg,
                        curve Rg, curve Rxs1 (possibly also Rxs2) and volume
    """

    if len(theor_data) > 0:
        if 'rxs2' in param:
            summ = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\t{4:7.4f}\t".format(
                theor_data['model_rg'],
                theor_data['curve_rg'],
                theor_data['curve_rxs1'],
                theor_data['curve_rxs2'],
                theor_data['volume'])
        else:
            summ = "{0:7.4f}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\t".format(
                theor_data['model_rg'],
                theor_data['curve_rg'],
                theor_data['curve_rxs1'],
                theor_data['volume'])

        for dataset in theor_data['rfac']:
            summ += "{0:7.4f}\t{1:7.4f}\t".format(dataset[0], 1.0 / dataset[1])
    else:
        if 'rxs2' in param:
            summ = "NA\tNA\tNA\tNA\tNA\tNA\tNA\t"
        else:
            summ = "NA\tNA\tNA\tNA\tNA\tNA\t"

    return summ


def create_comparison_header_text(data_files, metric):
    """
    Create the spacing for needed to provide two columns (Rfactor/Chi^2 and
    scale)     for each input experimental curve in the summary output header.
    The header     has the following format (in the final version it is tab
    separated) except without the line numbering:

    0 Path to input PDBs

    1 Neutron                                                                   X-ray

    2 Model Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * neutron curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves

    This function provides the spacing for one set of curves (either x-ray or
    neutron).

    @type data_files:   list
    @param data_files:  List of input files names containing experimental
                        scattering curves
    @rtype:             list
    @return:            List of strings containing the spacing and column
                        headings required for the three lines described above.
    """

    #if data_files is not None:
    if data_files:
        no_files = len(data_files)
#        comparison_head = ['\t'.join(['\t'] * no_files),
#                           '\t\t'.join(data_files),
#                           '\t'.join([metric + '\tscale'] * no_files)]
        comparison_head = ['\t'*(2*no_files),
                           '\t\t'.join(data_files) + '\t\t',
                           '\t'.join([metric + '\tscale'] * no_files)]
    else:
        comparison_head = ['\t\t', '\t\t', metric + '\tscale']

    return comparison_head


def write_summary_header(in_pdb, in_neut, in_xray, param, output, chi2=False):
    """
    Print the header for the summary data from sphere model analysis.
    When properties of the sphere model are references dry models are
    associated with neutrons and hydrated (wet) ones with x-rays.
    The header has the following format (in the final version it is tab
    separated) except without the line numbering (also a Rxs2 cloumn is added
    after Rxs1 if a range is supplied in param):

    0 Path to input PDBs

    1 Neutron                                                                    X-ray

    2 Model Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * neutron curves Rg_model Rg_curve Rxs1_curve Volume (Rfactor scale) * xray curves

    Each experimental curve input requires two columns - Rfactor and scale
    Where the latter is the scaling factor used to match the theoretical and
    experimental curves in the Rfactor calculation.

    @type in_pdb:    string
    @param in_pdb:   Path to directory containg input PDB files
    @type in_neut:   list
    @param in_neut:  List of input experimental data files (neutron)
    @type in_xray:   list
    @param in_xray:  List of input experimental data files (xray)
    @type param:     dictionary
    @param param:    Dictionary containing parameters to use when creating
                     models and analysing curves.
    @type output:    file
    @param output:   File to use to print output
    @type chi2:      boolean
    @param chi2:     Use Chi^2 instead of R factor as curve comparison metric
    """
    # Print header - note the inputs and print as first header row
    output.write("Input PDB path: {0:s}\n".format(in_pdb))

    if chi2:
        metric = 'Chi2'
    else:
        metric = 'Rfactor'

    # Create appropriate spacing for headings of the comparison metric
    # (Rfactor/Chi^2) related columns for both neutron and xray input.
    neut_rfact_head = create_comparison_header_text(in_neut, metric)
    xray_rfact_head = create_comparison_header_text(in_xray, metric)

    if 'rxs2' in param:
        # Header 0
        output.write("\tNeutron\t\t\t\t" + neut_rfact_head[0] + "\tX-ray\n")
        # Header 1
        output.write(
            "\t\t\t\t\t\t" +
            neut_rfact_head[1] +
            "\t\t\t\t\t" +
            xray_rfact_head[1] +
            "\n")
        # Individual column headings for header 2
        basic_head = "Rg_model\tRg_curve\tRxs1_curve\tRxs2_curve\tVolume\t"
    else:
        # Header 0
        output.write("\tNeutron\t\t\t" + neut_rfact_head[0] + "\tX-ray\n")
        # Header 1
        output.write(
            "\t\t\t\t\t" +
            neut_rfact_head[1] +
            "\t\t\t\t" +
            xray_rfact_head[1] +
            "\n")
        # Individual column headings for header 2
        basic_head = "Rg_model\tRg_curve\tRxs1_curve\tVolume\t"

    # Header 2
    col_head = "Model\t" + basic_head + \
        neut_rfact_head[2] + '\t' + basic_head + xray_rfact_head[2] + "\n"
    output.write(col_head)

    return


def perform_sas_analysis_pdb(
        pdb_file,
        neut_data,
        xray_data,
        param,
        out_paths,
        chi2=False):
    """
    Create sphere models from PDBs and use to generate theoretical scattering
    curves and compare these to experimental inputs

    @type  pdb_file:   string
    @param pdb_file:   Filename of the input PDB
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
    @type  chi2:       boolean
    @param chi2:       Use Chi^2 instead of R factor as curve comparison metric

    @rtype:            dictionary, dictionary
    @return:           Dictionaries containing the following key/value pairs
                       for the dry model/neutron and wet model/xray
                       comparisons:

                       model_rg: Radius of gyration calculated directly
                       from sphere model.

                       curve_rg: Radius of gyration calculated from the
                       theoretical scattering curve derived from the sphere
                       model.

                       curve_rxs1: Cross-section calculated from the
                       theoretical scattering curve derived from the sphere
                       model.

                       rfac: R factor comparing experimental and
                       theoretical scattering curves.

                       volume: Sphere model volume
    """

    cutoff = param['sphere']['cutoff']
    box_side = param['sphere']['boxside']
    box_side3 = param['sphere']['boxside3']
    radius = param['sphere']['radius']

    pdb_basename = os.path.basename(pdb_file)
    pdb_id = os.path.splitext(pdb_basename)[0]

    # Read PDB
    res_freq, atom_coords = pdb.read_pdb_atom_data(pdb_file)

    if len(atom_coords) > 0:

        # PDB to dry sphere model
        dry_spheres, x_axis, y_axis, z_axis = sphere.create_sphere_model(
            atom_coords, cutoff, box_side)

        # If neutron data provided compare with curve computed from dry sphere
        # model
        if len(neut_data) != 0:

            neut_theor = analyse_sphere_model(
                dry_spheres,
                neut_data,
                radius,
                param,
                chi2,
                True)

            # Write curve to file
            curve_file = os.path.join(out_paths['scn'], pdb_id + '.scn')
            curve.output_sas_curve(neut_theor['curve'], curve_file)

            # Write model to file
            model_file = os.path.join(out_paths['dry_model'], pdb_id + '.pdb')
            pdb.write_sphere_pdb(dry_spheres, radius, model_file)

            neut_theor['volume'] = box_side3 * len(dry_spheres)

        else:
            neut_theor = {}

        # If x-ray data provided compare with curve computed from wet sphere
        # model
        if len(xray_data) != 0:
            # Hydrate model - > wet sphere model
            wet_spheres = sphere.hydrate_sphere_model(
                dry_spheres,
                param['hydrate']['positions'],
                box_side,
                param['hydrate']['cutoff'],
                xaxis=x_axis,
                yaxis=y_axis,
                zaxis=z_axis)

            xray_theor = analyse_sphere_model(
                wet_spheres,
                xray_data,
                radius,
                param,
                chi2)

            # Write curve to file
            curve_file = os.path.join(out_paths['scx'], pdb_id + '.scx')
            curve.output_sas_curve(xray_theor['curve'], curve_file)

            # Write model to file
            model_file = os.path.join(out_paths['wet_model'], pdb_id + '.pdb')
            pdb.write_sphere_pdb(wet_spheres, radius, model_file)

            xray_theor['volume'] = box_side3 * len(wet_spheres)

        else:
            xray_theor = {}

        return neut_theor, xray_theor


def output_expt_summary(neut_data, xray_data, output_path, title):
    """
    Ouput the Rg and Rxs for all input experimental data files to the file:
    output_path/title_expt.sum

    @type  neut_data:    list
    @param neut_data:    List of curves read in and analysed by SCT to produce
                         scattering curve dictionaries.
    @type  xray_data:    list
    @param xray_data:    List of curves read in and analysed by SCT to produce
                         scattering curve dictionaries.
    @type  output_path:  string
    @param output_path:  Path to use for the output file
    @type  title:        string
    @param title:        Title to use as the start of the output file.
    """

    # Output summary analysis of the experimental data curves
    expt_name = os.path.join(output_path, title + '_expt.sum')
    expt_data = open(expt_name, 'w')

    all_curves = neut_data + xray_data

    if 'curve_rxs2' in all_curves[0]:
        expt_data.write("Filename\tRg\tRxs1\tRxs2\n")
    else:
        expt_data.write("Filename\tRg\tRxs1\n")

    for dat_curve in all_curves:
        if 'curve_rxs2' in dat_curve:
            expt_data.write(
                "{0:s}\t{1:7.4f}\t{2:7.4f}\t{3:7.4f}\n".format(
                    dat_curve['file'],
                    dat_curve['curve_rg'],
                    dat_curve['curve_rxs1'],
                    dat_curve['curve_rxs2']))
        else:
            expt_data.write(
                "{0:s}\t{1:7.4f}\t{2:7.4f}\n".format(
                    dat_curve['file'],
                    dat_curve['curve_rg'],
                    dat_curve['curve_rxs1']))
    expt_data.close()

    return

def valid_analyses(data):
    """
    Check that all Rg and Rxs for input experimental data are valid

    @type  data:    list
    @param data:    List of curves read in and analysed by SCT to produce
                    scattering curve dictionaries.
    @rtype:         boolean
    @return         Were all specified metrics correctly calculated
    """

    valid = True

    for data_curve in data:

        for metric in ['curve_rg','curve_rxs1','curve_rxs2']:
            if metric in data_curve:
                if np.isnan(data_curve[metric]):
                    valid = False
                    break

    return valid

def process_expt_data(
        neutron_files,
        neutron_unit,
        xray_files,
        xray_unit,
        output_path,
        output_title,
        param):
    """
    Read in lists of neutron and xray experimental data files. Create
    directories to hold theoretical models and curves used to compare against
    the experimental data. Write the summary data for each data file to a file.

    Returns the experimental data (converted in to Angstrom if required and
    including Rg and Rxs1 values).

    @type  neutron_files:  list
    @param neutron_files:  List of files containing Q, I(Q) data from neutron
                           experiments.
    @type  neutron_unit:   string
    @param neutron_unit:   Length unit used in neutron data ('a' or 'nm')
    @type  xray_files:     list
    @param xray_files:     List of files containing Q, I(Q) data from xray
                           experiments.
    @type  xray_unit:      string
    @param xray_unit:      Length unit used in xray data ('a' or 'nm')
    @type  output_path:    string
    @param output_path:    Path in which to create output directories
    @type  output_title:   string
    @param output_title:   Name to use as the start of the output file
    @type  param:          dictionary
    @param param:          Dictionary containing standard SCT parameters
    @rtype:                list, list, dictionary
    @return:               Two lists containing dictionaries of containing the
                           experimental data ('data'), Rg ('curve_rg'), Rxs1 ('curve_rxs1')
                           and filename ('file'). One dictionary containing the
                           paths created for the neutron models and data
                           ('dry_model' and 'scn') and the same for xrays
                           ('wet_model', 'scx').
    """
    # Read in experimental curves and calculate Rg and Rxs
    # Setup output directories for theoretical curves and sphere models
    neut_data, xray_data, out_paths = read_expt_data(neutron_files,
                                                     neutron_unit,
                                                     xray_files,
                                                     xray_unit,
                                                     output_path,
                                                     param)

    # Output summary analysis of the experimental data curves
    output_expt_summary(neut_data, xray_data, output_path, output_title)

    return neut_data, xray_data, out_paths
