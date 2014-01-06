#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculate the scattering length per unit volume for a protein or glycoprotein
from its amino acid and carbohydrate composition.

This is a reimplementation of the sluv tool originally created by
Stephen J. Perkins in 1979. Documentation of the methods used can be found in:

[SJP1] Perkins, S. J. (1986). Protein volumes and hydration effects: The
calculations of partial specific volumes, neutron scattering matchpoints
and 280-nm absorption coefficients for proteins and glycoproteins from
amino acid sequences. Eur. J. Biochem. 157, 169-180

[SJP2] Perkins, S. J. (2001). X-ray and neutron scattering analyses of
hydration shells: a molecular interpretation based on sequence predictions
and model fits. Biophysical Chemistry 93, 129–139

"""
# David W Wright - 21/10/2013

import sys
import yaml
import argparse
from sct_seq import *

# Load parameters into module global variables


# TODO: make path depend on an environment variable or something else sensible

# Load scattering and mass parameters:
# bH, bD, mass, no_electron, no_exchange_H, no_exchange_peptide_H,
# solvent -[BDDO, BHHO, EHHO], constants -[avagadro]
param_file = file('/usr0/usr3/people/davidw/dev/sjp_sas/share/sluv_parameters.yml', 'r')
params = yaml.load(param_file)

# Load the different volume datasets
vol_file = file('/usr0/usr3/people/davidw/dev/sjp_sas/share/aa_volumes.yml', 'r')
res_vols = yaml.load(vol_file)

bDH_diff = params['solvent']['BOD']-params['solvent']['BOH']

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description = 'Calculate the scattering length per unit volume for a '
        'protein or glycoprotein from its amino acid and carbohydrate '
        'composition.'
        )

    parser.add_argument('-i','--input_file', nargs='?', type=str,
        help = 'Path to the input composition file', required=True)

    parser.add_argument('-t','--input_type', choices = ['fas','pdb', 'yml'],
        help = 'Input file format (pdb, fasta or sluv yaml)', default = 'yml')

    parser.add_argument('-o','--output_file', nargs='?', type=str,
        help = 'Path to the output file', default=None)

    parser.add_argument('-m','--mode', choices = ['classic', 'model', 'auc', 'project'],
        default = 'classic', help = 'Type of analysis to run')

    args = parser.parse_args()
    return args

def print_resid_data(resid_list, res_freqs, vol_methods):
    """Print the volume, mass and scattering parameters for each passed resid
    and volume method"""

    # Use the module global params and res_vols data

    for resid in resid_list:
        line = '{0:>5s} {1:>4d}'.format(resid, res_freqs[resid])

        for dataset in vol_methods:
            line = line + ' {0:5.1f}'.format(res_vols[dataset]['residue'][resid])
        line = line + '   {0:>3.3g} {1:>4.3g} {2:>6.3f} {3:>6.3f}'.format(
                params['no_electron'][resid], params['mass'][resid],
                params['bH'][resid], params['bD'][resid])
        output.write(line + '\n')

def print_basic_description(resid_freqs, output):
    """Print table of volume, mass and scattering parameters for each passed
    resid and volume method"""

    # Sort the names of the volume datasets for output
    vol_datasets = sorted(res_vols.iterkeys())

    title = create_volume_title("RESID  TOT", " ", vol_datasets, 'aa')
    title = title + "   MWT ELEC B(H2O) B(D2O)\n\n"
    output.write(title)

    print_resid_data(non_polar, resid_freqs, vol_datasets)

    output.write("\n\n")

    print_resid_data(polar, resid_freqs, vol_datasets)

    title = create_volume_title("\n          ", " ", vol_datasets, 'carb')
    output.write(title + "\n\n")

    print_resid_data(monosaccharides, resid_freqs, vol_datasets)

    output.write("\n\n")
    output.write("Total " + str(sum(resid_freqs.itervalues())) + "\n\n")

def create_volume_title(start_string, deliminator, vol_datasets, res_type):

    title = start_string
    for dataset in vol_datasets:
        title = deliminator.join([title, res_vols[dataset]['title'][res_type]])

    title = title + "\n"

    return title


def print_summary_data(resids, resid_freqs, output):
    """Print volume and absorption/scattering data for specified residues"""

    no_res = sum_res_no(resids, resid_freqs)

    if no_res == 0:
        output.write("No such residues in the provided input\n")
        return

    # compute total mass and D/H scattering lengths
    mass = sum_mass(resids, resid_freqs)
    bH_tot = sum_b(resids, resid_freqs, False)
    bD_tot = sum_b(resids, resid_freqs, True)

    # whole = is this calculation for the entire glyco protein
    # If it is we will be printing different values
    whole = False

    if len(resids) != len(all_residues):

        total_mass = sum_mass(all_residues, resid_freqs)
        frac_mass = 100*mass/total_mass
        output.write("Molecular Weight:  {0:5.0f}  Fraction of Total:  {1:3.2f}  Residues:  {2:4d}\n".format(
            mass, frac_mass, no_res))

    else:

        whole = True

        output.write("Molecular Weight:  {0:5.0f}\n".format(mass))

        abs_coeffs = calc_absorption_coeffs(resid_freqs, mass)

        output.write("Absorption coefficient (280 nM):  {0:7.3f}\n".format(abs_coeffs[0]))
        output.write("Absorption coefficient x 1.03:    {0:7.3f}\n".format(abs_coeffs[1]))
        output.write("Absorption coefficient x 1.06:    {0:7.3f}\n".format(abs_coeffs[2]))

        # Calculate the contribution of the hydration layer to the
        # scattering length (b)
        vol_diff, oh_diff = calc_hydration_effect(resid_freqs)

        no_prot_res = sum_res_no(amino_acids, resid_freqs)
        hydra_per_res = oh_diff / no_prot_res

        hydra_delta = params['solvent']['vol_bound'] * oh_diff

        bH_tot_hydr = bH_tot + (params['solvent']['BOH'] * oh_diff)
        bD_tot_hydr = bD_tot + (params['solvent']['BOD'] * oh_diff)


    output.write("Total b in      H2O:  {0:8.3f}  D2O:  {1:8.3f}\n".format(bH_tot, bD_tot))
    output.write("Total b on M in H2O:  {0:8.6f}  D2O:  {1:8.6f}\n".format(
                bH_tot/mass, bD_tot/mass))

    output.write("Scattering density of water H2O:  {0:7.6f}  D2O:  {1:7.6f}\n".format(
                params['solvent']['BHHO'],params['solvent']['BDDO']))

    total_electrons = sum_electrons(resids, resid_freqs)
    output.write("Total no. electrons:\t\t{0:10.0f}\n".format(total_electrons))

    output.write("Electron density of water:\t{0:10.6f}\n".format(
                params['solvent']['EHHO']))

    # Sort the names of the volume datasets for output
    vol_datasets = sorted(res_vols.iterkeys())

    output.write(create_volume_title("                             ","   ",vol_datasets,'aa'))
    vol_line =    "Volume                       "
    spec_v_line = "Specific Volume              "
    match_line =  "Match Point                  "
    scat_line =   "Scattering Density at MPt    "
    elect_line =  "Electron Density             "

    if whole:
        hyd_vol_line = "Volume                       "
        hyd_match_line = "Match Point                  "

    # Create lines containing data for output
    for dataset in vol_datasets:
        tot_volume = sum_volume(resids, resid_freqs, dataset)
        vol_line += ' {0:7.0f}'.format(tot_volume)

        specific_volume = spec_volume(resids, resid_freqs, dataset)
        spec_v_line += ' {0:7.4f}'.format(specific_volume)

        match_point =  calc_match_point(tot_volume, bH_tot, bD_tot)
        match_line += ' {0:7.2f}'.format(match_point)

        scat_density = calc_mpt_scattering_density(match_point)
        scat_line += ' {0:7.5f}'.format(scat_density)

        elect_density = sum_electrons(resids, resid_freqs) / tot_volume
        elect_line  += ' {0:7.5f}'.format(elect_density)

        # Additional lines about protein hydration for the full protein
        if whole:
            hydr_vol = tot_volume + hydra_delta
            hyd_vol_line += ' {0:7.0f}'.format(hydr_vol)
            hydr_match_point = calc_match_point(hydr_vol, bH_tot_hydr,bD_tot_hydr)
            hyd_match_line += ' {0:7.2f}'.format(hydr_match_point)


    output.write(vol_line + '\n')
    output.write(spec_v_line + '\n')
    output.write(match_line + '\n')
    output.write(scat_line + '\n')
    output.write(elect_line + '\n')

    if whole:
        output.write("********* HYDRATION OF TOTAL GLYCOPROTEIN BY OH GROUPS *********************************\n")
        output.write("Difference in CHO75 and CON85 Volumes:  {0:7.0f}  Total of equivalent bound H2O: {1:7.0f}\n".format(vol_diff, oh_diff))
        output.write("Average H20 per AA Residue:  {0:7.2f}\n".format(hydra_per_res))
        output.write("Total b in      H2O:  {0:8.3f}  D2O:  {1:8.3f}\n".format(bH_tot_hydr, bD_tot_hydr))
        output.write(create_volume_title("                             ","   ",vol_datasets,'aa'))
        output.write(hyd_vol_line + '\n')
        output.write(hyd_match_line + '\n')

def print_exchange_data(resid_freqs, peptide_only, output):
    """Print D/H exchange data - exchange fraction + scattering lengths"""

    # Print a section title for the output table
    vol_datasets = sorted(res_vols.iterkeys())
    output.write( create_volume_title("                            ", "   ",
                              vol_datasets, 'aa'))

    # H scattering length is a constant
    bH_tot = sum_b(all_residues, resid_freqs, False)

    # Increase the D exchanged percentage
    # Recalculate the D scattering length for each different percentage
    for ii in [x/10.0 for x in range(1, 11)]:

        bD_tot = 0.0

        for resid in all_residues:


            # bD_tot = sum over all residues:
            # bD(all D) - (diff in b for D to H * fraction exchanged * no exchangable sites)
            # peptide only = only calculate exchange in peptide hydrogens
            if peptide_only:

                bD_tot += params['bD'][resid] * resid_freqs[resid] - (
                    params['no_exchange_peptide_H'][resid] * bDH_diff  * ii * resid_freqs[resid])

            else:

                bD_tot += params['bD'][resid] * resid_freqs[resid] - (
                    params['no_exchange_H'][resid] * bDH_diff  * ii * resid_freqs[resid])

            line = "Exc {0:3.1f}  Tot b: {1:5.0f}  {2:5.0f}".format(ii, bH_tot, bD_tot)

            for dataset in vol_datasets:
                tot_volume = sum_volume(all_residues, resid_freqs, dataset)
                match_point = calc_match_point(tot_volume, bH_tot, bD_tot)
                line += ' {0:7.2f}'.format(match_point)

        output.write(line + '\n')

def calc_model_wet_volume(res_freq):

    # Hydration is estimated to be 0.3 g water per 1 g protein
    protein_mass = sum_mass(all_residues, res_freq)
    water_mass = 0.3 * protein_mass
    # Mass of single water molecule = 18
    no_water = water_mass / 18
    # Bound water occupies smaller volume than in bulk
    water_vol = no_water * params['solvent']['vol_bound']

    return water_vol

def calc_hydration_effect(res_freq):
    """Calculates the effect of hydration of protein volume

    Returns:
        Difference in total hydrated and un-hydrated volumes
        Equivalent to the above in bound H2O
    """

    vol_chothia = sum_volume(all_residues, res_freq, 'chothia1975')
    vol_consensus = sum_volume(all_residues, res_freq, 'perkins1986b')
    volume_diff = vol_chothia - vol_consensus

    water_bound_diff = params['solvent']['vol_free'] - params['solvent']['vol_bound']

    oh_diff = volume_diff / water_bound_diff

    return volume_diff, oh_diff

def sum_res_no(resids, res_freq):
    """Get the total number of residues for given selection of residue names """
    no = 0
    for resid in resids:
        no += resid_freq[resid]
    return no

def sum_mass(resids, res_freq):
    """Calculate the total mass of the input residue constitution"""

    mass = 0.0
    for resid in resids:
        mass += params['mass'][resid] * res_freq[resid]

    return mass

def sum_volume(resids, res_freq, dataset):
    """Calculate the total volume of the input residue constitution"""

    volume = 0.0
    for resid in resids:
        volume += res_vols[dataset]['residue'][resid] * res_freq[resid]

    return volume

def spec_volume(resids, res_freq, dataset):
    """Calculate the specific volume of the input residue constitution"""

    volume = sum_volume(resids, res_freq, dataset)
    mass = sum_mass(resids, res_freq)

    return volume * params['constants']['avagadro'] / mass

def calc_absorption_coeffs(res_freq, mass):
    """Calculate absorption coefficients"""

    coeff = 10 * (150.0 * res_freq['CYS'] + 1340.0 * res_freq['TYR']
            + 5550.0 * res_freq['TRP']) / mass

    return [coeff, 1.03 * coeff, 1.06 * coeff]

def sum_b(resids, res_freq, heavy_water):
    """Calculate total scattering length of selected residues"""
    b_tot = 0.0
    for resid in resids:
        if heavy_water:
            b_tot += params['bD'][resid] * res_freq[resid]
        else:
            b_tot += params['bH'][resid] * res_freq[resid]

    return b_tot

def sum_electrons(resids, res_freq):
    """Calculate the total number of elections in the selected residues"""
    electrons = 0
    for resid in resids:
        electrons += params['no_electron'][resid] * res_freq[resid]

    return electrons

def calc_match_point(volume, bH_tot, bD_tot):
    """Calculate matchpoint given a proteins volume + B/H scattering lengths"""

    spec_bH = bH_tot / volume
    spec_bD = bD_tot / volume

    match_point = (spec_bH - params['solvent']['BHHO']) * 100
    match_point = match_point / (
        params['solvent']['BDDO'] - params['solvent']['BHHO'] - spec_bD + spec_bH)

    return match_point

def calc_mpt_scattering_density(match_point):
    """Calculate the scattering density at a given matchpoint"""

    scat_den = (match_point * (
        params['solvent']['BDDO'] - params['solvent']['BHHO']
        ) / 100.0) + params['solvent']['BHHO']

    return scat_den

def classic_output(res_freq, output):
    """Print volumes and scattering properties in a format similar to original sluv"""

    output.write("SLUV2 2013 by David W. Wright and Stephen J. Perkins\n")
    output.write("Based on SLUV written by Stephen J. Perkins shortly after the dawn of time.\n\n")

    # Print frequencies and parameters for all residues
    print_basic_description(res_freq)

    output.write("******************** TOTAL GLYCOPROTEIN ************************************************\n")
    print_summary_data(all_residues, res_freq, output)
    output.write("******************** AA RESIDUES ONLY **************************************************\n")
    print_summary_data(amino_acids, res_freq, output)
    output.write("******************** NONPOLAR AA RESIDUES **********************************************\n")
    print_summary_data(non_polar, res_freq, output)
    output.write("******************** POLAR AA RESIDUES *************************************************\n")
    print_summary_data(polar, res_freq, output)
    output.write("******************** CARBOHYDRATE RESIDUES *********************************************\n")
    print_summary_data(monosaccharides, res_freq, output)

    output.write("******************** EXCHANGEABLE PEPTIDE HYDROGENS ************************************\n")
    print_exchange_data(res_freq, True, output)
    output.write("******************** TOTAL OF EXCHANGEABLE HYDROGENS ***********************************\n")
    print_exchange_data(res_freq, False, output)

def auc_output(res_freq, output):
    """Print weight, absorption coeff and specific volume (Consensus)"""

    mass = sum_mass(all_residues, res_freq)
    output.write( "Molecular Weight: {0:7.0f}\n".format(mass))

    abs_coeffs = calc_absorption_coeffs(res_freq, mass)
    output.write( "Absorption Coefficient x 1.03: {0:7.3f}\n".format(abs_coeffs[1]))

    specific_vol = spec_volume(all_residues, res_freq, 'perkins1986b')
    output.write( "Specific Volume (Perkins 1986 - Consensus): {0:7.4f}\n".format(specific_vol))

def modelling_output(res_freq, output):
    """Print volumes used for modelling purposes"""

    volume = sum_volume(all_residues, res_freqs, 'chothia1975')

    output.write("Volume (Chothia 1975 - Crystal Structures): {0:7.0f}\n".format(volume))

    volume = sum_volume(all_residues, res_freqs, 'perkins1986b')

    output.write("Volume (Perkins 1986 - Consensus): {0:7.0f}\n".format(volume))

def main():

    args = parse_arguments()

    # Get amino acid/carbohydrate occurence frequencies from file
    # Can be sluv yaml file, pdb or a fasta file
    if args.input_type == 'yml':
        protein_file = file(args.input_file, 'r')
        protein_res_freq = yaml.load(protein_file)
    elif args.input_type == 'pdb':
        protein_res_freq = pdb_res_freq(args.input_filename)
    elif args.input_type == 'fas':
        protein_res_freq = fasta_res_freq(args.input_filename)

    if args.output_file != None:
        output = open(args.output_file,'w')
    else:
        output = sys.stdout

    # Print out the data in a format chosenfrom the command line arguments
    if args.mode == 'classic':
        classic_output(protein_res_freq, output)
    if args.mode in ('project','auc'):
        auc_output(protein_res_freq, output)
    if args.mode in ('project','model'):
        modelling_output(protein_res_freq, output)


if __name__ == "__main__":
    main()
