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

import sys
import yaml
import argparse

import sct


def parse_arguments():
    """
    Parse command line arguments and ensure correct combinations present
    """

    parser = argparse.ArgumentParser(
        description='Calculate the scattering length per unit volume for a '
        'protein or glycoprotein from its amino acid and carbohydrate '
        'composition.'
    )

    parser.add_argument(
        '-i',
        '--input_file',
        nargs='?',
        type=str,
        help='Path to the input composition file',
        required=True)

    parser.add_argument(
        '-t',
        '--input_type',
        choices=[
            'fas',
            'pdb',
            'yml'],
        help='Input file format (pdb, fasta or sluv yaml)',
        default='pdb')

    parser.add_argument('-o', '--output_file', nargs='?', type=str,
                        help='Path to the output file', default=None)

    parser.add_argument(
        '-m',
        '--mode',
        choices=[
            'classic',
            'model',
            'auc',
            'project'],
        default='classic',
        help='Type of analysis to run')

    args = parser.parse_args()
    return args


def print_resid_data(resid_list, res_freq, vol_methods, out):
    """
    Print the volume, mass and scattering parameters for each passed resid
    and volume method

    @type  res_list:     list
    @param res_list:     List of residues to be output (rows)
    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @type  vol_methods:  list
    @param vol_methods:  List of the volume datasets that should be output
                         (columns).
                         Default list: cohn1943, zamyatin1972, zamyatin1984,
                         chothia1975, perkins1986a, richards1974 & perkins1986b.
    @type  out:          file object
    @param out:          File object to output the values to.
    """

    # Use params and res_vols data read in by sct.seq

    for resid in resid_list:
        line = '{0:>5s} {1:>4d}'.format(resid, res_freq[resid])

        for dataset in vol_methods:

            line = line + \
                ' {0:5.1f}'.format(sct.seq.res_vols[dataset]['residue'][resid])

        line += '   {0:>3.3g} {1:>4.3g}'.format(
            sct.seq.params['no_electron'][resid],
            sct.seq.params['mass'][resid])

        # Add scattering lengths (normal and deuterated) to line
        line += ' {0:>6.3f} {1:>6.3f}'.format(sct.seq.params['bH'][resid],
                                              sct.seq.params['bD'][resid])

        out.write(line + '\n')


def print_basic_description(res_freq, out):
    """
    Print table of volume, mass and scattering parameters for each passed
    resid and volume method

    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @type  out:          file object
    @param out:          File object to output the values to.
    """

    # Sort the names of the volume datasets for output
    # res_vols data read in by sct.seq
    vol_datasets = sorted(sct.seq.res_vols.iterkeys())

    # Title contains column headings (including sorted volume dataset names)
    title = create_volume_title("RESID  TOT", " ", vol_datasets, 'aa')
    title = title + "   MWT ELEC B(H2O) B(D2O)\n\n"
    out.write(title)

    # Print out the data, with residues grouped by type (separated by two blank
    # lines).

    print_resid_data(sct.seq.non_polar, res_freq, vol_datasets, out)

    out.write("\n\n")

    print_resid_data(sct.seq.polar, res_freq, vol_datasets, out)

    title = create_volume_title("\n          ", " ", vol_datasets, 'carb')
    out.write(title + "\n\n")

    print_resid_data(sct.seq.monosaccharides, res_freq, vol_datasets, out)

    out.write("\n\n")

    # Total number of residues
    out.write("Total " + str(sum(res_freq.itervalues())) + "\n\n")


def create_volume_title(start_string, deliminator, vol_datasets, res_type):
    """
    Create a title line for table including all chosen volume dataset names

    @type  start_string:  string
    @param start_string:
    @type  deliminator:   string
    @param deliminator:
    @type  vol_datasets:  list
    @param vol_datasets:  List of the volume datasets that should be output
                          (columns).
                          Default list: cohn1943, zamyatin1972, zamyatin1984,
                          chothia1975, perkins1986a, richards1974 &
                          perkins1986b.
    @type  res_type:      string
    @param res_type:      'aa' for amino acids or 'carb' for monosaccharides.
                          This selects the title origin - some datasets use
                          monosaccharide values from a second source.
    @rtype:               string
    @return:              Title line for volume table (including all chosen
                          volume dataset names).
    """

    title = start_string
    for dataset in vol_datasets:
        title = deliminator.join(
            [title, sct.seq.res_vols[dataset]['title'][res_type]])

    title = title + "\n"

    return title


def print_summary_data(resids, res_freq, out):
    """
    Print volume and absorption/scattering data for specified residues

    @type  resids:    list
    @param resids:    List of three letter residue codes.
    @type  res_freq:  dictionary
    @param res_freq:  Dictionary for residue type frequency. Three letter
                      residue code for the key and residue type frequency as
                      the values.
    @type  out:       file object
    @param out:       File object to use for output.
    """

    no_res = sct.seq.sum_res_no(resids, res_freq)

    if no_res == 0:
        out.write("No such residues in the provided input\n")
        return

    # compute total mass and D/H scattering lengths
    mass = sct.seq.sum_mass(resids, res_freq)
    bH_tot = sct.seq.sum_b(resids, res_freq, False)
    bD_tot = sct.seq.sum_b(resids, res_freq, True)

    # whole = is this calculation for the entire glyco protein
    # If it is we will be printing different values
    whole = False

    if len(resids) != len(sct.seq.all_residues):

        total_mass = sct.seq.sum_mass(sct.seq.all_residues, res_freq)

        frac_mass = 100 * mass / total_mass

        out.write(
            "Molecular Weight:  {0:5.0f}  Fraction of Total:  {1:3.2f}  Residues:  {2:4d}\n".format(
                mass,
                frac_mass,
                no_res))

    else:

        whole = True

        out.write("Molecular Weight:  {0:5.0f}\n".format(mass))

        abs_coeffs = sct.seq.calc_absorption_coeffs(res_freq, mass)

        out.write(
            "Absorption coefficient (280 nM):  {0:7.3f}\n".format(
                abs_coeffs[0]))
        out.write(
            "Absorption coefficient x 1.03:    {0:7.3f}\n".format(
                abs_coeffs[1]))
        out.write(
            "Absorption coefficient x 1.06:    {0:7.3f}\n".format(
                abs_coeffs[2]))

        # Calculate the contribution of the hydration layer to the
        # scattering length (b)
        vol_diff, oh_diff = sct.seq.calc_hydration_effect(res_freq)

        no_prot_res = sct.seq.sum_res_no(sct.seq.amino_acids, res_freq)
        hydra_per_res = oh_diff / no_prot_res

        hydra_delta = sct.seq.params['solvent']['vol_bound'] * oh_diff

        bH_tot_hydr = bH_tot + (sct.seq.params['solvent']['BOH'] * oh_diff)
        bD_tot_hydr = bD_tot + (sct.seq.params['solvent']['BOD'] * oh_diff)

    out.write(
        "Total b in      H2O:  {0:8.3f}  D2O:  {1:8.3f}\n".format(
            bH_tot,
            bD_tot))
    out.write("Total b on M in H2O:  {0:8.6f}  D2O:  {1:8.6f}\n".format(
        bH_tot / mass, bD_tot / mass))

    out.write(
        "Scattering density of water H2O:  {0:7.6f}  D2O:  {1:7.6f}\n".format(
            sct.seq.params['solvent']['BHHO'],
            sct.seq.params['solvent']['BDDO']))

    total_electrons = sct.seq.sum_electrons(resids, res_freq)
    out.write("Total no. electrons:\t\t{0:10.0f}\n".format(total_electrons))

    out.write("Electron density of water:\t{0:10.6f}\n".format(
        sct.seq.params['solvent']['EHHO']))

    # Sort the names of the volume datasets for output
    vol_datasets = sorted(sct.seq.res_vols.iterkeys())

    out.write(
        create_volume_title(
            "                             ",
            "   ",
            vol_datasets,
            'aa'))
    vol_line = "Volume                       "
    spec_v_line = "Specific Volume              "
    match_line = "Match Point                  "
    scat_line = "Scattering Density at MPt    "
    elect_line = "Electron Density             "

    if whole:
        hyd_vol_line = "Volume                       "
        hyd_match_line = "Match Point                  "

    # Create lines containing data for output
    for dataset in vol_datasets:
        tot_volume = sct.seq.sum_volume(resids, res_freq, dataset)
        vol_line += ' {0:7.0f}'.format(tot_volume)

        specific_volume = sct.seq.spec_volume(resids, res_freq, dataset)
        spec_v_line += ' {0:7.4f}'.format(specific_volume)

        match_point = calc_match_point(tot_volume, bH_tot, bD_tot)
        match_line += ' {0:7.2f}'.format(match_point)

        scat_density = calc_mpt_scattering_density(match_point)
        scat_line += ' {0:7.5f}'.format(scat_density)

        elect_density = sct.seq.sum_electrons(resids, res_freq) / tot_volume
        elect_line += ' {0:7.5f}'.format(elect_density)

        # Additional lines about protein hydration for the full protein
        if whole:
            hydr_vol = tot_volume + hydra_delta
            hyd_vol_line += ' {0:7.0f}'.format(hydr_vol)
            hydr_match_point = calc_match_point(
                hydr_vol,
                bH_tot_hydr,
                bD_tot_hydr)
            hyd_match_line += ' {0:7.2f}'.format(hydr_match_point)

    out.write(vol_line + '\n')
    out.write(spec_v_line + '\n')
    out.write(match_line + '\n')
    out.write(scat_line + '\n')
    out.write(elect_line + '\n')

    if whole:
        out.write(
            "********* HYDRATION OF TOTAL GLYCOPROTEIN BY OH GROUPS *********************************\n")
        out.write(
            "Difference in CHO75 and CON85 Volumes:  {0:7.0f}  Total of equivalent bound H2O: {1:7.0f}\n".format(
                vol_diff,
                oh_diff))
        out.write(
            "Average H20 per AA Residue:  {0:7.2f}\n".format(hydra_per_res))
        out.write(
            "Total b in      H2O:  {0:8.3f}  D2O:  {1:8.3f}\n".format(
                bH_tot_hydr,
                bD_tot_hydr))
        out.write(
            create_volume_title(
                "                             ",
                "   ",
                vol_datasets,
                'aa'))
        out.write(hyd_vol_line + '\n')
        out.write(hyd_match_line + '\n')


def print_exchange_data(res_freq, peptide_only, out):
    """
    Print D/H exchange data - exchange fraction + scattering lengths

    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @type  peptide_ony:  boolean
    @param peptide_ony:  Choice of whether peptide only or sidechain hydrogens
                         are considered as possible exhcanging with the solvent.
    @type  out:          file object
    @param out:          File object to use for output.
    """

    # Print a section title for the output table
    vol_datasets = sorted(sct.seq.res_vols.iterkeys())
    out.write(create_volume_title("                            ", "   ",
                                  vol_datasets, 'aa'))

    # H scattering length is a constant
    bH_tot = sct.seq.sum_b(sct.seq.all_residues, res_freq, False)

    # Increase the D exchanged percentage
    # Recalculate the D scattering length for each different percentage
    for ii in [x / 10.0 for x in range(1, 11)]:

        bD_tot = 0.0

        for resid in sct.seq.all_residues:

            # bD_tot = sum over all residues:
            # bD(all D) - (diff in b for D to H * fraction exchanged * no exchangable sites)
            # peptide only = only calculate exchange in peptide hydrogens
            if peptide_only:

                bD_tot += sct.seq.params['bD'][resid] * res_freq[resid] - (
                    sct.seq.params['no_exchange_peptide_H'][resid] * sct.seq.bDH_diff * ii * res_freq[resid])

            else:

                bD_tot += sct.seq.params['bD'][resid] * res_freq[resid] - (
                    sct.seq.params['no_exchange_H'][resid] * sct.seq.bDH_diff * ii * res_freq[resid])

            line = "Exc {0:3.1f}  Tot b: {1:5.0f}  {2:5.0f}".format(
                ii,
                bH_tot,
                bD_tot)

            for dataset in vol_datasets:

                tot_volume = sct.seq.sum_volume(sct.seq.all_residues,
                                                res_freq, dataset)

                match_point = calc_match_point(tot_volume, bH_tot, bD_tot)

                line += ' {0:7.2f}'.format(match_point)

        out.write(line + '\n')


def calc_match_point(volume, bH_tot, bD_tot):
    """
    Calculate matchpoint given a proteins volume + H/D scattering lengths

    @type  volume:  float
    @param colume:  Protein volume
    @type  bH_tot:  float
    @param bH_tot:  Total scattering length of protein
    @type  bD_tot:  float
    @param bD_tot:  Total scattering length of deuterated protein
    @rtype:         float
    @return:        Calculated match point -  ratio H2O:D2O where the scatter
                    from the solute is equal to that of the solvent
    """

    spec_bH = bH_tot / volume
    spec_bD = bD_tot / volume

    match_point = (spec_bH - sct.seq.params['solvent']['BHHO'])

    match_point = match_point / \
        (sct.seq.params['solvent']['BDDO'] - sct.seq.params['solvent']['BHHO'] - spec_bD + spec_bH)

    return match_point * 100


def calc_mpt_scattering_density(match_point):
    """
    Calculate the scattering density at a given matchpoint

    @type  match_point:  float
    @param match_point:  Scattering matchpoint
    @rtype:              float
    @return:             Calculated scattering density
    """

    scat_den = (match_point * (
        sct.seq.params['solvent']['BDDO'] - sct.seq.params['solvent']['BHHO']
    ) / 100.0) + sct.seq.params['solvent']['BHHO']

    return scat_den


def classic_output(res_freq, out):
    """
    Print volumes and scattering properties in a format similar to original sluv

    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @type  out:          file object
    @param out:          File object to use for output.
    """

    out.write("SLUV2 2013 by David W. Wright and Stephen J. Perkins\n")
    out.write(
        "Based on SLUV written by Stephen J. Perkins shortly after the dawn of time.\n\n")

    # Print frequencies and parameters for all residues
    print_basic_description(res_freq, out)

    out.write(
        "******************** TOTAL GLYCOPROTEIN ************************************************\n")
    print_summary_data(sct.seq.all_residues, res_freq, out)
    out.write(
        "******************** AA RESIDUES ONLY **************************************************\n")
    print_summary_data(sct.seq.amino_acids, res_freq, out)
    out.write(
        "******************** NONPOLAR AA RESIDUES **********************************************\n")
    print_summary_data(sct.seq.non_polar, res_freq, out)
    out.write(
        "******************** POLAR AA RESIDUES *************************************************\n")
    print_summary_data(sct.seq.polar, res_freq, out)
    out.write(
        "******************** CARBOHYDRATE RESIDUES *********************************************\n")
    print_summary_data(sct.seq.monosaccharides, res_freq, out)

    out.write(
        "******************** EXCHANGEABLE PEPTIDE HYDROGENS ************************************\n")
    print_exchange_data(res_freq, True, out)
    out.write(
        "******************** TOTAL OF EXCHANGEABLE HYDROGENS ***********************************\n")
    print_exchange_data(res_freq, False, out)


def auc_output(res_freq, out):
    """
    Print weight, absorption coeff and specific volume (Consensus)

    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @type  out:          file object
    @param out:          File object to use for output.
    """

    mass = sct.seq.sum_mass(sct.seq.all_residues, res_freq)
    out.write("Molecular Weight: {0:7.0f}\n".format(mass))

    abs_coeffs = sct.seq.calc_absorption_coeffs(res_freq, mass)
    out.write(
        "Absorption Coefficient x 1.03: {0:7.3f}\n".format(
            abs_coeffs[1]))

    specific_vol = sct.seq.spec_volume(
        sct.seq.all_residues,
        res_freq,
        'perkins1986b')
    out.write(
        "Specific Volume (Perkins 1986): {0:7.4f}\n".format(specific_vol))


def modelling_output(res_freq, out):
    """
    Print volumes used for modelling purposes

    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @type  out:          file object
    @param out:          File object to use for output.
    """

    volume = sct.seq.sum_volume(sct.seq.all_residues, res_freq, 'perkins1986a')
    wet_volume = sct.seq.calc_hydration_volume(res_freq) + volume
    out.write(
        "Volume (Perkins 1986 - Amino Acid Crystals): {0:7.0f}\n".format(volume))
    out.write("Estimated Hydrated Volume: {0:7.0f}\n".format(wet_volume))


def main():

    args = parse_arguments()

    # Get amino acid/carbohydrate occurence frequencies from file
    # Can be sluv yaml file, pdb or a fasta file
    if args.input_type == 'yml':
        protein_file = file(args.input_file, 'r')
        protein_res_freq = yaml.load(protein_file)
    elif args.input_type == 'pdb':
        try:
            protein_res_freq, coords = sct.pdb.read_pdb_atom_data(args.input_file)
        except:
            print "Input does not appear to be a valid PDB."
    elif args.input_type == 'fas':
        protein_res_freq = sct.seq.fasta_res_freq(args.input_file)

    if args.output_file is not None:
        out = open(args.output_file, 'w')
    else:
        out = sys.stdout

    # Print out the data in a format chosenfrom the command line arguments
    if args.mode == 'classic':
        classic_output(protein_res_freq, out)
    if args.mode in ('project', 'auc'):
        auc_output(protein_res_freq, out)
    if args.mode in ('project', 'model'):
        modelling_output(protein_res_freq, out)


if __name__ == "__main__":
    main()
