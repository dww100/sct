#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module that deal with sequence properties within SCT.
Loads in yaml files containing scattering parameters and residue volumes
Provides the following variables:
    polar:           list of polar amino acids (aa) - 3 letter codes
    non_polar:       list of non-polar aa - 3 letter codes
    amino_acids:     list of all aa - 3 letter codes
    monosaccharides: list of monosaccharides - 3 letter codes
    all_residues:    list combining amino_acids and monosaccharides - 3 letter codes
    aa1:             list of amino acid 1 letter codes
    aa1to3:          dictionary providing 1 to 3 letter aa code conversion

    params:          dictionary of scattering parameters read in from file
    res_vols:        dictionary of residue volumes read in from file

    bDH_diff:        difference in scattering length for heavy and light water
"""

import collections
import yaml
import os

share_dir = os.path.join(__file__.rsplit(os.sep,1)[0], 'share')

# Lists of the residues which are handled by SCT programs
polar = ['ARG','ASN','ASP','GLN','GLU','HIS','LYS','SER','THR']
non_polar = ['ALA','CYS','GLY','ILE','LEU','MET','PHE','PRO','TRP','TYR','VAL']
amino_acids = polar + non_polar
monosaccharides = ['FUC','GAL','GLC','MAN','NAG','NGA','SIA']
all_residues = amino_acids + monosaccharides

# Create dictionary to convert one letter to three letter amino acid codes
# Trick taken from the pymol wiki
aa1 = list("ARNDCQEGHILKMFPSTWYV")
aa1to3 = dict(zip(aa1,sorted(amino_acids)))

# Load parameters into module global variables

# TODO: make path depend on an environment variable or something else sensible

# Load scattering and mass parameters:
# bH, bD, mass, no_electron, no_exchange_H, no_exchange_peptide_H,
# solvent -[BDDO, BHHO, EHHO], constants -[avagadro]
param_file = file(os.path.join(share_dir,'sluv_parameters.yml'), 'r')
params = yaml.load(param_file)

# Load the different volume datasets
vol_file = file(os.path.join(share_dir,'aa_volumes.yml'), 'r')
res_vols = yaml.load(vol_file)

# Difference in scattering length from heavy to light water
bDH_diff = params['solvent']['BOD'] - params['solvent']['BOH']

def residue_freq_dict():
    """
    Create dictionary holding the frequencies of all residues

    @rtype:    dictionary
    @return:   Initialized dictionary for all residue types recognized by SCT
    (three letter residue codes are used as the key). Values for each residue
    type are set to 0.
    """

    res_freq = {}
    for res_name in all_residues:
        res_freq[res_name] = 0

    return res_freq

def parse_fasta(filename):
    """
    Parse fasta files. Adapted from:
    http://www.petercollingridge.co.uk/python-bioinformatics-tools/fasta-parser

    @type  filename: string
    @param filename: Path to fasta file to be read.
    @rtype:   list, dictionary
    @return:  1. Sequence identifiers found in order they appear in the file

              2. Dictionary where keys are sequence identifiers and values are
              sequences as strings (single letter amino acid codes used).
    """

    # We are going to store ordered sequence names and the sequences themselves
    order = []
    sequences = {}

    with open(filename) as f:
        for line in f:
            # Title line for sequences in fasta files start with >
            if line.startswith('>'):
                name = line[1:].rstrip('\n')
                name = name.replace('_', ' ')
                order.append(name)
                sequences[name] = ''
            elif len(order) > 0:
                # If we have seen a title but not in this line
                # add the contents of the line to the currently named sequence
                # Note: * = chain ending character in fasta so is removed
                sequences[name] += line.rstrip('\n').rstrip('*')

    return order, sequences

def fasta_res_freq(filename):
    """
    Read fasta file and return dictionary of residue frequencies.

    @type  filename: string
    @param filename: Path to fasta file to be read.
    @rtype:   dictionary
    @return:  Dictionary for all residue types recognized by SCT
              (three letter residue codes are used as the key). Values are the
              frequency with which each residue is found in input fasta file.
    """

    # Initialize a dictionary of all recognized residues (frequency = 0)
    res_freq = residue_freq_dict()

    seq_names, sequences = parse_fasta(filename)

    # Loop through all sequences found in the fasta file
    for seq_name in seq_names:
        freq_list = collections.Counter(sequences[seq_name])
        for aa in freq_list:
            # Convert one to three letter amino acid code
            res_id = aa1to3(aa)
            # Add sequence totals to overall tally for each residue
            res_freq[res_id] += freq_list[aa]

    return res_freq

def calc_model_wet_volume(res_freq):
    """
    Estimate hydrated volume of a protein from its amino acid composition

    @type  res_freq: dictionary
    @param res_freq: Dictionary for residue type frequency. Three letter residue
                     code for the key and residue type frequency as the values.
    @rtype:          float
    @return:         Estimated total volume of a protein containing the residue
                     frequencies input + a hydration shell of water.
    """

    # Hydration is estimated to be 0.3 g water per 1 g protein
    protein_mass = sum_mass(all_residues, res_freq)
    water_mass = 0.3 * protein_mass
    # Mass of single water molecule = 18
    no_water = water_mass / 18
    # Bound water occupies smaller volume than in bulk
    water_vol = no_water * params['solvent']['vol_bound']

    return water_vol

def calc_hydration_effect(res_freq):
    """
    Calculates the effect of hydration on protein volume

    @type  res_freq: dictionary
    @param res_freq: Dictionary for residue type frequency. Three letter residue
                     code for the key and residue type frequency as the values.
    @rtype:          float, float
    @return:         Difference in total hydrated and un-hydrated volumes.

                     Equivalent number of bound H2O to the volume difference.
    """

    # Dry volume
    vol_chothia = sum_volume(all_residues, res_freq, 'chothia1975')
    # Wet volume
    vol_consensus = sum_volume(all_residues, res_freq, 'perkins1986b')

    volume_diff = vol_chothia - vol_consensus

    water_bound_diff = params['solvent']['vol_free'] - params['solvent']['vol_bound']

    oh_diff = volume_diff / water_bound_diff

    return volume_diff, oh_diff

def sum_res_no(resids, res_freq):
    """
    Get the total number of residues for selected residue codes

    @type  resids:   list
    @param resids:   List of three letter residue codes.
    @type  res_freq: dictionary
    @param res_freq: Dictionary for residue type frequency. Three letter residue
                     code for the key and residue type frequency as the values.
    @rtype:          integer
    @return:         Total number of residues of the types given in resids.
    """

    no = 0
    for resid in resids:
        no += res_freq[resid]

    return no

def sum_mass(resids, res_freq):
    """
    Calculate the total mass of the input residue constitution

    @type  resids:   list
    @param resids:   List of three letter residue codes
    @type  res_freq: dictionary
    @param res_freq: Dictionary for residue type frequency. Three letter residue
                     code for the key and residue type frequency as the values.
    @rtype:          float
    @return:         Total mass of residues of the types given in resids (i.e.
                     sum[mass of residue type * frequency] over all input
                     residue types).
    """

    mass = 0.0
    for resid in resids:
        mass += params['mass'][resid] * res_freq[resid]

    return mass

def sum_volume(resids, res_freq, dataset):
    """
    Calculate the total volume of the input residue constitution

    @type  resids:   list
    @param resids:   List of three letter residue codes
    @type  res_freq: dictionary
    @param res_freq: Dictionary for residue type frequency. Three letter residue
                     code for the key and residue type frequency as the values.
    @type  dataset:  string
    @param dataset:  Choice from the read in residue volume datasets - default
                     list: cohn1943, zamyatin1972, zamyatin1984, chothia1975,
                     perkins1986a, richards1974 & perkins1986b.
    @rtype:          float
    @return:         Total volume of residues of the types given in resids -
                     sum(volume of residue type * frequency) over all input
                     residue types.
    """

    volume = 0.0
    for resid in resids:
        volume += res_vols[dataset]['residue'][resid] * res_freq[resid]

    return volume

def spec_volume(resids, res_freq, dataset):
    """
    Calculate the specific volume of the input residue constitution

    @type  resids:   list
    @param resids:   List of three letter residue codes
    @type  res_freq: dictionary
    @param res_freq: Dictionary for residue type frequency. Three letter residue
                     code for the key and residue type frequency as the values.
    @type  dataset:  string
    @param dataset:  Choice from the read in residue volume datasets - default
                     list: cohn1943, zamyatin1972, zamyatin1984, chothia1975,
                     perkins1986a, richards1974 & perkins1986b.
    @rtype:          float
    @return:         Specific volume of the input residue constitution -
                     sum(mass of residue type * frequency) * avagadros number /
                     sum(mass of residue type * frequency). Where sums are over
                     all input residue types.
    """

    volume = sum_volume(resids, res_freq, dataset)
    mass = sum_mass(resids, res_freq)

    return volume * params['constants']['avagadro'] / mass

def calc_absorption_coeffs(res_freq, mass):
    """
    Calculate absorption coefficients

    @type  res_freq: dictionary
    @param res_freq: Dictionary for residue type frequency. Three letter residue
                     code for the key and residue type frequency as the values.
    @type  mass:     float
    @param mass:     Total mass of the protein
    @rtype:          float
    @return:         Absorption coefficient
    """

    coeff = 10 * (150.0 * res_freq['CYS'] + 1340.0 * res_freq['TYR']
            + 5550.0 * res_freq['TRP']) / mass

    return [coeff, 1.03 * coeff, 1.06 * coeff]

def sum_b(resids, res_freq, heavy_water):
    """
    Calculate total scattering length of selected residues

    @type  resids:       list
    @param resids:       List of three letter residue codes
    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @type  heavy_water:  boolean
    @param            :  Are we calculating for deuterated system?
    @rtype:              float
    @return:             Total scattering length for the selected residues.
    """

    b_tot = 0.0

    for resid in resids:
        if heavy_water:
            b_tot += params['bD'][resid] * res_freq[resid]
        else:
            b_tot += params['bH'][resid] * res_freq[resid]

    return b_tot

def sum_electrons(resids, res_freq):
    """
    Calculate the total number of elections in the selected residues

    @type  resids:       list
    @param resids:       List of three letter residue codes
    @type  res_freq:     dictionary
    @param res_freq:     Dictionary for residue type frequency. Three letter
                         residue code for the key and residue type frequency as
                         the values.
    @rtype:              integer
    @return:             Total number of electrons in the selected residues
    """

    electrons = 0
    for resid in resids:
        electrons += params['no_electron'][resid] * res_freq[resid]

    return electrons
