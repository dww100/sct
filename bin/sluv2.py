#!/bin/env python
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

import yaml
import argparse

# Load parameters into module global variables


# TODO: make path depend on an environment variable or something else sensible

# Load scattering and mass parameters:
# bH, bD, mass, no_electron, no_exchange_H, no_exchange_peptide_H, 
# solvent -[BDDO, BHHO, EHHO], constants -[avagadro]
param_file = file('../share/sluv_parameters.yml', 'r')
params = yaml.load(param_file)

# Load the different volume datasets
vol_file = file('../share/aa_volumes.yml', 'r')
res_vols = yaml.load(vol_file)

polar = ['ARG','ASN','ASP','GLN','GLU','HIS','LYS','SER','THR']
non_polar = ['ALA','CYS','GLY','ILE','LEU','MET','PHE','PRO','TRP','TYR','VAL']
amino_acids = polar + non_polar
carbohydrate = ['FUC','GAL','GLC','MAN','NAG','NGA','SIA']
all_residues = amino_acids + carbohydrate

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description = 'Calculate the scattering length per unit volume for a '
        'protein or glycoprotein from its amino acid and carbohydrate '
        'composition.'
        )

    parser.add_argument('-i','--input_file', nargs='?', type=str,  
        help = 'Path to the input composition file', required=True)
        
    parser.add_argument('-o','--output_file', nargs='?', type=str,  
        help = 'Path to the output file', default=None)
   
    parser.add_argument('-m','--mode', choices = ['classic', 'model', 'auc'], 
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
        print line

def print_basic_description(resid_freqs):
    """Print table of volume, mass and scattering parameters for each passed 
    resid and volume method""" 

    # Sort the names of the volume datasets for output
    vol_datasets = sorted(res_vols.iterkeys())    
            
    title = create_volume_title("\nRESID  TOT", " ", vol_datasets, 'aa')
    title = title + "   MWT ELEC B(H2O) B(D2O)\n"
    print title

    print_resid_data(non_polar, resid_freqs, vol_datasets)

    print "\n"

    print_resid_data(polar, resid_freqs, vol_datasets)

    title = create_volume_title("\n          ", " ", vol_datasets, 'carb')
    print title + "\n"

    print_resid_data(carbohydrate, resid_freqs, vol_datasets)

    print "\n"    
    print "Total " + str(sum(resid_freqs.itervalues())) + "\n"

def create_volume_title(start_string, deliminator, vol_datasets, res_type):
    
    title = start_string
    for dataset in vol_datasets:
        title = deliminator.join([title, res_vols[dataset]['title'][res_type]])
        
    return title
    
           
def print_summary_data(resids, resid_freqs):
    """Print all data for specified part of the system
    
    Once I figure out exactly what this will do this will need updating"""

    no_res = sum_res_no(resids, resid_freqs)

    if no_res == 0:
        print "No such residues in the provided input"
        return
            
    mass = sum_mass(resids, resid_freqs)
    
    whole = False
    
    if len(resids) != len(all_residues):
        
        total_mass = sum_mass(all_residues, resid_freqs)
        frac_mass = 100*mass/total_mass
        print "Molecular Weight:  {0:5.0f}  Fraction of Total:  {1:3.2f}  Residues:  {2:4d}".format(
            mass, frac_mass, no_res)
            
    else:
        
        whole = True
        
        print "Molecular Weight:  {0:5.0f}".format(mass)
            
        abs_coeffs = calc_absorption_coeffs(resid_freqs, mass)
    
        print "Absorption coefficient (280 nM):  {0:7.3f}".format(abs_coeffs[0])
        print "Absorption coefficient x 1.03:    {0:7.3f}".format(abs_coeffs[1])
        print "Absorption coefficient x 1.06:    {0:7.3f}".format(abs_coeffs[2])
        
    bH_tot = sum_b(resids, resid_freqs, False)        
    bD_tot = sum_b(resids, resid_freqs, True)
    
    if whole:
        vol_diff, oh_diff = calc_hydration_effect(resid_freqs)

        no_prot_res = sum_res_no(amino_acids, resid_freqs)        
        hydra_per_res = oh_diff / no_prot_res
        
        hydra_delta = params['solvent']['vol_bound'] * oh_diff
        
        bH_tot_hydr = bH_tot + (params['solvent']['BOH'] * oh_diff)
        bD_tot_hydr = bD_tot + (params['solvent']['BOD'] * oh_diff)
    
    
    print "Total b in      H2O:  {0:8.3f}  D2O:  {1:8.3f}".format(bH_tot, bD_tot)
    print "Total b on M in H2O:  {0:8.6f}  D2O:  {1:8.6f}".format(
                bH_tot/mass, bD_tot/mass)
    
    print "Scattering density of water H2O:  {0:7.6f}  D2O:  {1:7.6f}".format(
                params['solvent']['BHHO'],params['solvent']['BDDO'])
                            
    total_electrons = sum_electrons(resids, resid_freqs)
    print "Total no. electrons:\t\t{0:10.0f}".format(total_electrons)

    print "Electron density of water:\t{0:10.6f}".format(
                params['solvent']['EHHO'])
       
    # Sort the names of the volume datasets for output
    vol_datasets = sorted(res_vols.iterkeys())

    print create_volume_title("                             ","   ",vol_datasets,'aa')
    vol_line =    "Volume                       "
    spec_v_line = "Specific Volume              "
    match_line =  "Match Point                  "
    scat_line =   "Scattering Density at MPt    "
    elect_line =  "Electron Density             "
    
    if whole:
        hyd_vol_line = "Volume                       "
        hyd_match_line = "Match Point                  "
        
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
            
            if whole:
                hydr_vol = tot_volume + hydra_delta
                hyd_vol_line += ' {0:7.0f}'.format(hydr_vol)
                hydr_match_point = calc_match_point(hydr_vol, bH_tot_hydr,bD_tot_hydr)
                hyd_match_line += ' {0:7.2f}'.format(hydr_match_point)
            
            
    print vol_line
    print spec_v_line
    print match_line
    print scat_line
    print elect_line
    
    if whole:
        print "********* HYDRATION OF TOTAL GLYCOPROTEIN BY OH GROUPS *********************************"
        print "Difference in CHO75 and CON85 Volumes:  {0:7.0f}  Total of equivalent bound H2O: {1:7.0f}".format(vol_diff, oh_diff)
        print "Average H20 per AA Residue:  {0:7.2f}".format(hydra_per_res)
        print "Total b in      H2O:  {0:8.3f}  D2O:  {1:8.3f}".format(bH_tot_hydr, bD_tot_hydr)
        print create_volume_title("                             ","   ",vol_datasets,'aa')
        print hyd_vol_line
        print hyd_match_line

def calc_hydration_effect(resid_freqs):

    vol_chothia = sum_volume(all_residues, resid_freqs, 'chothia1975')
    vol_consensus = sum_volume(all_residues, resid_freqs, 'perkins1986b')
    volume_diff = vol_chothia - vol_consensus
    
    water_bound_diff = params['solvent']['vol_free'] - params['solvent']['vol_bound']
    
    oh_diff = volume_diff / water_bound_diff

    return volume_diff, oh_diff
                                
def sum_res_no(resids, resid_freqs):
    """Get the total number of residues for given selection of residue names """
    no = 0
    for resid in resids:
        no += resid_freqs[resid]
    return no
      
def sum_mass(resids, resid_freqs):
    """Calculate the total mass of the input residue constitution"""

    mass = 0.0
    for resid in resids:
        mass += params['mass'][resid] * resid_freqs[resid]
        
    return mass
    
def sum_volume(resids, resid_freqs, dataset):
    """Calculate the total volume of the input residue constitution"""
    
    volume = 0.0
    for resid in resids:
        volume += res_vols[dataset]['residue'][resid] * resid_freqs[resid]    

    return volume

def spec_volume(resids, resid_freqs, dataset):
    """Calculate the specific volume of the input residue constitution"""
    
    volume = sum_volume(resids, resid_freqs, dataset)
    mass = sum_mass(resids, resid_freqs)
    
    return volume * params['constants']['avagadro'] / mass

def calc_absorption_coeffs(resid_freqs, mass):
    """Calculate absorption coefficients"""
    
    coeff = 10 * (150.0 * resid_freqs['CYS'] + 1340.0 * resid_freqs['TYR'] 
            + 5550.0 * resid_freqs['TRP']) / mass

    return [coeff, 1.03 * coeff, 1.06 * coeff]

def sum_b(resids, resid_freqs, heavy_water):
    """Calculate total scattering length of selected residues"""
    b_tot = 0.0
    for resid in resids:
        if heavy_water:
            b_tot += params['bD'][resid] * resid_freqs[resid]
        else:
            b_tot += params['bH'][resid] * resid_freqs[resid]
    
    return b_tot

def sum_electrons(resids, resid_freqs):
    """Calculate the total number of elections in the selected residues"""
    electrons = 0
    for resid in resids:
        electrons += params['no_electron'][resid] * resid_freqs[resid]    

    return electrons

def calc_match_point(volume, bH_tot, bD_tot):
#def calc_match_point(resids, resid_freqs, vol_dataset):

#    volume = sum_volume(resids, resid_freqs, vol_dataset)
#    bH_tot = sum_b(resids, resid_freqs, False)
#    bD_tot = sum_b(resids, resid_freqs, True)
    
    spec_bH = bH_tot / volume
    spec_bD = bD_tot / volume
    
    match_point = (spec_bH - params['solvent']['BHHO']) * 100
    match_point = match_point / (
        params['solvent']['BDDO'] - params['solvent']['BHHO'] - spec_bD + spec_bH)
                
    return match_point

def calc_mpt_scattering_density(match_point):
    
    scat_den = (match_point * (
        params['solvent']['BDDO'] - params['solvent']['BHHO']
        ) / 100.0) + params['solvent']['BHHO']
        
    return scat_den
      
def main():

    args = parse_arguments()
        
    # Get amino acid/carbohydrate occurence ferquencies from file
    protein_file = file(args.input_file, 'r')
    protein_res_freq = yaml.load(protein_file)

    # Print frequencies and parameters for all residues
    print_basic_description(protein_res_freq)
    
    print "******************** TOTAL GLYCOPROTEIN ************************************************" 
    print_summary_data(all_residues, protein_res_freq)
    print "********************* AA RESIDUES ONLY *************************************************"
    print_summary_data(amino_acids, protein_res_freq)
    print "******************* NONPOLAR AA RESIDUES ***********************************************" 
    print_summary_data(non_polar, protein_res_freq)
    print "******************** POLAR AA RESIDUES *************************************************" 
    print_summary_data(polar, protein_res_freq)
    print "****************** CARBOHYDRATE RESIDUES ***********************************************"
    print_summary_data(carbohydrate, protein_res_freq)
         

#    print "************* EXCHANGEABLE PEPTIDE HYDROGENS **************"


if __name__ == "__main__":
    main()

