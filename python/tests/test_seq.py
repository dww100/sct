# -*- coding: utf-8 -*-
"""
Tests for all functions in the sct.seq module
Created on Thu Sep 11 10:15:31 2014

@author: David W. Wright
"""

import sct
import os
import yaml
from nose.tools import *

def test_residue_lists():

    n_polar = len(sct.seq.polar)
    n_non_polar = len(sct.seq.non_polar)
    n_monosaccharides = len(sct.seq.monosaccharides)
    
    assert_equals(n_polar, 9)
    assert_equals(n_non_polar, 11)
    assert_equals(n_monosaccharides, 7)
    assert_equals(len(sct.seq.amino_acids), n_polar + n_non_polar)
    assert_equals(len(sct.seq.all_residues), 
                  n_polar + n_non_polar + n_monosaccharides)
                  
def test_residue_freq_dict():
    
    freq = sct.seq.residue_freq_dict()
    assert_equals(len(freq),len(sct.seq.all_residues))
    
    for res in freq:
        assert_in(res, sct.seq.all_residues)
        assert_equals(freq[res], 0)

def create_fasta(filename):
    
    fasta_text = '''>header1
ARNDCQEGHILKMFPSTWYV
>header2
LL
'''

    f = open(filename,'w')
    f.write(fasta_text)
    f.close()

def test_parse_fasta():
    
    create_fasta('tmp.fasta')
    
    order, seq = sct.seq.parse_fasta('tmp.fasta')

    assert_equals(order[0],'header1')
    assert_equals(seq['header1'],'ARNDCQEGHILKMFPSTWYV')
    assert_equals(order[1],'header2')
    assert_equals(seq['header2'],'LL')

    os.remove('tmp.fasta')

def check_fasta_freq(freq):
    
    for res in freq:
        if res == 'LEU':
            assert_equals(freq[res], 3)
        elif res in sct.seq.monosaccharides:
            pass
        else:
            assert_equals(freq[res], 1)

def test_fasta_res_freq():
    
    create_fasta('tmp.fasta')
    
    freq = sct.seq.fasta_res_freq('tmp.fasta')

    check_fasta_freq(freq)
    
    os.remove('tmp.fasta')
    

def create_seq():

    residues = ['ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'LYS', 'SER', 'THR', 
               'ALA', 'CYS', 'GLY', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 
               'TYR', 'VAL', 'FUC', 'GAL', 'GLC', 'MAN', 'NAG', 'NGA', 'SIA']
    freq = {}
    for res in residues:
        freq[res] = 1
        
    return freq
    

def create_seq_yaml(filename):
    
    freq = create_seq()

    with open(filename, 'w') as outfile:
        outfile.write(yaml.dump(freq, default_flow_style=True) )
    

def test_seq_file_to_freq():
    
    create_fasta('tmp.fasta')
    freq = sct.seq.seq_file_to_freq('tmp.fasta', 'fas')
    
    check_fasta_freq(freq)
    
    create_seq_yaml('tmp.yml')
    freq = sct.seq.seq_file_to_freq('tmp.yml', 'yml')
    
    for res in freq:
        assert_equals(freq[res], 1)
    
def test_sum_mass():

    freq = create_seq()    
    
    mass = sct.seq.sum_mass(sct.seq.all_residues, freq)
    
    assert_equals(mass, 3702.0)

def test_sum_res_no():

    freq = create_seq()    
    assert_equals(sct.seq.sum_res_no(sct.seq.all_residues,freq), 27)

    
def test_sum_volume():

    freq = create_seq()
    
    vol = sct.seq.sum_volume(sct.seq.all_residues,freq, 'perkins1986b')
    
    assert_equals(round(vol,1), 4323.5)
    
    
def test_spec_volume():
    
    freq = create_seq()
    
    spec_vol = sct.seq.spec_volume(sct.seq.all_residues,freq, 'perkins1986b')
    
    assert_equals(round(spec_vol,4), 0.7034)
    
def test_calc_absorption_coeffs():

    freq = create_seq()
    mass = sct.seq.sum_mass(sct.seq.all_residues, freq)
    
    c1, c2, c3 = sct.seq.calc_absorption_coeffs(freq, mass)
    
    assert_equals(round(c1,4), 19.0167)
    assert_equals(round(c2,4), 19.5873)
    assert_equals(round(c3,4), 20.1578)
    
def test_sum_b():
    
    freq = create_seq()
    
    b = sct.seq.sum_b(sct.seq.all_residues, freq, True)
    assert_equals(round(b,4), 147.0878)
    b = sct.seq.sum_b(sct.seq.all_residues, freq, False)    
    assert_equals(round(b,4), 84.6078)

    
def test_sum_electrons():
    
    freq = create_seq()
    
    no_e = sct.seq.sum_electrons(sct.seq.all_residues, freq)

    assert_equals(round(no_e,0), 1975)
    

def test_calc_hydration_volume():
    
    freq = create_seq()
    vol = sct.seq.calc_hydration_volume(freq)
    
    assert_equals(round(vol,2),1511.65)


def test_calc_hydration_effect():
    
    freq = create_seq()
    vol_diff, oh_diff = sct.seq.calc_hydration_effect(freq)
    
    assert_equals(round(vol_diff,1), 152.4)
    assert_equals(round(oh_diff,1), 28.2)

