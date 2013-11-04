#!/usr/bin/env python
import argparse
import collections
import sys

polar = ['ARG','ASN','ASP','GLN','GLU','HIS','LYS','SER','THR']
non_polar = ['ALA','CYS','GLY','ILE','LEU','MET','PHE','PRO','TRP','TYR','VAL']
amino_acids = polar + non_polar
monosaccharides = ['FUC','GAL','GLC','MAN','NAG','NGA','SIA']
all_residues = amino_acids + monosaccharides

# Create dictionary to convert one letter to three letter amino acid codes
# Trick taken from the pymol wiki
aa1 = list("ARNDCQEGHILKMFPSTWYV")
aa1to3 = dict(zip(aa1,sorted(amino_acids)))

def residue_freq_dict():

    # Create dictionary holding the frequencies of each amino acid
    res_freq = {}
    for res_name in all_residues:
        res_freq[res_name] = 0
        
    return res_freq

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""
   
    parser = argparse.ArgumentParser(
        description= 'Convert fasta file to yaml amino acid frequency file\n')   
    parser.add_argument('-i','--infile', nargs='?', type=str, 
        dest='inputFilename', help = 'Path to the input fasta file',
        required=True)
    parser.add_argument('-o','--outfile', nargs='?', type=str, 
        dest='outputFilename', help = 'Path to the output file', 
        required=True)
    
    return parser.parse_args()

def parse_fasta(filename):
    """Parse fasta files. Adapted from:
    http://www.petercollingridge.co.uk/python-bioinformatics-tools/fasta-parser
    """
    order = []
    sequences = {}
    
    with open(filename) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].rstrip('\n')
                name = name.replace('_', ' ')
                order.append(name)
                sequences[name] = ''
            elif len(order) > 0:
                sequences[name] += line.rstrip('\n').rstrip('*')

    return order, sequences
    
def pdb_res_line_parse(line):
    
    data = {}
    
    record = line[0:6].strip()
    if record in ['ATOM','HETATM']:
        res_id = line[17:20].strip()
        if res_id in all_residues:
            data['atom_no'] = int(line[6:11])
            data['atom_name'] = line[12:16].strip()
            data['res_id'] = res_id
            data['chain'] = line[21:22].strip()
            data['res_no'] = int(line[22:26])
            data['x'] = float(line[30:38])
            data['y'] = float(line[38:46])
            data['z'] = float(line[46:54])
    
    return data
    
def pdb_res_freq(filename):
            
    res_freq = residue_freq_dict()

    last_res_no = 0    
            
    with open(filename) as f:
        for line in f:
            residue_data = pdb_res_line_parse(line)
            if 'res_no' in residue_data:
                print residue_data
                if residue_data['res_no'] != last_res_no:
                    res_freq[residue_data['res_id']] += 1
                    last_res_no = residue_data['res_no']
                
    return res_freq
    
def fasta_res_freq(filename):
    
    res_freq = residue_freq_dict()
    
    seq_names, sequences = parse_fasta(filename)

    for seq_name in seq_names:
        freq_list = collections.Counter(sequences[seq_name])
        for aa in freq_list:
            res_id = aa1to3(aa)
            res_freq[res_id] += freq_list[aa]
            
    return res_freq
    
def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""
   
    parser = argparse.ArgumentParser(
        description= 'Convert fasta or pdb file to yaml residue frequency file\n')   
    parser.add_argument('-i','--infile', nargs='?', type=str, 
        dest='input_filename', help = 'Path to the input fasta/pdb file',
        required=True)
    parser.add_argument('-p', '--pdb', action='store_true', default=False,
        help = 'Flag to indicate input is a PDB rather than fasta file')
    parser.add_argument('-o','--outfile', nargs='?', type=str, 
        dest='output_filename', help = 'Path to the output file')
    
    return parser.parse_args()
    
def main():

    args = parse_arguments()
    if args.pdb:
        res_freq = pdb_res_freq(args.input_filename)
    else:
        res_freq = fasta_res_freq(args.input_filename)

    if args.output_file != None:
        sys.stdout = open(args.output_file,'w')    

    for res_name in amino_acids:
        print res_name + ': ' + str(res_freq[res_name]) + '\n'
    for res_name in monosaccharides:
        print res_name + ': ' + str(res_freq[res_name]) + '\n'

if __name__ == "__main__":
    main()