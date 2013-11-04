#!/usr/bin/env python
import argparse
import collections
import sys

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

def residue_freq_dict():
    """Create dictionary holding the frequencies of all residues"""
    
    res_freq = {}
    for res_name in all_residues:
        res_freq[res_name] = 0
        
    return res_freq

def parse_fasta(filename):
    """Parse fasta files. Adapted from:
    http://www.petercollingridge.co.uk/python-bioinformatics-tools/fasta-parser
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
    
def pdb_res_line_parse(line):
    """Parse a single line from a PDB file into a dictionary according to 
    the standard PDB column definitions
    
    Output: dictionary containing:
    atom_no, atom_name, res_id, chain, res_no, coords
    coords is a list of the x, y & z coordinates"""
    
    data = {}

    # Only interested in ATOM and HETATM lines not REMARKs etc.    
    # The type of record on the line is determined by the first 6 characters
    record = line[0:6].strip()    
    if record in ['ATOM','HETATM']:

        # Split data on the line according to column definitions for PDB files
        # Ignore residues that we can't handle in SCT
        # TODO: This should perhaps give a warning        
        
        # res_id is a three letter residue code
        res_id = line[17:20].strip()                        
        if res_id in all_residues:
            data['atom_no'] = int(line[6:11])
            data['atom_name'] = line[12:16].strip()
            data['res_id'] = res_id
            data['chain'] = line[21:22].strip()
            data['res_no'] = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])          
            data['coords'] = [x, y, z]
    
    return data
    
def pdb_res_freq(filename):
    """ Read PDB file and return dictionary of residue frequencies"""
            
    # Initialize a dictionary of all recognized residues (frequency = 0)
    res_freq = residue_freq_dict()

    # Track the last residues seen
    last_res_no = 0    
            
    with open(filename) as f:
        for line in f:
            residue_data = pdb_res_line_parse(line)
            if 'res_no' in residue_data:                
                # If there is a change in residue increment appropriate tally
                if residue_data['res_no'] != last_res_no:
                    res_freq[residue_data['res_id']] += 1
                    last_res_no = residue_data['res_no']
                
    return res_freq
    
def fasta_res_freq(filename):
    """ Read fasta file and return dictionary of residue frequencies"""
    
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
    
def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""
   
    parser = argparse.ArgumentParser(
        description= 'Convert fasta or pdb file to yaml residue frequency file\n')   
    parser.add_argument('-i','--input_filename', nargs='?', type=str, 
        dest='input_filename', help = 'Path to the input fasta/pdb file',
        required=True)
    parser.add_argument('-p', '--pdb', action='store_true', default=False,
        help = 'Flag to indicate input is a PDB rather than fasta file')
    parser.add_argument('-o','--output_filename', nargs='?', type=str, 
        dest='output_filename', help = 'Path to the output file')
    
    return parser.parse_args()
    
def main():

    args = parse_arguments()
    
    # Read input
    # Input should either be PDB or fasta format - chosen by input flag
    if args.pdb:
        res_freq = pdb_res_freq(args.input_filename)
    else:
        res_freq = fasta_res_freq(args.input_filename)

    if args.output_file != None:
        sys.stdout = open(args.output_file,'w')    

    # Output pdb frequencies in a yaml style format that can be read by SLUV2
    for res_name in amino_acids:
        print res_name + ': ' + str(res_freq[res_name]) + '\n'
    for res_name in monosaccharides:
        print res_name + ': ' + str(res_freq[res_name]) + '\n'

if __name__ == "__main__":
    main()