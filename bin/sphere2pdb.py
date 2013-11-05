#!/usr/bin/env python

"""Convert file containing sphere co-ordinates and radii into a PDB
Each sphere is set to be a C1 atom"""

import argparse

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""
   
    parser = argparse.ArgumentParser(
        description= 'Convert atomic pdb to sphere model\n')
        
    parser.add_argument('-i','--input_filename', nargs='?', type=str, 
        dest='input_filename', help = 'Path to the input file',
        required=True)

    parser.add_argument('-o','--output_filename', nargs='?', type=str, 
        dest='output_filename', default=None,
        help = 'Path to the output file', required=True)
   
    return parser.parse_args()

def parse_sphere_line(line):
    """Parse line from sphere file, return dict with list of 'coords' and 'radius'"""
    
    data = {}
    
    data['coords'] = [float(line[0:10]), float(line[11:20]), float(line[21:30])]
    data['radius'] = float(line[31:40])
    
    return data
    
def create_pdb_atom(res_no, res_id, atom_no, atom_type, coords, **kwargs):
    """Create ATOM line for output in a PDB file (output has a \n)"""
    
    # Columns for temperatur factor (beta) and occupancy (occ) often used to
    # store other data
    beta = kwargs.get('beta', 0.0)
    occ = kwargs.get('occ', 0.0)
    chain = kwargs.get('chain', 'A')
    
    
    line = "ATOM  {0:5d} {1:4s} {2:3s}".format(atom_no, atom_type, res_id)
    line += " {0:1s}{1:4d}    ".format(chain, res_no)
    line += "{0:8.3f}{1:8.3f}{2:8.3f}".format(coords[0],coords[1],coords[2])
    line += "{0:6.2f}{1:6.2f}\n".format(occ, beta)
    
    return line

def main():
    
    args = parse_arguments()

    out_file = open(args.output_filename, 'w')

    # We will count the number of spheres read
    count = 1

    # Read input from spheres file and output each sphere with a different 
    # atom and residue number
    with open(args.input_filename) as f:
        for line in f:
            sphere = parse_sphere_line(line)
            # As in SJP original codes use SER residue type and C1 for atom type
            # Store the sphere radius at the temperature factor (beta)
            pdb_line = create_pdb_atom(count, 'SER', count, 'C1', 
                                      sphere['coords'], beta = sphere['radius'])
            out_file.write(pdb_line)
            count += 1
        
    out_file.close()
        
if __name__ == "__main__":
    main()
    