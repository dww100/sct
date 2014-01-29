#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Automatic optimization of the box side used to create sphere models from
atomistic PDBs in SCT tools.
The target is the theoretical volume of the PDB structure calculated from sluv2
Takes a PDB, a range of box sizes to try, a cutoff for the number of atoms
per box and a percentage tolerance for the optimization as inputs.

Within the Perkins lab this replaces the do_cubeside script.
"""

import argparse
import sys

import sct

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description= 'Optimize the box size used in creation of sphere models\n'
        )

    parser.add_argument('-i','--input_filename', nargs='?', type=str,
        dest='input_filename', help = 'Path to the input file',
        required=True)

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
        dest='output_filename', default=None,
        help = 'Path to the output file')

    parser.add_argument('-b','--box_range', nargs=2, type=float,
        default=[1.1,10.0], help = 'Minimum and maximum box sides to try')

    parser.add_argument('-t','--tolerance', nargs='?', type=float,
        default=0.01, help = 'Tolerance for box size determination')

    parser.add_argument('-c','--cutoff', nargs='?', type=int,
        default=4, help = 'Atom number cutoff for sphere models')

    return parser.parse_args()

def get_opt_box_side(filename, cutoff, side_min, side_max, tolerance):
    """
    Read PDB, get optimal sphere model box_side (+ deviation).
    The idea is to reproduce theoretical volume estimated from the sequence as
    in sluv.

    @type  filename:  string
    @param filename:  Path to PDB file to optimize box_side on
    @type  side_min:  float
    @param side_min:  Minimum side length to consider in the optimization.
    @type  side_max:  float
    @param side_min:  Maximum side length to consider in the optimization.
    @type  tolerance: float
    @param tolerance: Percentage tolerance before optimization is said to
                      have converged.
    """

    # Read in the residues frequencies (to calculate target volume) and
    # atomic coordinates from input PDB
    res_freq, atom_coords = sct.pdb.read_pdb_atom_data(filename)

    # Get the target volume from sluv (based on sequence)
    #volume = sct.seq.sum_volume(sct.seq.all_residues, res_freq, 'chothia1975')
    volume = sct.seq.sum_volume(sct.seq.all_residues, res_freq, 'perkins1986b')

    best_side, dev = sct.sphere.optimize_box_side(cutoff, atom_coords, volume,
                                       side_min, side_max, tolerance)

    return best_side, dev

def main ():

    args = parse_arguments()

    best_side, dev = get_opt_box_side(args.input_filename, args.cutoff,
                            args.box_range[0], args.box_range[1], args.tolerance)

    if args.output_filename == None:
        output = sys.stdout
    else:
        output = open(args.output_filename,'w')

    output.write("# Optimized box_side for {0:s}, using cutoff {1:d}\n".format(
                 args.input_filename, args.cutoff))
    output.write("# Deviation from target volume was {0:f}\n".format(dev))
    output.write("box_side: {0:f}\n".format(best_side))

    output.close()

if __name__ == "__main__":
    main()
