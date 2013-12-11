#!/usr/bin/env python

"""Automatic optimization of the box side used to create sphere models from
atomistic PDBs in SCT tools.
The target is the theoretical volume of the PDB structure calculated from sluv2
Takes a PDB, a range of box sizes to try, a cutoff for the number of atoms
per box and a percentage tolerance for the optimization as inputs.

Within the Perkins lab this replaces the do_cubeside script."""

import argparse
import pdb2sphere as p2s
import sluv2
import sys
import scipy.optimize as optimize

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

def residual2_box(box_side, cutoff, atom_coords, targ_volume):
    """Compute squared residual of sphere model to the theroetical target volume"""

    # Create sphere model
    sphere_coords, x_axis, y_axis, z_axis = p2s.create_sphere_model(atom_coords,
                                                               cutoff, box_side)
    # Compute volume of the model
    vol_spheres = len(sphere_coords) * box_side**3

    # Return squared residual
    return (vol_spheres - targ_volume)**2

def optimize_side(cutoff, coords, target_vol, side_min, side_max, tolerance):
    """Get optimal box_side to reproduce target_volume (and deviation)"""

    # Format the bounds within which to run the optimization
    side_bounds = (side_min, side_max)

    # Minimize the squared residuals between the sphere model ang target volume
    # Uses minimize_scalar from scipy.optimize
    opt = optimize.minimize_scalar(residual2_box, args=(cutoff, coords, target_vol),
                             bounds=side_bounds, method='bounded',
                             options={'xtol' : tolerance})

    # Return the optimized box_side and residual
    return opt['x'], opt['fun']**0.5

def get_opt_side(filename, cutoff, side_min, side_max, tolerance):
    """Read PDB, get optimal sphere model box_side (+ deviation)
    The idea is to reproduce theoretical volume from sluv"""

    # Read in the residues frequencies (to calculate target volume) and
    # atomic coordinates from input PDB
    res_freq, atom_coords = p2s.read_pdb_atom_data(filename)

    # Get the target volume from sluv (based on sequence)
    volume = sluv2.sum_volume(sluv2.all_residues, res_freq, 'chothia1975')

    best_side, dev = optimize_side(cutoff, atom_coords, volume,
                              side_min, side_max, tolerance)

    return best_side, dev

def main ():

    args = parse_arguments()

    best_side, dev = get_opt_side(args.input_filename, args.cutoff,
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
