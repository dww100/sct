#!/usr/bin/env python

"""Add hydration shell of spheres around an existing sphere model"""

import argparse
import sys
import pdb2sphere as p2s
import scipy.optimize as optimize

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""
    parser = argparse.ArgumentParser(
        description= 'Add hydration shell of spheres around an existing sphere model\n')

    parser.add_argument('-i','--input_filename', nargs='?', type=str,
        dest='input_filename', help = 'Path to the input file',
        required=True)

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
        dest='output_filename', default=None,
        help = 'Path to the output file', required=True)

    parser.add_argument('-n','--hydration_no', nargs='?', type=int,
        default=27, help = 'No. spheres to add as a hydration shell (1-27)')

    parser.add_argument('-c','--cutoff', nargs='?', type=int,
        default=5.0, help = 'Cutoff used to reduce the number of water spheres')

    return parser.parse_args()

def parse_sphere_line(line):
    """Parse line from sphere file, return dict with list of 'coords' and 'radius'"""

    data = {}

    data['coords'] = [float(line[0:10]), float(line[11:20]), float(line[21:30])]
    data['radius'] = float(line[31:40])

    return data

def hydrate_sphere(coords, hydration_pos, radius):
    """Add spheres to selected positions coordinate
    Inputs:
        * coordinates of sphere in dry model
        * list of positions to include in wet sphere model:
          Position 1 = original sphere position,
          2 to 27 positions on cube centred on original sphere
        * radius for each sphere
    Outputs
        * list containing coordinates of all spheres in wet model
    """

    sep = radius * 2.0
    hyd_coords = []

    hyd_coords.append(coords)

    if 1 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1], coords[2]])

    if 2 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1], coords[2]])

    if 3 in hydration_pos:
        hyd_coords.append([coords[0], coords[1], coords[2] + sep])

    if 4 in hydration_pos:
        hyd_coords.append([coords[0], coords[1], coords[2] - sep])

    if 5 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] + sep, coords[2]])

    if 6 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] - sep, coords[2]])


    if 7 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1], coords[2] + sep])

    if 8 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1], coords[2] - sep])

    if 9 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] + sep, coords[2] + sep])

    if 10 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] + sep, coords[2] - sep])

    if 11 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1] - sep, coords[2]])

    if 12 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1] - sep, coords[2]])


    if 13 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1], coords[2] - sep])

    if 14 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1], coords[2] + sep])

    if 15 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] - sep, coords[2] - sep])

    if 16 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] - sep, coords[2] + sep])

    if 17 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1] + sep, coords[2]])

    if 18 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1] + sep, coords[2]])


    if 19 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1] + sep, coords[2] + sep])

    if 20 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1] - sep, coords[2] + sep])

    if 21 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1] + sep, coords[2] - sep])

    if 22 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1] - sep, coords[2] - sep])

    if 23 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1] + sep, coords[2] - sep])

    if 24 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1] - sep, coords[2] + sep])


    if 25 in hydration_pos:
        hyd_coords.append([coords[0] - sep, coords[1] + sep, coords[2] + sep])

    if 26 in hydration_pos:
        hyd_coords.append([coords[0] + sep, coords[1] - sep, coords[2] - sep])

    return hyd_coords

def hydrate_spheres(coord_list, hydration_pos, radius):
    """Add hydration layer spheres to the spheres at input coordinates"""

    wet_spheres = []

    for coords in coord_list:
        wet_spheres += hydrate_sphere(coords, hydration_pos, radius)

    return wet_spheres

def read_mono_spheres(filename):
    """Read in sphere file and return coordinates and radius
    Note: assumed that we have only one radius of sphere here"""

    spheres = []

    with open(filename) as f:
        for line in f:
            sphere = parse_sphere_line(line)
            spheres.append(sphere['coords'])
            radius = sphere['radius']

    return spheres, radius

def hydrate_model(dry_spheres, hydration_no, box_side, water_cutoff, **kwargs):
    """Create hydrated model in 4 steps (returning the appropriate sphere
    coordinates):
    1. Add spheres to the dry model in positions assigned through hydration_pos
    2. Filter out excess spheres using create_sphere_model with water_cutoff
    (this cutoff should be high ~ 10-12 spheres per box)
    3. Add dry spheres back to the model. Some will have been removed alongside
    the excess water in the last step
    4. Filter out overlapping spheres using create_sphere_model with cutoff set
    to 1."""

    if ('xaxis' in kwargs) and ('yaxis' in kwargs) and ('zaxis' in kwargs):
        old_axes = True
        x_axis = kwargs['xaxis']
        y_axis = kwargs['yaxis']
        z_axis = kwargs['zaxis']
    else:
        old_axes = False

    radius = box_side / 2.0

    # List of positions to include in wet sphere model:
    # Position 0 = original sphere position,
    # Positions 1 to 26 positions on cube centred on original sphere
    hydration_pos = xrange(0, hydration_no + 1)

    # Create hydrated model
    wet_spheres = hydrate_spheres(dry_spheres, hydration_pos, radius)

    if old_axes:
        wet_spheres, x_axis, y_axis, z_axis = p2s.create_sphere_model(wet_spheres,
                                                                      water_cutoff,
                                                                      box_side,
                                                                      xaxis = x_axis,
                                                                      yaxis = y_axis,
                                                                      zaxis = z_axis)
    else:
        wet_spheres, x_axis, y_axis, z_axis = p2s.create_sphere_model(wet_spheres,
                                                                      water_cutoff,
                                                                      box_side)

    wet_spheres = wet_spheres + dry_spheres

    # Remove un-necessary/overlapping spheres from model
    # Use grid system from pdb2sphere with a cutoff of 1 for the number of
    # 'atoms' per box required to add a sphere

    if old_axes:
        wet_spheres, x_axis, y_axis, z_axis = p2s.create_sphere_model(wet_spheres,
                                                                      1,
                                                                      box_side,
                                                                      xaxis = x_axis,
                                                                      yaxis = y_axis,
                                                                      zaxis = z_axis)

    else:
        wet_spheres, x_axis, y_axis, z_axis = p2s.create_sphere_model(wet_spheres,
                                                                      1,
                                                                      box_side)

    return wet_spheres

def residual2_cut(cutoff, box_side, dry_spheres, hydration_no, targ_volume):
    """Compute squared residual of sphere model to the theroetical target volume"""

    # Create sphere model
    wet_model = hydrate_model(dry_spheres, hydration_no, box_side, cutoff)

    # Compute volume of the model
    vol_spheres = len(wet_model) * box_side**3

    # Return squared residual
    return (vol_spheres - targ_volume)**2

def optimize_cut(box_side, coords, hydration_no, target_vol, cut_min, cut_max, tolerance):
    """Get optimal box_side to reproduce target_volume (and deviation)"""

    # Format the bounds within which to run the optimization
    cut_bounds = (cut_min, cut_max)

    # Minimize the squared residuals between the sphere model and target volume
    # Uses minimize_scalar from scipy.optimize
    opt = optimize.minimize_scalar(residual2_cut, args=(box_side, coords, hydration_no, target_vol),
                                   bounds=cut_bounds, method='bounded',
                                   options={'xtol' : tolerance})

    # Return the optimized box_side and residual
    return int(round(opt['x'])), opt['fun']**0.5

def main ():

    args = parse_arguments()

    # Select number of solvent positions to occupy before filtering
    # Position 0 = original sphere position,
    # Positions 1 to 26 positions on cube centred on original sphere
    hydration_no = 26

    # Read dry sphere model from file
    dry_spheres, radius = read_mono_spheres(args.input_filename)

    # Create hydrated model
    wet_spheres = hydrate_model(dry_spheres, hydration_no, radius, args.cutoff)

    out = open(args.output_filename, 'w')
    p2s.write_spheres(wet_spheres, radius, out)
    out.close()

if __name__ == "__main__":
    main()
