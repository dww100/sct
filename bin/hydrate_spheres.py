#!/usr/bin/env python

"""Add hydration shell of spheres around an existing sphere model"""

import argparse
import sys
import pdb2sphere as p2s

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

    return parser.parse_args()

def parse_sphere_line(line):
    """Parse line from sphere file, return dict with list of 'coords' and 'radius'"""

    data = {}

    data['coords'] = [float(line[0:10]), float(line[11:20]), float(line[21:30])]
    data['radius'] = float(line[31:40])

    return data

def hydrate_sphere(coords, hydration_pos, radius):

    sep = radius * 2.0
    hyd_coords = []

    hyd_coords.append[coords]

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

def main ():

    args = parse_arguments()

    hydration_pos = xrange(0, args.hydration_no)

    out = open(args.output_filename, 'w')

    with open(args.input_filename) as f:
        for line in f:
            sphere = parse_sphere_line(line)
            wet_spheres = hydrate_sphere(sphere['coords'],
                                         hydration_pos, sphere['radius'])
            p2s.write_spheres(wet_spheres, sphere['radius'], out)

    out.close()

if __name__ == "__main__":
    main()
