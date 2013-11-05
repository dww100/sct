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

def write_hydrated_model(coords, hydration_pos, radius, out):

    sep = radius * 2.0

    p2s.write_sphere_line(coords[0], coords[1], coords[2], radius, out)

    if 1 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1], coords[2], radius, out)

    if 2 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1], coords[2], radius, out)

    if 3 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1], coords[2] + sep, radius, out)

    if 4 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1], coords[2] - sep, radius, out)

    if 5 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1] + sep, coords[2], radius, out)

    if 6 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1] - sep, coords[2], radius, out)


    if 7 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1], coords[2] + sep, radius, out)

    if 8 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1], coords[2] - sep, radius, out)

    if 9 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1] + sep, coords[2] + sep, radius, out)

    if 10 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1] + sep, coords[2] - sep, radius, out)

    if 11 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1] - sep, coords[2], radius, out)

    if 12 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1] - sep, coords[2], radius, out)


    if 13 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1], coords[2] - sep, radius, out)

    if 14 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1], coords[2] + sep, radius, out)

    if 15 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1] - sep, coords[2] - sep, radius, out)

    if 16 in hydration_pos:
        p2s.write_sphere_line(coords[0], coords[1] - sep, coords[2] + sep, radius, out)

    if 17 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1] + sep, coords[2], radius, out)

    if 18 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1] + sep, coords[2], radius, out)


    if 19 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1] + sep, coords[2] + sep, radius, out)

    if 20 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1] - sep, coords[2] + sep, radius, out)

    if 21 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1] + sep, coords[2] - sep, radius, out)

    if 22 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1] - sep, coords[2] - sep, radius, out)

    if 23 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1] + sep, coords[2] - sep, radius, out)

    if 24 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1] - sep, coords[2] + sep, radius, out)


    if 25 in hydration_pos:
        p2s.write_sphere_line(coords[0] - sep, coords[1] + sep, coords[2] + sep, radius, out)

    if 26 in hydration_pos:
        p2s.write_sphere_line(coords[0] + sep, coords[1] - sep, coords[2] - sep, radius, out)

def main ():

    args = parse_arguments()

    hydration_pos = xrange(0, args.hydration_no)

    out = open(args.output_filename, 'w')

    with open(args.input_filename) as f:
        for line in f:
            sphere = parse_sphere_line(line)
            write_hydrated_model(sphere['coords'], hydration_pos, sphere['radius'], out)

    out.close()

if __name__ == "__main__":
    main()
