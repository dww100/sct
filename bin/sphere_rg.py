#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
import hydrate_spheres

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description= 'Convert atomic pdb to sphere model\n')

    parser.add_argument('-i','--input_filename', nargs='?', type=str,
        dest='input_filename', help = 'Path to the input sphere model',
        required=True)

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
        default=None, help = 'Path to the output file')

    return parser.parse_args()

def sphere_model_rg(coords, radius):

    box_rg2 = radius**2

    no_spheres = len(coords)
    coords = np.array(coords)
    x_mean = np.mean(coords[ : , 0 ])
    y_mean = np.mean(coords[ : , 1 ])
    z_mean = np.mean(coords[ : , 2 ])

    rad_sq = 0.0

    for ii in range(0, no_spheres):
        d = (coords[ii, 0] - x_mean)**2 + (coords[ii , 1] - y_mean)**2 + (coords[ii, 2] - z_mean)**2
        rad_sq += d

    return np.sqrt(rad_sq / no_spheres + box_rg2)

def main():

    args = parse_arguments()

    coords, radius = hydrate_spheres.read_mono_spheres(args.input_filename)

    r_gyr = sphere_model_rg(coords, radius)

    print r_gyr

    if args.output_filename != None:
        output = open(args.output_filename,'w')
    else:
        output = sys.stdout

    output.write("Radius of gyration: {0:7.4f}\n".format(r_gyr))

if __name__ == "__main__":
    main()
