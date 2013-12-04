#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
import scipy.spatial.distance as dist
import pdb2sphere as p2s

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description= '\n')

    parser.add_argument('-i','--input_filename', nargs='?', type=str,
        dest='input_filename', help = 'Path to the input PDB file',
        required=True)

    parser.add_argument('-o','--output_filename', nargs='?', type=str,
        dest='output_filename', default=None,
        help = 'Path to the output file')

    return parser.parse_args()

def sphere_squared_form_factor(q, r):

    qr = q * r
    qr6 = qr**6

    numerator = (3 * (np.sin(qr) - qr * np.cos(qr)))**2

    return numerator / qr6


def calc_r_hist(coords, r_min, r_max, no_bins):

    pair_dists = dist.pdist(coords, 'euclidean')
    hist, bin_edges = np.histogram(pair_dists, bins = no_bins, range=(r_min, r_max))
    #hist, bin_edges = np.histogram(pair_dists, bins = no_bins)
    last_used = len([i for i in hist if i > 0])

    return hist, bin_edges, last_used

def smear_curve(curve):
    pass


args = parse_arguments()

# Faked values for testing
# r_res = r resolution = DIG in original SCT
r_res = 0.4
radius = 3.77

no_points = 100
q_max = 0.16
q_delta = q_max / no_points

q = []
for n in range(0, no_points):
    q.append((n + 1) * q_delta)

res_freq, coords = p2s.read_pdb_atom_data(args.input_filename)
natom = len(coords)

smear = False

no_bins = 401
hist, bin_edges, last_used = calc_r_hist(coords, 0, (no_bins + 1) * r_res, no_bins)

sigma = []

for i in range(0, no_points):
    sigma_tmp = 0.0
    for j in range(0, no_bins):
        qr = q[i] * bin_edges[j + 1]
        sin_term = np.sin(qr)/ qr
        sigma_tmp += hist[j] * sin_term
    sigma.append(sigma_tmp)


inv_n = 1.0 / natom
inv_n2 = 1.0 / natom**2

scat = []

for i in range(0,no_points):
    scat.append(sphere_squared_form_factor(q[i], radius)
                * (inv_n + 2.0 * inv_n2 * sigma[i]))

if smear:
    scat = smear_curve(scat)

if args.output_filename != None:
    output = open(args.output_filename,'w')
else:
    output = sys.stdout

for i in range(0,no_points):
    output.write("{0:7.4f} {1:7.4f}\n".format(q[i], scat[i]))
