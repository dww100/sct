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

    parser.add_argument('-q','--q_max', nargs='?', type=float,
        default=0.16, help = 'Maximum q value to make curve out to')

    parser.add_argument('-r','--radius', nargs='?', type=float,
        default=3.77, help = 'sphere radius')

    parser.add_argument('-m','--r_max', nargs='?', type=float,
        default=160, help = 'Maximum r value to use in histogram')

    parser.add_argument('-b','--n_bins', nargs='?', type=int,
        default=0.16, help = 'No. bins to use in histogram of r')

    parser.add_argument('-p','--n_points', nargs='?', type=int,
        default=100, help = 'No. points in output curve')

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


r_res = args.r_max / args.n_bins

q_delta = args.q_max / args.n_points

q = []
for n in range(0, args.n_points):
    q.append((n + 1) * q_delta)

res_freq, coords = p2s.read_pdb_atom_data(args.input_filename)
natom = len(coords)

smear = False

r_min = 0
r_max = args.r_max + r_res

hist, bin_edges, last_used = calc_r_hist(coords, r_min, r_max, args.n_bins)

sigma = []

for i in range(0, args.n_points):
    sigma_tmp = 0.0
    for j in range(0, args.no_bins):
        qr = q[i] * bin_edges[j + 1]
        sin_term = np.sin(qr)/ qr
        sigma_tmp += hist[j] * sin_term
    sigma.append(sigma_tmp)


inv_n = 1.0 / natom
inv_n2 = 1.0 / natom**2

scat = []

for i in range(0, args.n_points):
    scat.append(sphere_squared_form_factor(q[i], args.radius)
                * (inv_n + 2.0 * inv_n2 * sigma[i]))

if smear:
    scat = smear_curve(scat)

if args.output_filename != None:
    output = open(args.output_filename,'w')
else:
    output = sys.stdout

for i in range(0, args.n_points):
    output.write("{0:7.4f} {1:7.4f}\n".format(q[i], scat[i]))
