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
        default=0.16, help = 'Maximum q value in output curve')

    parser.add_argument('-r','--radius', nargs='?', type=float,
        default=3.77, help = 'Sphere radius')

    parser.add_argument('-b','--n_bins', nargs='?', type=int,
        default=400, help = 'No. bins to use in histogram of r')

    parser.add_argument('-p','--n_points', nargs='?', type=int,
        default=100, help = 'No. points in output curve')

    parser.add_argument('-s','--smear', action='store_true',
        help = 'Apply smearing to curve')

    parser.add_argument('-e','--spread', nargs='?', type=float,
        default=0.1, help = 'Wavelength spread used to calculate smearing')

    parser.add_argument('-w','--wavelength', nargs='?', type=float,
        default=6.0, help = 'Wavelength used to calculate smearing')

    parser.add_argument('-d','--divergence', nargs='?', type=float,
        default=0.016, help = 'Beam divergence used to calculate smearing')

    return parser.parse_args()

def sphere_squared_form_factor(q, r):

    qr = q * r
    qr6 = qr**6

    numerator = (3 * (np.sin(qr) - qr * np.cos(qr)))**2

    return numerator / qr6


def calc_r_hist(coords, no_bins):

    pair_dists = dist.pdist(coords, 'euclidean')
    hist, bin_edges = np.histogram(pair_dists, bins = no_bins)
    last_used = len([i for i in hist if i > 0])

    return hist, bin_edges, last_used

def spheres_to_sas_curve(coords, radius, q_max, n_points, **kwargs):

    n_bins = kwargs.get('rbins', 400)

    q_delta = q_max / n_points

    q = []
    for n in range(1, n_points + 1):
        q.append(n * q_delta)

    natom = len(coords)

    hist, bin_edges, last_used = calc_r_hist(coords, n_bins)

    sigma = []

    for i in range(0, n_points):
        sigma_tmp = 0.0
        for j in range(0, n_bins):
            qr = q[i] * bin_edges[j + 1]
            sin_term = np.sin(qr)/ qr
            sigma_tmp += hist[j] * sin_term
        sigma.append(sigma_tmp)

    inv_n = 1.0 / natom
    inv_n2 = 1.0 / natom**2

    scat = []

    for i in range(0, n_points):
        scat.append(sphere_squared_form_factor(q[i], radius)
                    * (inv_n + 2.0 * inv_n2 * sigma[i]))

    return q, scat

def output_sas_curve(q, i, output):

    for ndx,x in enumerate(q):
        output.write("{0:7.4f} {1:7.4f}\n".format(x, i[ndx]))

def smear_curve(i, q, q_delta, wavelength, spread, divergence):

    # Calculation of the sigma for the smearing gaussian
    # Based on CHAUVIN but 0.25/D^2 removed and 8LN2 not 2LN2
    # spread and alpha taken from Cusack JMB 1981 145, 539-541

    inv_wave_no = wavelength / (2.0 * np.pi)
    con = 4.0 * ((spread * inv_wave_no)**2 )
    bon = divergence**2
    aon = inv_wave_no * np.sqrt(8.0 * np.log(2.0))

    n_q = len(q)
    n_qq = (2 * n_q) + 101
    mid_qq = n_q + 50

    vv_break = np.ceil(mid_qq / 2.0)

    # Calculation of lookup table for sigma

    qq = np.empty(n_qq)
    sigma = np.empty(n_qq)
    f = np.empty(n_qq)

    qq[mid_qq] = 0.0
    f[mid_qq] = 1.0

    for ndx in range(1, mid_qq + 1):
        neg_ndx = mid_qq - ndx
        pos_ndx = mid_qq + ndx
        qq[neg_ndx] = -ndx * q_delta
        qq[pos_ndx] = ndx * q_delta
        sigma[pos_ndx] = np.sqrt(con * qq[pos_ndx]**2 + bon) / aon
        sigma[neg_ndx] = sigma[pos_ndx]
        if ndx <= n_q:
            f[pos_ndx] = i[ndx - 1]
            f[neg_ndx] = i[ndx - 1]
        else:
            f[pos_ndx] = i[-1]
            f[neg_ndx] = i[-1]

    sigma[mid_qq] = 2.0 * sigma[mid_qq + 1] - sigma[mid_qq + 2]

    g = np.empty(n_qq)
    asum  = np.zeros(n_qq)

    sqrt2pi = np.sqrt(2.0 * np.pi)

    for ndx1 in range(0, n_qq):
        # Calculation of Gaussian G
        for ndx2 in range(0, n_qq):
            # Note: difference from Chauvin and Cusack - ndx1 not ndx2 for aa, vv
            aa = 1.0 / (sigma[ndx1] * sqrt2pi)
            g[ndx2] = 0.0
            vv = (qq[ndx1] - qq[ndx2]) / sigma[ndx1]
            vv = vv**2
            if vv < vv_break:
                g[ndx2] = aa * np.exp(-vv / 2.0)

        for ndx2 in range(0, n_qq):
            asum[ndx1] += f[ndx2]*g[ndx2]

        asum[ndx1] *= q_delta

    for ndx in range(mid_qq + 1, mid_qq + n_q + 1):
        ndx2 = ndx - mid_qq - 1
        i[ndx2] = asum[ndx]

def main():

    args = parse_arguments()

    res_freq, coords = p2s.read_pdb_atom_data(args.input_filename)

    q, scat = spheres_to_sas_curve(coords, args.radius, args.q_max, args.n_points, rbins = args.n_bins)

    if args.smear:
        q_delta = args.q_max / args.n_points
        smear_curve(scat, q, q_delta, args.wavelength, args.spread, args.divergence)

    if args.output_filename != None:
        output = open(args.output_filename,'w')
    else:
        output = sys.stdout

    output_sas_curve(q, scat, output)

if __name__ == "__main__":
    main()
