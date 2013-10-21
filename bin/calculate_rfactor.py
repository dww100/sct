#!/bin/env python

# David W. Wright - 08 October 2013
"""Script to compute the R factor of a comparison between two scattering curves.
This will usually be one calculated and one experimental data set.

R factor is used by analogy with crystallography where:
R = sum (abs(F_expt - F_calc)) / sum (abs(F_expt))
"""

import argparse
import numpy
from sjp_util import sjp_util

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description = 'Compute R factor from a comparison of two scattering curves.'
        )
    parser.add_argument(
        '-c','--calc', nargs='?', type=str, dest='calc_curve', 
        help = 'Path to the input calculated curve', required=True)
    parser.add_argument(
        '-e','--expt', nargs='?', type=str, dest='expt_curve', 
        help = 'Path to the input experimental curve', required=True)
    parser.add_argument(
        '-q','--qrange', nargs=2, type=float, 
        help = 'Minimum and maximum Q values used in curve fitting', 
        required=True)
    parser.add_argument(
        '-o','--outfile', nargs='?', type=str, dest='out_file', 
        help = 'Path to the output file', required=False)
    parser.add_argument(
        '--header', action='store_true', 
        help = 'Output a header (no data analysis performed)')   
        
    args = parser.parse_args()
    return args

def load_scatter_curve(filename, q_min, q_max):
    """Load I, Q data from file and filter for values between q_min and q_max"""
    
    scatter_data = numpy.loadtxt(args.calc_curve)
    qrange_mask = (scatter_data[:,0] >= q_min) & (scatter_data[:,0] <= q_max)
    
    return scatter_data[qrange_mask]
    
def match_scatter_curves(target_data, source_data):
    """Match I values from one data set to the Q values of another.
    
    Input data is two I vs Q scattering curves.
    Output is I values from the source data set matched to the Q values of 
    the target set.
    """

    # Create list of calculated I values matched to the nearest experimental Q value
    # Remember that the arrays in python start at 0, those in Fortran at 1
    last_calc = len( calc_data )
    last_expt = len( expt_data )

    # Initialize array to hold the calculated I values matched to experimental q values
    matched_calc_I = numpy.zeros( last_expt, dtype=float )
    
    # Use the old fortran routine to match the data sets by q value
    # matched_no is the number of datapoints which contain matched data
    matched_no = sjp_util.qrange_match(expt_data[:,0], calc_data[:,0], 
                                        calc_data[:,1], last_expt, 
                                        last_calc, matched_calc_I)
    
    return matched_no, matched_calc_I

# Interpret command line arguments
args = parse_arguments()

if (not args.header) :

    q_min = args.qrange[0]
    q_max = args.qrange[1]

    # Load the data from the two input files
    calc_data = numpy.loadtxt(args.calc_curve)
    expt_data = numpy.loadtxt(args.expt_curve)

    # Filter inputs to retain only values for Q values between q_min and q_max
    idx=(calc_data[:,0] >= q_min) & (calc_data[:,0] <= q_max)
    calc_data = calc_data[idx]

    idx=(expt_data[:,0] >= q_min) & (expt_data[:,0] <= q_max)
    expt_data = expt_data[idx]

    # Create list of calculated I values matched to nearest experimental Q value
    # Remember that the arrays in python start at 0, those in our Fortran at 1
    last_calc = len(calc_data)
    last_expt = len(expt_data)

    matched_calc_I = numpy.zeros(last_expt,dtype=float)
    matched_no = sjp_util.qrange_match(expt_data[:,0], calc_data[:,0], 
        calc_data[:,1], last_expt, last_calc, matched_calc_I)

    # Get average I for experimental and calculated values over matched Q range

    expt_avg = numpy.mean( expt_data[0:matched_no - 1, 1] )
    calc_avg = numpy.mean( matched_calc_I[0:matched_no - 1] )

    # Initial guess of the concentration:
    # ratio of experimental and calculated average intensities
    con = expt_avg / calc_avg

    r_factor = sjp_util.calc_rfactor( expt_data[:,0], expt_data[:,1], 
        matched_calc_I, matched_no, q_min, q_max, con, False)

    # 1/con is the scaling factor needed to multiply experimental I values 
    # to compare with calculated data
    scale = 1.0 / con

    output_data = '{0:s}\t{1:s}\t{2:0.5f}\t{3:0.5f}\t{4:0.5f}\t{5:0.5f}'.format(
        args.calc_curve, args.expt_curve, q_min, q_max, scale, r_factor)

    if args.out_file == None:
        print output_data
    else:
        with open(args.out_file, "a") as fle:
            fle.write(output_data)
else:

    header = 'Calculated file\tExperimental file\tQ min\tQ max\tScaling factor\tR factor'
    if args.out_file == None:
        print header
    else:
        with open(args.out_file, "a") as fle:
            fle.write(header)