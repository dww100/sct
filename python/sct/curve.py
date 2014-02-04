#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Scattering curve related functions used in the SCT suite of programs
"""

# Copyright 2014 David W. Wright

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import yaml
import numpy as np
import matplotlib as plt

from sjp_util import sjp_util

def load_scatter_curve(filename, q_min, q_max):
    """
    Load magnitude of scattering vector, q, and intensity, I, data from file and
    filter for values between q_min and q_max.
    Note: Assuming that a/ there is no header and b/ q and I are the first two
    columns here.

    @type  filename:  string
    @param filename:  Path to file containing scattering curve (q, I data).
    @type  q_min:     float
    @param q_min:     Minimum value of the magnitude of the scattering vector,
                      q, to keep in the output dataset.
    @type  q_max:     float
    @param q_max:     Minimum value of the magnitude of the scattering vector,
                      q, to keep in the output dataset.
    @rtype:           numpy array
    @return:          Two dimensional array of scattering vector, q, and
                      intensity, I, values.
    """

    scatter_data = np.loadtxt(filename)
    qrange_mask = (scatter_data[:,0] >= q_min) & (scatter_data[:,0] <= q_max)

    return scatter_data[qrange_mask]

def process_qrange_file(filename):
    """
    Load yaml file containing the ranges for all analyses.

    @type  filename:  string
    @param filename:  Path to file containing q_min and q_max values
    @rtype:           dictionary
    @return:          Dictionary containing minimum and maximum values of q
    """

    f = file(filename,'r')
    q_ranges = yaml.load(f)

    return q_ranges

def match_scatter_curves(target_data, source_data):
    """
    Match two scatter curves. Get intensity, I, values from one data set at the
    q (scaterring vector magnitude) values present in the other. Input is two q
    vs I scattering curves. Output is I values from the source data set matched
    to the q values of the target set.

    @type  target_data:  numpy array
    @param target_data:  Target scattering vector magnitude, q, and intensity,
                         I, dataset.
    @type  source_data:  numpy array
    @param source_data:  Source scattering vector magnitude, q, and intensity,
                         I, dataset.
    rtype:               numpy array
    return:              Intensity, I, data from source_data matched to the
                         scattering vector magnitudes, q, found in target_data.
    """

    # Create list of calculated I values matched to the nearest experimental q
    # Remember that the arrays in python start at 0, those in Fortran at 1
    last_source = len(source_data)
    last_target = len(target_data)

    # Initialize array to hold the calculated I values matched to
    # experimental Q values
    matched_I = np.zeros(last_target, dtype=float)

    # Use the old fortran routine to match the data sets by q value
    # matched_no is the number of datapoints which contain matched data
    matched_no = sjp_util.qrange_match(target_data[:,0], source_data[:,0],
                                       source_data[:,1], last_target,
                                       last_source, matched_I)

    matched_I.resize(matched_no)

    return matched_I

def calculate_rfactor(target_data, source_data, q_min, q_max):
    """
    Compute R factor comparing two scattering curves. Input is two
    q vs I scattering curves and the min/max q values to use to compare them.
    The target (experimental) curve is scaled to match the source (theoretical)
    one. This is because the theoretical curve is based on a calculation of
    I/Io. Output is the R factor and the scaling factor needed to match I from
    the target scattering curve to the source data.

    R factor is used by analogy with crystallography where:
    R = sum (abs(F_expt - F_calc)) / sum (abs(F_expt))

    @type  target_data:  numpy array
    @param target_data:  Target scattering vector magnitude, q, and intensity,
                         I, dataset.
    @type  source_data:  numpy array
    @param source_data:  Source scattering vector magnitude, q, and intensity,
                         I, dataset.
    @type  q_min:        float
    @param q_min:        Minimum value of the magnitude of the scattering
                         vector, q, to use in matching the curves.
    @type  q_max:        float
    @param q_max:        Minimum value of the magnitude of the scattering
                         vector, q, to use in matching the curves.
    @rtype:              float, float
    @return:             1. R factor comparing target_data and source_data.

                         2. Scaling factor needed to superimpose target_data and
                         source_data.
    """

    matched_source_I = match_scatter_curves(target_data, source_data)

    # Get average I for experimental and calculated values over matched q range
    matched_no = len(matched_source_I)
    expt_avg = np.mean( target_data[0:matched_no, 1] )
    calc_avg = np.mean( matched_source_I )

    # Initial guess of the concentration:
    # ratio of experimental and calculated average intensities
    con = expt_avg / calc_avg

    # Call fortran code to calculate the R factor
    rfactor = sjp_util.calc_rfactor( target_data[:,0], target_data[:,1],
        matched_source_I, matched_no, q_min, q_max, con, False)

    # 1/con is the scaling factor needed to multiply experimental I values
    # to compare with calculated data
    return rfactor, 1.0/con

def sas_curve_fit(x, y, calc_type):
    """
    @type  x:          numpy array
    @param x:          Function of the scattering vector magnitude, q, used for
                       the x-axis of the plot.
    @type  y:          numpy array
    @param y:          Function of the scattered intensity, I, used for the
                       y-axis of the plot.
    @type  calc_type:  string
    @param calc_type:  Specification of the type of calculation to perform on
                       the output of the linear fit:
                       - rxs1/rxs2: r = sqrt(2 * gradient), i = None
                       - rg/default: r = sqrt(3 * gradient), i = exp(intercept)
    @rtype:            dictionary
    @returns:          Dictionary containing the following keys/value pairs:
                       - fit: The output of the numpy polyfit (gradient and
                         intercept)
                       - r: The r value (Rg/Rsx?) calculated from the gradient
                       - i: The Io value calculated from the intercept
    """

    # Linear fit to the input x and y values
    fit_coeffs = np.polyfit(x, y, 1)

    result = {}

    result['fit'] = fit_coeffs

    if (calc_type == 'rxs1') or (calc_type == 'rxs2'):
        # Cross section Rxs1/Rxs2 calculated from fits of q^2 vs ln(I*q)
        # Rxs?^2 = 2 * gradient
        result['r'] = np.sqrt(2 * abs(fit_coeffs[0]))
        result['i'] = None
    else:
        # Assume a standard Guinier fit of q^2 vs ln(I):
        # Rg^2 = 3 * gradient
        # Io = exp(intercept)
        result['r'] = np.sqrt(3 * abs(fit_coeffs[0]))
        result['i'] = np.exp(fit_coeffs[1])

    return result

def graph_sas_curve(filename, x, y, title_text, x_lab, y_lab,
                  x_min, x_max, y_min, y_max, **kwargs):
    """
    Outputs graph of x and y to a pdf file. x and y are intended to be functions
    of the magnitude of scattering vector, q, and intensity, I. Values that are
    computed from the graph (outputs)
    and range of R? * qfit written on graph

    @type  x:             numpy array
    @param x:             Function of the scattering vector magnitude, q, used for the
                          x-axis of the plot.
    @type  y:             numpy array
    @param y:             Function of the scattered intensity, I, used for the y-axis
                          of the plot.
    @type  x_lab:         string
    @param x_lab:         Label for the x-axis of the plot.
    @type  y_lab:         string
    @param y_lab:         Label for the y-axis of the plot.
    @type  x_min:         float
    @param x_min:         Minimum value to plot on x-axis.
    @type  x_max:         float
    @param x_max:         Maximum value to plot on x-axis.
    @type  y_min:         float
    @param y_min:         Minimum value to plot on x-axis.
    @type  y_max:         float
    @param y_max:         Maximum value to plot on x-axis.
    @keyword fit_coeffs:  Fit coefficients for a linear fit from a numpy polyfit
                          of the input x and y values.
    @keyword outputs:     Text of values that should be added to plot (e.g. Rg).
    @keyword rq_range:    Text detailing the R*q values over the range of any
                          fit performed on the x, y values.
    @keyword mask:        Numpy style mask to select points used in fit so that
                          these can be highlighted in the plot.
    """

    fit_coeffs = kwargs.get('fitcoeffs', None)
    outputs = kwargs.get('outputs', None)
    rq_range = kwargs.get('rqrange', None)
    mask = kwargs.get('mask', None)

    # Plot the input x, y values and if provided a fit line
    # output is used to display values calculated from fit (Rg, Rxs?, R?*q over
    # fit values) on the plot

    plt.figure(figsize=(8, 6), dpi=300)

    ax = plt.subplot(111, xlabel=x_lab, ylabel=y_lab, title=title_text,
                xlim=(x_min, x_max), ylim=(y_min,y_max))

    plt.scatter(x, y, s=5)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(15)

    ax.tick_params(axis='both', which='major', labelsize=12)

    # Plot linear fit and highlight points used in its construction
    if fit_coeffs != None:

        # Plot the fit line along the whole x range shown in teh plot
        fitLine = np.poly1d(fit_coeffs)
        x_points = np.linspace(x_min, x_max, 300)
        plt.plot(x_points, fitLine(x_points))

        # Highlight points used in the fit
        if mask == None:
            plt.scatter(x, y, s=30)
        else:
            plt.scatter(x[mask], y[mask], s=30)

        plt.annotate(outputs + '\n' + rq_range, xy=(0.45, 0.85),
                     xycoords='axes fraction')

    plt.savefig(filename)

def smear_sas_curve(curve, q_delta, wavelength, spread, divergence):
    """
    Adds smearing to input theoretical  q vs I curve. Based on the formula
    published in Perkins, S. J. & Weiss, H. (1983) J. Mol. Biol. 168, 847–866.
    Note: the input I values are altered by the function in place.

    @type  curve:       numpy array
    @param curve:       Two dimensional array containing scattered vector
                        magnitudes, q, and intensities, I.
    @type  q_delta:     float
    @param q_delta:     Separation between values of q.
    @type  wavelength:  float
    @param wavelength:  Wavelength of incident beam
    @type  spread:      float
    @param spread:      Spread of beam
    @type  divergence:  float
    @param divergence:  Divergence of beam
    """

    # Calculation of the sigma for the smearing gaussian
    # Based on CHAUVIN but 0.25/D^2 removed and 8*ln(2) not 2*ln(2)
    # spread and alpha taken from Cusack JMB 1981 145, 539-541

    inv_wave_no = wavelength / (2.0 * np.pi)
    con = 4.0 * ((spread * inv_wave_no)**2 )
    bon = divergence**2
    aon = inv_wave_no * np.sqrt(8.0 * np.log(2.0))

    n_q = len(curve)
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
            f[pos_ndx] = curve[:,1][ndx - 1]
            f[neg_ndx] = curve[:,1][ndx - 1]
        else:
            f[pos_ndx] = curve[:,1][-1]
            f[neg_ndx] = curve[:,1][-1]

    sigma[mid_qq] = 2.0 * sigma[mid_qq + 1] - sigma[mid_qq + 2]

    g = np.empty(n_qq)
    asum  = np.zeros(n_qq)

    sqrt2pi = np.sqrt(2.0 * np.pi)

    # Convolution of f and g
    for ndx1 in range(0, n_qq):
        # Calculation of Gaussian, g
        for ndx2 in range(0, n_qq):
            # Note: difference from Chauvin and Cusack - ndx1 not ndx2 for aa, vv
            aa = 1.0 / (sigma[ndx1] * sqrt2pi)
            g[ndx2] = 0.0
            vv = (qq[ndx1] - qq[ndx2]) / sigma[ndx1]
            vv = vv**2
            if vv < vv_break:
                g[ndx2] = aa * np.exp(-vv / 2.0)

        for ndx2 in range(0, n_qq):
            asum[ndx1] += f[ndx2] * g[ndx2]

        asum[ndx1] *= q_delta

    for ndx in range(mid_qq + 1, mid_qq + n_q + 1):
        ndx2 = ndx - mid_qq - 1
        curve[:,1][ndx2] = asum[ndx]

def output_sas_curve(curve, out):
    """
    Prints out formated q and I columns to output.

    @type  curve:  numpy array
    @param curve:  Two dimensional array containing scattered vector
                   magnitudes, q, and intensities, I.
    @type  out:    file object
    @param out:    File object used for output of q/I curve.
    """

    for qi_pair in curve:
        out.write("{0:7.4f} {1:7.4f}\n".format(qi_pair[0], qi_pair[1]))


def get_curve_descriptors(curve, rg_min, rg_max, rxs_min, rxs_max):
    """
    Calculate the Rg and Rxs1 from linear fits to the functions of the input
    curve values (q, I). In the case of Rg the Guinier fit is to q^2 vs ln(I),
    for Rxs1 it is to q^2 vs ln(I*q).

    @type  curve:    numpy array
    @param curve:    Two dimensional array containing scattered vector
                     magnitudes, q, and intensities, I.
    @type  rg_min:   float
    @param rg_min:   Minimum q value to use in Guinier fit to compute Rg
    @type  rg_max:   float
    @param rg_max:   Maximum q value to use in Guinier fit to compute Rg
    @type  rxs_min:  float
    @param rxs_min:  Minimum q value to use in linear fit to compute Rxs1
    @type  rxs_max:  float
    @param rxs_min:  Maximum q value to use in linear fit to compute Rxs1
    @rtype:          float, float
    @return:         1. Rg value calculated from curve fit
                     2. Rxs value calculated from curve fit
    """

    # Create mask to select range of q values for Rg fitting
    rg_mask = (curve[:,0] > rg_min) & (curve[:,0] < rg_max)
    # Create mask to select range of q values for Rxs fitting
    rxs_mask = (curve[:,0] > rxs_min) & (curve[:,0] < rxs_max)

    # Fitting is performed on:
    # q^2 vs ln(I) for Rg
    # q^2 vs ln(I*q) for Rxs1
    x = curve[:,0]**2
    y_rg = np.log(curve[:,1])
    y_rxs = np.log(curve[:,1] * curve[:,0])

    rg_result = sas_curve_fit(x[rg_mask], y_rg[rg_mask], 'rg')
    rxs_result = sas_curve_fit(x[rxs_mask], y_rxs[rxs_mask], 'rxs1')

    return rg_result['r'], rxs_result['r']
