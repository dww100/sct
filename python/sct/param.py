#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to parse and check the SCT parameter YAML file
"""

# Copyright 2014 University College London

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
import sys


def output_error(message, filename):
    """
    Output error message when reading filename has failed and exit.

    @type  message:   string
    @param message:   Error message to be output.
    @type  filename:  string
    @param filename:  Path to the input file
    """

    print 'Error encountered whilst reading parameters from ' + filename + ':'
    print message
    sys.exit(1)


def read_parameter_file(filename):
    """
    Read the input YAML file

    @type  filename:  string
    @param filename:  Path to the input file
    @rtype:           dictionary, string
    @return:          1. Dictionary containing the user provided parameters for
                      SCT.

                      2. Error message (None if no error)
    """

    err = None

    try:
        param_file = file(filename)
        params = yaml.load(param_file)
    except IOError:
        err = 'Unable to read file.'
        params = None

    return params, err


def check_wide_plot(params):
    """
    Check that the parameters needed for a wide angle plot are present in
    params. Namely, params['wide']['qmin'] and params['wide']['qmax']

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    try:
        float(params['wide']['qmin'])
        float(params['wide']['qmax'])
    except:
        err = "Incomplete or invalid 'wide' parameters provided."

    return err


def check_rg_variant(variant, params):
    """
    Check that the parameters needed for a Rg type calculation are present in
    params. Namely, params[variant]['fitmin'] and params[variant]['fitmax'].
    The acceptable variants are: 'rsx1', 'rsx2' and 'rg'.

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    try:
        fitmin = float(params[variant]['fitmin'])
        fitmax = float(params[variant]['fitmax'])
    except:
        err = "Incomplete or invalid '" + variant + "' parameters provided."

    if fitmin > fitmax:
        err = variant + " fitmin must be greater than fitmax"
    elif (fitmin < params['rfac']['qmin']) or (fitmin > params['rfac']['qmax']):
        err = variant + \
            " fitmin must inside the rfac calculation range (between input rfac qmin and qmax)"
    elif (fitmax < params['rfac']['qmin']) or (fitmax > params['rfac']['qmax']):
        err = variant + \
            " fitmin must inside the rfac calculation range (between input rfac qmin and qmax)"

    return err


def check_rg(params):
    """
    Check that the parameters needed for a Rg calculation are present in
    params. Namely, params['rg']['fitmin'] and params['rg']['fitmax'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = check_rg_variant('rg', params)

    return err


def check_rxs1(params):
    """
    Check that the parameters needed for a Rxs1 calculation are present in
    params. Namely, params['rxs1']['fitmin'] and params['rxs1']['fitmax'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = check_rg_variant('rxs1', params)

    return err


def check_rxs2(params):
    """
    Check that the parameters needed for a Rxs1 calculation are present in
    params. Namely, params['rxs1']['fitmin'] and params['rxs1']['fitmax'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = check_rg_variant('rxs2', params)

    return err


def check_rg_plot_variant(variant, params):
    """
    Check that the parameters needed for a Rg type plot are present in
    params. Namely, params[variant]['qmin'] and params[variant]['qmax'],
    params[variant]['fitmin'] and params[variant]['fitmax'].
    The acceptable variants are: 'rsx1', 'rsx2' and 'rg'.

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    try:
        float(params[variant]['qmin'])
        float(params[variant]['qmax'])
        err = check_rg_variant(variant, params)
    except:
        err = "Incomplete or invalid '" + variant + "' parameters provided."

    return err


def check_rg_plot(params):
    """
    Check that the parameters needed for an Rg plot are present in
    params. Namely, params['rg']['qmin'] and params['rg']['qmax'],
    params['rg']['fitmin'] and params['rg']['fitmax'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = check_rg_plot_variant('rg', params)

    return err


def check_rxs1_plot(params):
    """
    Check that the parameters needed for an Rxs1 plot are present in
    params. Namely, params['rxs1']['qmin'] and params['rxs1']['qmax'],
    params['rxs1']['fitmin'] and params['rxs1']['fitmax'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = check_rg_plot_variant('rxs1', params)

    return err


def check_rxs2_plot(params):
    """
    Check that the parameters needed for an Rxs2 plot are present in
    params. Namely, params['rxs2']['qmin'] and params['rxs2']['qmax'],
    params['rxs2']['fitmin'] and params['rxs2']['fitmax'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = check_rg_plot_variant('rxs2', params)

    return err


def check_rfac(params):
    """
    Check that the parameters needed for a R factor calculation are present in
    params. Namely, params['rfac']['qmin'] and params['rfac']['qmax'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    try:
        float(params['rfac']['qmin'])
        float(params['rfac']['qmax'])
    except:
        err = "Incomplete or invalid 'rfac' parameters provided."

    return err


def check_curve(params):
    """
    Check that the parameters needed for a theoretical curve calculation are
    present in params. Namely:
    params['curve']['qmax'], params['curve']['npoints'],
    params['curve']['radbins'], params['curve']['smear'],
    params['curve']['wavelength'], params['curve']['spread']
    and params['curve']['divergence'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    try:
        float(params['curve']['qmax'])
        int(params['curve']['npoints'])
        int(params['curve']['radbins'])
        bool(params['curve']['smear'])
        float(params['curve']['wavelength'])
        float(params['curve']['spread'])
        float(params['curve']['divergence'])
    except:
        err = "Incomplete or invalid 'curve' parameters provided."

    return err


def check_sphere(params):
    """
    Check that the parameters needed for a sphere model creation are
    present in params. Namely: params['sphere']['cutoff'] and
    params['sphere']['boxside'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    try:
        int(params['sphere']['cutoff'])
        float(params['sphere']['boxside'])

    except:
        err = "Incomplete or invalid 'sphere' parameters provided."

    return err


def check_hydrate(params):
    """
    Check that the parameters needed to hydrate a sphere model are
    present in params. Namely: params['hydrate']['cutoff'] and
    params[hydrate]['boxside'].

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    try:
        int(params['hydrate']['positions'])
        float(params['hydrate']['cutoff'])

    except:
        err = "Incomplete or invalid 'hydrate' parameters provided."

    return err


def check_parameters(params, needed):
    """
    Check that all the parameters needed to perform the calculations specified
    in the input list (needed) are in the parameter dictionary (params).
    Possible calculations are:
    'rg', 'rg_plot' - Radius of gyration (calculation alone and with plot also)
    'rxs1', 'rxs1_plot' - Rxs1 (calculation alone and with plot also)
    'rxs2', 'rxs2_plot' - Rxs2 (calculation alone and with plot also)
    'sphere' - creation of a sphere model
    'hydrate' - hydration of sphere model for x-ray comparisons
    'curve' - calculation of a theoretical curve
    'rfac' - calculation of the R factor

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @type  needed:  list
    @param needed:  List of strings describing the calculations for which
                    parameters will be required.
    @rtype:         string
    @return:        Error message (None if no error)
    """

    err = None

    if params is None:
        err = 'No parameters could be read from the input YAML file.'
    else:
        for param_set in needed:
            fn = 'check_' + param_set
            err = globals().get(fn)(params)

    return err


def parse_parameter_file(parameter_file, needed):
    """
    Load the parameter file abd check that all the parameters needed to perform
    the calculations specified in the input list (needed) are in the parameter
    dictionary (params).
    Possible calculations are:
    'rg', 'rg_plot' - Radius of gyration (calculation alone and with plot also)
    'rxs1', 'rxs1_plot' - Rxs1 (calculation alone and with plot also)
    'rxs2', 'rxs2_plot' - Rxs2 (calculation alone and with plot also)
    'sphere' - creation of a sphere model
    'hydrate' - hydration of sphere model for x-ray comparisons
    'curve' - calculation of a theoretical curve
    'rfac' - calculation of the R factor

    @type  params:  dictionary
    @param params:  Dictionary containing the user provided parameters for
                    SCT.
    @type  needed:  list
    @param needed:  List of strings describing the calculations for which
                    parameters will be required.
    @rtype:         string
    @return:        Error message (None if no error)
    """
    param, err = read_parameter_file(parameter_file)

    if err is not None:
        output_error(err, parameter_file)

    err = check_parameters(param, needed)

    if err is not None:
        output_error(err, parameter_file)

    if 'sphere' in needed:
        param['sphere']['boxside3'] = param['sphere']['boxside'] ** 3
        # Model spheres just fit into grid box (no overlap) therefore:
        param['sphere']['radius'] = param['sphere']['boxside'] / 2.0

    if 'curve' in needed:
        param['curve']['q_delta'] = param['curve'][
            'qmax'] / param['curve']['npoints']

    return param
