calculate_rfactor - Calculate R factor Between Two Curves
=================================================

Compares pre-calculated theorerical scattering curves to
experimental x-ray and neutron scattering curves.

Useage:
-------

Run using the command:

~~~~~~~
calculate_rfactor.py [-h] -c CALC_PATHS [CALC_PATHS ...] -e EXPT_CURVES
                     [EXPT_CURVES ...] -p [PARAMETER_FILE] -o
                     [OUT_FILE] [-u {nm,a}]
~~~~~~~

Arguments:
----------
~~~~~~~

-c CALC_PATHS [CALC_PATHS ...], --calc CALC_PATHS [CALC_PATHS ...]
                        Path to the input calculated curve or directory
                        containing mutiple curves (with scn, scx extensions)
-e EXPT_CURVES [EXPT_CURVES ...], --expt EXPT_CURVES [EXPT_CURVES ...]
                        List of input experimental curves
-p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                        Path to YAML format SCT parameter file
-o [OUT_FILE], --outfile [OUT_FILE]
                        Path to the output file
-u {nm,a}, --unit {nm,a}
                        Unit for Q in input experimental data (nanometer or Angstom)
~~~~~~~

Notes:
------

Output:
-------

A four column file, with the following data:

Path to the experimental curve, path to the theoretical curve, scale factor, R factor.

The scale factor provided is the I(0) value that can be used to match the I(Q)/I(0)
theoretical curve to the experimental data.
