% SCT: calculate_curve documentation
% David W. Wright
% 20 Feb 2014
calculate_curve - Calculate Scattering Curve from a Sphere Model
================================================================

Calculate scattering curve from sphere model

Useage:
-------

Run using the command:

~~~~~~~
calculate_curve.py [-h] -i [INPUT_FILENAME] [-o [OUTPUT_FILENAME]]
                          [-q [Q_MAX]] [-r [RADIUS]] [-b [N_BINS]]
                          [-p [N_POINTS]] [-s] [-e [SPREAD]] [-w [WAVELENGTH]]
                          [-d [DIVERGENCE]]
~~~~~~~

Arguments:

~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_FILENAME], --input_filename [INPUT_FILENAME]
                        Path to the input PDB file
  -o [OUTPUT_FILENAME], --output_filename [OUTPUT_FILENAME]
                        Path to the output file
  -q [Q_MAX], --q_max [Q_MAX]
                        Maximum q value in output curve
  -r [RADIUS], --radius [RADIUS]
                        Sphere radius
  -b [N_BINS], --n_bins [N_BINS]
                        No. bins to use in histogram of sphere separation
  -p [N_POINTS], --n_points [N_POINTS]
                        No. points in output curve
  -s, --smear           Apply smearing to curve
  -e [SPREAD], --spread [SPREAD]
                        Wavelength spread used to calculate smearing
  -w [WAVELENGTH], --wavelength [WAVELENGTH]
                        Wavelength used to calculate smearing
  -d [DIVERGENCE], --divergence [DIVERGENCE]
                        Beam divergence used to calculate smearing
~~~~~~~

Notes:
------

Defaults:

+ q_max = 0.16
+ radius = 3.77
+ n_bins = 400
+ n_points = 100
+ spread = 0.1
+ wavelength = 6.0
+ divergence = 0.016

Output:
-------

Two column text file - with $Q$ and $I/I_{0}$ columns (floating point numbers 
as strings formatted using 7.4f in Python).

$Q$ units are $\AA^{-1}$.
