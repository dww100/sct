% SCT: optimize_model_parameters documentation
% David W. Wright
% 20 Feb 2014
optimize_model_parameters - Optimize paramateres used for Sphere Modelling
==========================================================================

Take an input PDB and optimize the box size and cutoff used to create sphere
models.
The optimization is run to recreate the volume predicted from sluv2.
The box size determines the size of the grid used to create the sphere models
from a PDB.
The cutoff is used in the creation of the hydration layer.

These parameters can then be put into the input file (yaml) for the
analyse_pdb_models.py script to run on a series of models.


Useage:
-------

Run using the command:

~~~~~~~
optimize_model_params.py [-h] -i [INPUT_PDB] [-o [OUTPUT_FILE]] -p
                                [PARAMETER_FILE] [-m MIN_MAX MIN_MAX]
~~~~~~~

Arguments:

~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_PDB], --input_pdb [INPUT_PDB]
                        Path to the input PDB file
  -o [OUTPUT_FILE], --output_file [OUTPUT_FILE]
                        Path to output file
  -p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                        Path to a file containing input parameters
  -m MIN_MAX MIN_MAX, --min_max MIN_MAX MIN_MAX
                        Minimum and maximum values of box side
~~~~~~~

Notes:
------

Defaults:

+ min_max = 2.0 11.0

Output:
-------

The following information is output to a file or stdout:

~~~~~~~
Target dry volume: 
Target hydrated volume: 
Sphere model box side: 
Hydration cutoff: 
~~~~~~~
