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
optimize_model_params.py [-h] -i [INPUT_PDB] [-s [INPUT_SEQ]]
                         [-f {fas,yml}] [-o [OUTPUT_FILE]] -p
                         [PARAMETER_FILE] [-m MIN_MAX MIN_MAX]
                         [-w MIN_MAX MIN_MAX]
~~~~~~~

Arguments:
----------
~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_PDB], --input_pdb [INPUT_PDB]
                        Path to the input PDB file
  -s [INPUT_SEQ], --input_seq [INPUT_SEQ]
                        Path to input sequence file if different from PDB
  -f {fas,yml}, --input_format {fas,yml}
                        Input file format (fasta or sluv yaml)
  -o [OUTPUT_FILE], --output_file [OUTPUT_FILE]
                        Path to output file
  -p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                        Path to YAML format SCT parameter file
  -m MIN_MAX MIN_MAX, --min_max MIN_MAX MIN_MAX
                        Minimum and maximum values of box side
~~~~~~~

Notes:
------

Defaults:

+ min_max = 2.0 11.0

Unless an *input_seq* is specified then the sequence of the PDB is used.

Output:
-------

The following information is output to a file or stdout:

~~~~~~~
Target dry volume: 
Target hydrated volume: 
Sphere model box side: 
Hydration cutoff: 
~~~~~~~
