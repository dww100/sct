% SCT: get_boxside documentation
% David W. Wright
% 20 Feb 2014
get_boxside - Automatic Optimization of Box Side Used in Sphere Model Creation
==============================================================================

Automatic optimization of the box side used to create sphere models from
atomistic PDBs in SCT tools.
The opimization is designed to get the minimum deviation from the volume 
computer from sluv2.py.

Useage:
-------

Run using the command:

~~~~~~~
get_boxside.py [-h] -i [INPUT_PDB] [-s [INPUT_SEQ]] [-f {fas,yml}]
               [-o [OUTPUT_FILENAME]] [-b BOX_RANGE BOX_RANGE]
               [-t [TOLERANCE]] [-c [CUTOFF]]

~~~~~~~

Arguments:

~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_PDB], --input_pdb [INPUT_PDB]
                        Path to the input PDB file
  -s [INPUT_SEQ], --input_seq [INPUT_SEQ]
                        Path to input sequence file if different from PDB
  -f {fas,yml}, --input_format {fas,yml}
                        Input file format (fasta or sluv yaml)
  -o [OUTPUT_FILENAME], --output_filename [OUTPUT_FILENAME]
                        Path to the output file
  -b BOX_RANGE BOX_RANGE, --box_range BOX_RANGE BOX_RANGE
                        Minimum and maximum box sides to try
  -t [TOLERANCE], --tolerance [TOLERANCE]
                        Tolerance for box size determination
  -c [CUTOFF], --cutoff [CUTOFF]
                        Atom number cutoff for sphere models
~~~~~~~

Notes:
------

Defaults:

+ box_range = 1.1 10.0
+ tolerance = 0.01
+ cutoff = 4

Unless an *input_seq* is specified then the sequence of the PDB is used.

Output:
-------

~~~~~~~
# Optimized box_side for INPUT_FILENAME, using cutoff CUTOFF
# Deviation from target volume was X
box_side: Y
~~~~~~~

Where X is the deviation from the theoretical volume taken from the sequence 
and Y is the optimized value for the box_side.
