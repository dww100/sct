% SCT: sphere2pdb documentation
% David W. Wright
% 20 Feb 2014
sphere2pdb - Convert a Sphere Model to PDB Format
=================================================

Convert sphere models created by SCT (text files consisting of four 
columns (for the X, Y and Z coordinates plus sphere radius) to a PDB format.
Each sphere is assigned as a SERine C1 atom.

Useage:
-------

Run using the command:

~~~~~~~
sphere2pdb.py [-h] -i [INPUT_FILENAME] -o [OUTPUT_FILENAME]
~~~~~~~

Arguments:

~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_FILENAME], --input_filename [INPUT_FILENAME]
                        Path to the input file
  -o [OUTPUT_FILENAME], --output_filename [OUTPUT_FILENAME]
                        Path to the output file
~~~~~~~

Notes:
------

Output:
-------

A PDB file containing SERine C1 atoms located at the coordinates of the 
spheres passed as input (the temperature factor column is used for the 
sphere radius).

