pdb2sphere - Convert Atomistic PDB Structure into Coarse Grained Sphere Model
=============================================================================

Convert atomic pdb to sphere model

Useage:
-------

Run using the command:

~~~~~~~
pdb2sphere.py [-h] -i [INPUT_FILENAME]
              [-m {info,archive,model,project}] [-o [OUTPUT_FILENAME]]
              [-c [CUTOFF]] [-b [BOX_SIDE]] [-p [PARAMETER_FILE]]

~~~~~~~

Arguments:
----------
~~~~~~~
  -i [INPUT_FILENAME], --input_filename [INPUT_FILENAME]
                        Path to the input PDB file
  -m {info,archive,model,project}, --mode {info,archive,model,project}
                        Type of output to be produced
  -o [OUTPUT_FILENAME], --output_filename [OUTPUT_FILENAME]
                        Path to the output file
  -c [CUTOFF], --cutoff [CUTOFF]
                        Cutoff number of atoms to create sphere in grid box
  -b [BOX_SIDE], --box_side [BOX_SIDE]
                        Side length of grid boxes
  -p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                        Path to YAML format SCT parameter file
~~~~~~~

Notes:
------

If a parameter file is used it is prioritized over any parameters passed via 
command line flags.

An input file is required.

Units of boxside are Angstrom.

Options for the type of output to be produced:

- info = summary information on sphere model to stdout only
- archive = summary information and sphere model to file
- model = model to file (no summary information produced)
- project = model to file and summary information to stdout

Defaults:

- mode = project
- cutoff = 4
- boxside = 5.0

Output:
-------

Model:
Output as four column - x, y and z coordinates and sphere radius (floating 
point numbers string formatted using 10.2f from Python).

Summary information:

~~~~~~~
pdb2sphere: version 0.5 - 05 November 2013
No. Of Atoms: 
Total Of Amino Acid Residues 
One Side Of The Box: 
Cutoff For Creating A Ball: 
No Of Balls: 
Volume Of Cubes: 
Volume of Protein:
Max And Min X: 
Max And Min Y: 
Max And Min Z: 
No. grid points in X, Y, Z: 
~~~~~~~

