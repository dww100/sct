brksph2pdb - Convert a Sphere Model to PDB Format
=================================================

Convert sphere models created by *brktos* (text files consisting of four 
columns (column width of 10 characters, 2 decimal places for the X, Y and 
Z coordinates plus sphere radius) to a PDB format.
Each sphere is assigned as a SERine C1 atom.

Useage:
-------

Run using the command:

~~~~~~~
brksph2pdb
~~~~~~~

You will then me prompted to provide paths to the input (sphere model) and 
output (PDB) files.

Notes:
------
The name of this utility come from a historical habit to name sphere model 
files *.brk. 
This is however confusing as .brk is an archaic file ending for PDB files.
I can only apologize for the current name created as a disambiguation of 
the old name brk2pdb when it was annoying me more than it does now.

Output:
-------

A PDB file containing SERine C1 atoms located at the coordinates of the 
spheres passed as input (the temperature factor column is used for the 
sphere radius).

~~~~~~
ATOM      1  C1  SER  1474     110.000  20.000  35.000  2.50
~~~~~~
