% SCT: hypro documentation
% David W. Wright
% 20 Feb 2014

=================================================

Add a hydration layer to a sphere model created by *brktos* (text files consisting of four 
columns (column width of 10 characters, 2 decimal places for the X, Y and 
Z coordinates plus sphere radius, output is the same format).

Useage:
-------

Run using the command:

~~~~~~~
hypro
~~~~~~~

This will produce the following prompts, you should provide the inputs shown in 
*italics*.
Individual values to input are surrounded by square braces. 
When multiple values are requires they can either be separated by a space on 
the same line or entered on separate lines.
The meaning of the input values are explained in the **Notes** section.


Input file : *[input dry model]*

Output File  *[output hydrated model]*

How many spheres [1,7,13,19,27 Recomended] *[hydration no]*

Notes:
------

+ *input dry model* = path to the input dry sphere model file
+ *output hydrated model* = path to output hydrated sphere model file
+ *hydration no* = A number between 1 and 27 of spheres which will appear in 
the final model.
A choice of 1 puts one sphere in the position in the original model, 27 puts 
one there plus one at each of the corners and interstitial sites of a cube 
centred on the original sphere with a dimension of double the sphere radius.

A sensible route to a good hydrated model is to convert the hydrated sphere 
model to a PDB and use brktos to filter out spheres (with a low cut off) and 
then add back in the original file.
Duplicate spheres can be removed with another run through brktos.
The cut off value and number of spheres added should be altered to match the
appropriate hydrated volume calculated from sluv.
An equivalent method to this is employed in do_curve.sh, hydrate_spheres.py 
and sct_pdb_analysis.py.

Output:
-------

Coordinates: A text file consisting of four columns (column width of 10 
characters, 2 decimal places).
Columns are X, Y, Z coordinates and sphere radius.
The radius found in the input file is used for all spheres.
