% SCT: frcc documentation
% David W. Wright
% 20 Feb 2014
frcc - Calculate the Frictional Coefficient from a Sphere Model
===============================================================

Calculates the frictional coefficient from a hydrated sphere model.

Useage:
-------

Run using the command:

~~~~~~~
frcc < INPUT_MODEL
~~~~~~~

Notes:
------

INPUT_MODEL is a file containing a hydrated sphere model (output from hypro). 
This file must be a four column text file (column width of 10 characters, 2 decimal places).
Columns are X, Y, Z coordinates and sphere radius.

Output:
-------

A list of frictional cofficients components are printed to stdout. 
The only value used in the Perkins lab protocols is the following:

~~~~~~
PROPER FRICTIONAL COEFFICIENT      =
~~~~~~
