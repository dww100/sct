% SCT: brktos documentation
% David W. Wright
% 20 Feb 2014 
brktos - Convert Atomistic PDB Structure into Coarse Grained Sphere Model
=========================================================================

Creates a coarse grained sphere model from an input PDB.
A grid is constructed surrounding the atoms in the structure with a specified 
box size.
Boxes within this grid containing more atoms than a prescribed cutoff are 
assigned to be spheres in the coarse grained model.

Useage:
-------

Run using the command:

~~~~~~~
brktos
~~~~~~~

This will produce the following prompts, you should provide the inputs shown in 
*italics*.
Individual values to input are surrounded by square braces. 
When multiple values are requires they can either be separated by a space on 
the same line or entered on separate lines.
The meaning of the input values are explained in the **Notes** section.

Enter name of Brookhaven format file (PDB):  *[path to PDB file]*

Enter side of box, cutoff, output type, grid type *[box length] 
[cut off] [output type] [grid type]*

Enter filename for output coords *[filename]*

Enter Number of boxes for x, y, z-axes *[X boxes] [Y boxes] [Z boxes]* 

Enter minimum x, y, z coords for cubic grid *[X min] [Y min] [Z min]*

Notes:
------

The last two prompts are dependent on the choice of *grid type*.

- *path to PDB file* = string containing location of the input PDB
- *box length* = length of the side of the boxes used to create sphere model 
(real) - this length will be 2 * sphere radius.
- *cut off* = cut off number of atoms per box to create a sphere (integer)
- *output type* = integer 0 to 4 chosing the output options (see below)
    0. Coordinates + summary data to file
    1. Coordinates only to file
    2. Summary data to screen only
    3. Coordinates to file and summary data to screen
- *grid type* = 0 or 1 to select if you want to input grid parameters (it is 
strongly suggested that you do)
    0. User inputs the grid parameters (see below)
    1. Let the warped mind of SCT do as it wished with your grid
- *filename* = path to use for output file
- *X boxes*, *Y boxes*, *Z boxes* = number of boxes to use in the X, Y and Z 
dimensions of the grid used to create sphere models (integer), typically 500 
works well.
- *X min*, *Y min*, *Z min* = lowest (generally most negative) position on the 
X, Y and Z axes of the grid (real), typically -500 works well.

Output:
-------

Depending on the choice of *output type* two different outputs are possible:

1. Coordinates: A text file consisting of four columns (column width of 10 
characters, 2 decimal places).
Columns are X, Y, Z coordinates and sphere radius.

2. Summary data:

~~~~~~
No Of Atoms
Total Of Amino Acid Residues
One Side Of The Box
Cutoff For Creating A Ball
No Of Balls
Volume Of Cubes
Volume Of Protein
Max And Min X,Y,Z Are
Increments MX, MY, MZ Are
~~~~~~

When *output type* = 0 the summary data is appended to the coordinated in output file. 
This may cause problems when used as input for other SCT programs.
