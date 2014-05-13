aps - Calculation of Radius of Gyration from Sphere Models
==========================================================

Calculates the radius of gyration directly from the coordinates of a sphere model.
Also produces some other summary data.

Useage:
-------

Run using the command:

~~~~~~~
aps < INPUT_MODEL
~~~~~~~

Notes:
------
INPUT_MODEL is a file containing a sphere model. 
This file must be a four column text file (column width of 10 characters, 2 decimal places).
Columns are X, Y, Z coordinates and sphere radius.
With some options brktos will produce sphere model files with additional data at the bottom, 
this must be removed before they are used as input to *aps*.

Output:
-------

Text giving summary of input file (no. spheres, model volume, etc). and the radius of gyration in Angstrom. 
Full output contains the following lines:

~~~~~~
 COORDINATE READ-IN DONE : TOTAL IS  

 SEMIAXES OF RECT. PRISM    =       
 RADIUS OF GYRATION OF SINGLE RECT PRISM =  
 VOLUME OF SINGLE BOX =             
 TOTAL VOLUME OF ALL BOXES   =     

 NO OF COORDINATES =     
 SUM OF ALL X,Y,Z COORDS =     
 MEAN X,Y,Z COORDS ARE       
 RADIUS OF UNIT WEIGHT=     
 SUM OF RADII SQ      =   

 RADIUS OF GYRATION   =  
~~~~~~

