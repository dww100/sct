sct - Calculate Scattering Curve from a Sphere Model
====================================================

Create a theoretical scattering curve (I/I(0) vs $Q$) from a sphere model using the 
Debye equation.
Units of $Q$ are Angstrom$^{-1}$.
For theoretical neutron curves smearing can be added as described in 
S. J. Perkins & H. Weiss, (1983) J. Mol. Biol. 168, 847-866

Useage:
-------

Run using the command:

~~~~~~~
sct
~~~~~~~

This will produce the following prompts, you should provide the inputs shown in 
*italics*.
Individual values to input are surrounded by square braces. 
When multiple values are requires they can either be separated by a space on 
the same line or entered on separate lines.
The meaning of the input values are explained in the **Notes** section.

Enter name of file containing model coords *[sphere model]*

Enter name of file for the calc curve *[curve file]*

Enter name of file containing more info *[output file]*

TITLE ? *[title]*

DEFAULT PARAMETERS? Y=1 N=0 *[use defaults]*

**If you answer '0'** to the above question you will recieve the following 
prompts too:

+ PLOTS?     Y=1 N=0 *[plots]*
+ SMEARING?     Y=1 N=0 *[smearing]*
+ **If you answer '0'** to the smearing question you will recieve the following 
prompts too:
    - WAVELENGTH SPREAD? 0.13?
    - WAVELENGTH? 5.97?
    - BEAM DIVERGENCE? 0.016?  NOT 0 *[spread]* *[wavelength]* *[divergence]*

+ MAXIMUM Q VALUE? *[max Q]*
+ RADIUS TO BE LOW AND CHECKED RADIUS OF ONE SPHERE? *[rr]*


Notes:
------

- *sphere model* = input sphere model coordinates in a text file consisting of 
four columns (column width of 10 characters, 2 decimal places), as generated 
by *brktos*.
Columns are X, Y, Z coordinates and sphere radius

- *curve file* = path to output $I/I_o$ vs $Q$ file

- *output file* = path to file which will summary output concerning the curve 
calculation

- *plots* = old parameter used to invoke GUI, no longer available kept for 
historical compatability

- *smearing* = choice of whether to include smearing (for comparison with 
neutron experimental results)

- *spread* = spread of the neutron beam for smearing calculation 
($\Delta \lambda / \lambda$)

- *wavelength* = wavelength of neutron beam for smearing calculation 
($\lambda$)

- *divergence* = divergences of neutron beam for smearing calculation

- *max Q* = maximum value of Q to claculate the curve up to

- *rr* = radius of sphere used in form factor calculation for the Debye equation

Output:
-------

*curve file* contains two columns $I/I_{0}$ and $Q$.
$Q$ is in units of $\AA^{-1}$.
