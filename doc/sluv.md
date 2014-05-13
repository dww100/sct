sluv - Calculate Volume and Other Properties from Sequence
==========================================================

Calculate molecular weight the partial specific volume and the molecular 
extinction coefficient of proteins (and glycoproteins) from their residue 
composition.

Key references and sources of data are:

+ S. J. Perkins (1986). Protein volumes and hydration effects: The 
calculations of partial specific volumes, neutron scattering matchpoints and 
280-nm absorption coefficients for proteins and glycoproteins from amino acid 
sequences. Eur. J. Biochem. 157, 169-180

+ S. J. Perkins (2001). X-ray and neutron scattering analyses of hydration 
shells: a molecular interpretation based on sequence predictions and model 
fits. Biophysical Chemistry 93, 129-139


Useage:
-------

Run using the command:

~~~~~~~
sluv < INPUT_FILE
~~~~~~~

INPUT_FILE is a specifically formated file containing residue frequencies (see 
**Notes**).

Notes:
------

The input file for sluv contains residue frequencies in a particular order and 
should be generated from the [input template](sluv_input_template), all lines 
must be retained but edited as follows:

+ "Template 0" can be replaced with any title
+ Delete each 3-letter code amino acid/glycan code and replace it by the 
frequency of that residue in the sequence

Output:
-------

The program spews a lot of data out to stdout (it may well be adviseable to 
pipe this to a file for future reference). 
The most important parameters are found in the secton titled:

~~~~~~
************* TOTAL GLYCOPROTEIN ********************************************
~~~~~~

Important values provided include:

~~~~~~
MOLECULAR WEIGHT
EXTINCTION COEFF (280 NM)
~~~~~~

and the values in the CON85 columns (the consensus values from Perkins 1986) for

~~~~~~
VOLUME
SPECIFIC VOLUMES
~~~~~~



