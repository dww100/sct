% SCT: sluv2 documentation
% David W. Wright
% 20 Feb 2014
sluv2 - Calculate Volume and Other Properties from Sequence
==========================================================

Calculate the scattering length per unit volume for a protein or glycoprotein
from its amino acid and carbohydrate composition.

This is a reimplementation of the sluv tool originally created by
Stephen J. Perkins in 1979. Documentation of the methods used can be found in:

[SJP1] Perkins, S. J. (1986). Protein volumes and hydration effects: The
calculations of partial specific volumes, neutron scattering matchpoints
and 280-nm absorption coefficients for proteins and glycoproteins from
amino acid sequences. Eur. J. Biochem. 157, 169-180

[SJP2] Perkins, S. J. (2001). X-ray and neutron scattering analyses of
hydration shells: a molecular interpretation based on sequence predictions
and model fits. Biophysical Chemistry 93, 129–139

Useage:
-------

Run using the command:

~~~~~~~
sluv2.py [-h] -i [INPUT_FILE] [-t {fas,pdb,yml}] [-o [OUTPUT_FILE]]
                [-m {classic,model,auc,project}]
~~~~~~~

Arguments:

~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_FILE], --input_file [INPUT_FILE]
                        Path to the input composition file
  -t {fas,pdb,yml}, --input_type {fas,pdb,yml}
                        Input file format (pdb, fasta or sluv yaml)
  -o [OUTPUT_FILE], --output_file [OUTPUT_FILE]
                        Path to the output file
  -m {classic,model,auc,project}, --mode {classic,model,auc,project}
                        Type of analysis to run
~~~~~~~

Notes:
------

Output:
-------

'auc' mode selected:

~~~~~~
Molecular Weight:
Absorption Coefficient x 1.03:
Specific Volume (Perkins 1986 - Consensus):
~~~~~~

'model' mode selected:

~~~~~~
Volume (Chothia 1975 - Crystal Structures):

Volume (Perkins 1986 - Consensus):
~~~~~~

'project' mode selected:

The output from both 'auc' and 'model' is produced.

'classic mode selected':

The program spews out a lot of data (it is adviseable to set an output file so 
this is stored for future reference). 
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
