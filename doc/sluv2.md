% SCT: sluv2 documentation
% David W. Wright
% 20 Feb 2014
sluv2 - Calculate Volume and Other Properties from Sequence
==========================================================

Calculate the scattering length per unit volume for a protein or glycoprotein
from its amino acid and carbohydrate composition.

This is a reimplementation of the sluv tool originally created by
Stephen J. Perkins in 1979. Documentation of the methods used can be found in:

+ Perkins, S. J. (1986). Protein volumes and hydration effects: The
calculations of partial specific volumes, neutron scattering matchpoints
and 280-nm absorption coefficients for proteins and glycoproteins from
amino acid sequences. Eur. J. Biochem. 157, 169-180

+ Perkins, S. J. (2001). X-ray and neutron scattering analyses of
hydration shells: a molecular interpretation based on sequence predictions
and model fits. Biophysical Chemistry 93, 129-139

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

Whilst the PDB and FASTA file formats are external standards the sluv YAML 
format contains simple residue name and frequency pairs using the YAML format.
An [input template](sluv_in.yml) is provided with number of each residue set
to 0, replace this value for your system of interest or use 
_sct_get_sequence.py_ to generate the YAML file.

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

Volume (Perkins 1986 - Amino Acid Crystals):
Estimated Hydrated Volume:

~~~~~~

The estimated hydrated volume comes from assuming that 0.3 g of water is added 
for each g of protein.

'project' mode selected:

The output from both 'auc' and 'model' is produced.

'classic' mode selected:

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

and the values in the PER85 columns (the values taken from crystal structures in Perkins 1986) for

~~~~~~
VOLUME
SPECIFIC VOLUMES
~~~~~~
