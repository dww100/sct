sct_get_sequence - Convert Fasta or PDB File to YAML Residue Frequency File
===========================================================================

Convert fasta or pdb file to yaml residue frequency file.

Useage:
-------

Run using the command:

~~~~~~~
sct_get_sequence.py [-h] -i [INPUT_FILENAME] [-p]
                           [-o [OUTPUT_FILENAME]]
~~~~~~~

Arguments:

~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_FILENAME], --input_filename [INPUT_FILENAME]
                        Path to the input fasta/pdb file
  -p, --pdb             Flag to indicate input is a PDB rather than fasta file
  -o [OUTPUT_FILENAME], --output_filename [OUTPUT_FILENAME]
                        Path to the output file
~~~~~~~

Notes:
------

Output:
-------

A YAML file with lines of the following format 

[Resid]: *frequency*

where resid [Resid] is a three letter residue code and *frequency* is the 
number of occurences in the input sequence.
