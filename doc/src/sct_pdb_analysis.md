sct_pdb_analysis - Run the Full SCT Constrained Modelling Pipeline
=================================================

This script combines the functions of all of the other SCT programs to perform 
the full sphere modelling pipeline.
This involves:

Reading input experimental curves, separated into Neutron and X-ray categories.
Calculate Rg and Rxs1 from the data and write to a file
Create unhydrated sphere model
**IF** we have Neutron data:

  *  Calculate theoretical scattering curve
  *  Compare theoretical and experimental curves - R factor

**IF** we have X-ray data:

  *  Create hydrated sphere model
  *  Calculate theoretical scattering curve
  *  Compare theoretical and experimental curves - R factor or $\chi^{2}$

Write Rg, Rxs1 and comparison data to a file

Useage:
-------

Run using the command:

~~~~~~~
sct_pdb_analysis.py [-h] -i [INPUT_PATH] [-o [OUTPUT_PATH]]
                    -p [PARAMETER_FILE] [-x XRAY [XRAY ...]]
                    [-n NEUTRON [NEUTRON ...]] [-t [TITLE]]
                    [-xu {nm,a}] [-nu {nm,a}] [-ou {nm,a}]
                    [--chi2]
~~~~~~~

Arguments:
----------

~~~~~~~
-i [INPUT_PATH], --input_path [INPUT_PATH]
                        Path to the input PDB files
-o [OUTPUT_PATH], --output_path [OUTPUT_PATH]
                        Path in which to save output files
-p [PARAMETER_FILE], --parameter_file [PARAMETER_FILE]
                        Path to a file containing input parameters
-x XRAY [XRAY ...], --xray XRAY [XRAY ...]
                        Paths to files containing experimental x-ray
                        scattering curve
-n NEUTRON [NEUTRON ...], --neutron NEUTRON [NEUTRON ...]
                        Paths to files containing experimental neutron
                        scattering curve
-t [TITLE], --title [TITLE]
                        Title to use for summary output file
-xu {nm,a}, --xray_unit {nm,a}
                        Unit for Q in input x-ray data
-nu {nm,a}, --neutron_unit {nm,a}
                        Unit for Q in input neutron data
-ou {nm,a}, --output_unit {nm,a}
                        Unit for Q in output data
--chi2
                        Chose to use Chi squared rather than R factor
~~~~~~~

Notes:
------

Default TITLE = sct_output

Output:
-------

Summary file (named *TITLE_expt.sum*) of the input experimental curves with the following columns:

~~~~~~
Filename        Rg      Rxs1
~~~~~~

Summary file (named *TITLE.sum*) of the model data and comparison with experiment with the following header and columns:

~~~~~~
Input PDB path: name_edited_pdbs/
       Neutron                                                 X-ray
                                               neutron_file                                            xray_file
Model  Rg_model  Rg_curve  Rxs1_curve  Volume  Rfactor  scale  Rg_model  Rg_curve  Rxs1_curve  Volume  Rfactor  scale
~~~~~~

The Model column contains the name of the input PDB file. The headings Rg_model, Rg_curve, Rxs1_curve and Volume are repeated twice. 
The first set is for the unhydrated model used for comparison to the neutron data, the second for the hydrated model used for x-ray comparisons.
After each are the comparisons for each experimental data file (Rfactor/Chi2 and scale).
These headings are repeated for each input experimental curve.
The scale factor provided is the I(0) value that can be used to match the I(Q)/I(0) theoretical curve to the experimental data.

In the OUTPUT_PATH two sets of directories are created:
*xray/{curve,models}*
*neutron/{curve,models}*

These contain the sphere models and the theoretical curves.
