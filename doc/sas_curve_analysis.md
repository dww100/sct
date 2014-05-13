sas_curve_analysis - Script to analsye $Q$ and $I$ data from SAS experiments
============================================================================

Script to analsye $Q$ and $I$ data from SAS experiments
Can produce the following plots and data:
wide - Wide Angle: $Q$ vs $\ln(I)$
rg - Radius of gyration/Guinier: $Q^2$ vs $\ln(I)$
rxs1/rxs2 - Rxs1/Rsx2: $Q^2$ vs $\ln(I*Q)$

Useage:
-------

Run using the command:

~~~~~~~
sas_curve_analysis.py [-h]
    (-q QRANGESFILENAME | -r PLOTRANGE PLOTRANGE -f FITRANGE FITRANGE)
    (-i INPUTFILENAME | --header)
    [-s SKIPNO]
    [-o [OUTPUTDIR] [-a {wide,rg,rxs1,rxs2,all}]]
    [-y YRANGE YRANGE]
~~~~~~~

Arguments:
----------
~~~~~~~
  -i [INPUT_FILENAME], --infile [INPUT_FILENAME]
                        Path to the input file
  --header              Output a header alongside output data
  -p PARAMETER_FILE, --parameter PARAMETER_FILE
                        Path to YAML format SCT parameter file
  -r PLOTRANGE PLOTRANGE, --plotrange PLOTRANGE PLOTRANGE
                        Range of Q to plot
  -f FIT_RANGE FIT_RANGE, --fitrange FIT_RANGE FIT_RANGE
                        Range of Q to use in fitting
  -o [OUTPUT_DIR], --outdir [OUTPUT_DIR]
                        Path for output files
  -a {wide,rg,rxs1,rxs2,all}, --analysis {wide,rg,rxs1,rxs2,all}
                        Type of analysis to run
  -y Y_RANGE Y_RANGE, --yrange Y_RANGE Y_RANGE
                        Range of the y-axis of plots; ln(I) for 'wide' and
                        'rg', ln(I*Q) for 'rxs1' and 'rxs2'.
~~~~~~~

Notes:
------

If a parameter file is used it is prioritized over any parameters passed via 
command line flags.

Input is assumed to be columns of data with no header - first column
being $Q$ and the second $I$.

Output calculations will use the same length units as the input curve.

Output:
-------

Graphs:

A PDF file will be produce showing the appropriate functions of $Q$ and $I$ 
plotted against one another. 
Points used in any fitting will be shown in larger dots and with a blue line 
depicting the linear fit.
Calculated values will also be printed to the figure.

Values:

For 'rg':

Will produce four values:
[Rg] [Rg * q_min] [Rg * q_max] [I0]

For 'rxs1' or 'rxs2':

Will produce four values:
[Rxs?] [Rxs? * q_min] [Rxs? * q_max]

A header for each column can be produced using the --header option.
