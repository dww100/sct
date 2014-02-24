% SCT: rfacXN documentation
% David W. Wright
% 20 Feb 2014
rfacXN - Calculate R Factor Comparing Two Curves
=================================================

Compare two scattering curves (usually one theoretical and one experimental) 
and compute the R factor. 
R factor is used by analogy with crystallography where:
$R = sum (abs(F_expt - F_calc)) / sum (abs(F_expt))$

Useage:
-------

Run using the command:

~~~~~~~
rfacXN
~~~~~~~

This will produce the following prompts, you should provide the inputs shown in 
*italics*.
Individual values to input are surrounded by square braces. 
When multiple values are requires they can either be separated by a space on 
the same line or entered on separate lines.
The meaning of the input values are explained in the **Notes** section.

Provide the name of the file containing the calculated data: *[model file]*

Provide the name of the file containing the experimental data: *[expt file]*

Provide the output filename: *[output file]*

What is the maximum Q value to be used in the fit? *[Max Q]*

Notes:
------

+ *model file* = path to file containing Q and I values computed from a 
sphere model
+ *expt file* = path to file containing experimental Q and I values
+ *output file* = path to file in which output will be written 
+ *Max Q* = maximum Q value to include in the curve comparison

It is assumed that the input experimental curve has been truncated at low Q if 
necessary.

Output:
-------

During the run the following summary data is calculated from the input data:

~~~~~~
 Mean Observed=
 Mean Calculated=
 Initial Guess=
 Max Xptl Q=
 Max Calc Q=
~~~~~~

The output file has the following format:

~~~~~~
*[model file]* *[Max Q]* *[Scaling factor]* *[R factor]*
~~~~~~

Where *model file* and *Max Q* are the values input by the user and:

+ *R factor* = the R factor as a percentage
+ *Scaling factor* = the multiplicative factor used to scale the model curve 
to match the experimental Io
