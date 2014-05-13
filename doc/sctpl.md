sctpl - Calculate Radii of Gyration From Scattering Curves
=================================================

Useage:
-------

Run using the command:

~~~~~~~
sctpl
~~~~~~~

This will produce the following prompts, you should provide the inputs shown in 
*italics*.
Individual values to input are surrounded by square braces. 
When multiple values are requires they can either be separated by a space on 
the same line or entered on separate lines.
The meaning of the input values are explained in the **Notes** section.

Enter data type (model [m] or exptl [x]): *[curve type]*

Enter title of plot: *[title]*

Number of curves? [1] *[no. curves]*

Enter name of file containing data: *[input file]*

Smearing (y/n)? *[smearing]*

Type of plot:
 0 Ln(IQQ) / Q*Q: Thickness RG
 1 Ln(IQ) / Q*Q: Cross-sectional RG
 2 Ln(I) / Q*Q: Overall RG; Guinier plot
 3 Ln(I) / Q: Wide Angle
 4 I / Q
 5 I*Q*Q / Q: Volume plot

Choose one of the above [3]: *[plot type]*
Which Q range to plot (Q1, Q2; Q1 < Q2) *[Q min]* *[Q max]*

Which Q range to fit a line to (Q3, Q4; Q3 < Q4)? *[fit min]* *[fit max]*

Choose a plot style:
Join the dots, no symbols [0]: Symbols only [1]: Dots and symbols [2]. *[plot style]*

Take your pick:
1) Get another;
2) Add another;
3) Change axes' ranges;
4) Start again;
5) Manual R-factor;
6) Auto R-factor;
7) Change a given plot's attributes;
8) Save selected plots as ASCII;
9) Screen dump;
10) Redisplay a given plot;
11) Remove a plot;
12) Quit:
*[menu choice]*

Notes:
------

The original *sctpl* performed plotting, this is not a current feature but the 
options remain for compatibility with old scripts.
Plotting and fitting can now be performed using *sas_curve_analysis.py*.

*curve type* = choose of model [m] or exptl [x] curve

*title* = title for plot

*no. curves* = number of curves to analyse

*input file* = path to input scattering curve (with $I$ and $Q$ columns).

*smearing* = should the curve be smeared

*plot type* = choice from 'Type of plot' menu shown above.
This defines the radius of gyration type to be out put:

0 radius of gyration ($R_{g}$)
1 cross sectional $R_{g}$ (usually called $R_{xs1}$)
2 thickness $R_{g}$ (usually called $R_{xs2}$)

*Q min* = minimum Q value for plot

*Q max* = maximum Q value for plot

*fit min* = minimum Q value for fit

*fit max* = maximum Q value for fit

*plot style* = type of points used in plot

*menu choice* = choice from 'Take your pick' menu

Output:
-------

*sctpl* produces a single output file in the current directory called 
'Sctpl6_Summary'.

This has the following format:

~~~~~~
*input file*,*calc type* ,*title*, *fit Q range*, *fit points*, *result*, *result error*,*I0*, *I0 error*
~~~~~~

*calc type* =
+ RG = radius of gyration ($R_{g}$)
+ RXS = cross sectional $R_{g}$ (usually called $R_{xs1}$)
+ RTH = thickness $R_{g}$ (usually called $R_{xs2}$)

*fit Q range* = *Q min* -> *Q max*

*fit points* = points on teh input curve that correspond to *fit Q range*

*result* = result of the radius of gyration calculation in the units of the inverse of the input Q values

*result err* = uncertainty in the *result* from curve fitting

*I0* = extrapolated $I_{0}$ value

*I0 error* = uncertainty in the $I_{0}$ value from curve fitting

