% SCT: hydrate_spheres documentation
% David W. Wright
% 20 Feb 2014
hydrate_spheres - Add Hydration Layer to a Sphere Model
=================================================

Add hydration shell of spheres around an existing sphere model read from a SCT
format sphere file.

Useage:
-------

Run using the command:

~~~~~~~
hydrate_spheres.py [-h] -i [INPUT_FILENAME] -o [OUTPUT_FILENAME]
                          [-n [HYDRATION_NO]] [-c [CUTOFF]]
~~~~~~~

Arguments:

~~~~~~~
  -h, --help            show this help message and exit
  -i [INPUT_FILENAME], --input_filename [INPUT_FILENAME]
                        Path to the input file
  -o [OUTPUT_FILENAME], --output_filename [OUTPUT_FILENAME]
                        Path to the output file
  -n [HYDRATION_NO], --hydration_no [HYDRATION_NO]
                        No. spheres to add as a hydration shell (1-26)
  -c [CUTOFF], --cutoff [CUTOFF]
                        Cutoff used to reduce the number of water spheres
~~~~~~~

Notes:
------

Defaults:

+ hydration_no = 26
+ cutoff = 5

Output:
-------

Hydrated Model:
Output as four column - x, y and z coordinates and sphere radius (floating 
point numbers string formatted using 10.2f from Python).
