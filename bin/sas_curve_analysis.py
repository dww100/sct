#!/usr/bin/env python

"""Script to analsye Q and I data from SAS experiments
Can produce the following plots and data:
wide - Wide Angle:  Q vs ln(I)
rg - Radius of gyration/Guinier:     Q^2 vs ln(I)
rxs1/rxs2 - Rxs1/Rsx2:   Q^2 vs ln(I*Q)
At the minute we assume input has no header and the
first column is Q and second is I

David W Wright - 09/08/2013"""

import os
import argparse
from pylab import *
import yaml

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    use_message="""%(prog)s [-h]
    (-q QRANGESFILENAME | -r PLOTRANGE PLOTRANGE -f FITRANGE FITRANGE)
    (-i INPUTFILENAME | --header)
    [-s SKIPNO]
    [-o [OUTPUTDIR] [-a {wide,rg,rxs1,rxs2,all}]]
    [-y YRANGE YRANGE]
    """
    #Command line option parsing
    parser = argparse.ArgumentParser(
        description= """Basic processing of SAS data.
        Produces Rg, Io, Rxs1 and Rxs2 estimates and associated plots.
        """, usage=use_message)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i','--infile', nargs='?', type=str,
                       dest='input_filename', help = 'Path to the input file')
    group.add_argument('--header', action='store_true',
                       help = 'Output a header alongside output data')

    parser.add_argument('-s','--skip', nargs = '?', type = int, default = 0,
                       help = 'Number of lines in input file to skip')

    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument('-q','--q-ranges', type = str, dest='q_ranges_filename',
                       help = 'Path to file containing Q plot and fit ranges')
    group2.add_argument('-r','--plotrange', nargs = 2, type = float,
                       help = 'Range of Q to plot')

    parser.add_argument('-f','--fitrange', nargs = 2, type = float,
                       dest='fit_range', help = 'Range of Q to use in fitting')
    parser.add_argument('-o','--outdir', nargs='?', type=str, default = '.',
                       dest='output_dir', help = 'Path for output files')
    parser.add_argument('-a','--analysis',
                       choices = ['wide', 'rg', 'rxs1', 'rxs2','all'],
                       default = 'wide', dest='anal_type',
                       help = 'Type of analysis to run')
    parser.add_argument('-y','--yrange', nargs = 2, type = float,
                       default = None, dest='y_range',
                       help = "Range of the y-axis of plots; ln(I) for 'wide' and 'rg', ln(I*Q) for 'rxs1' and 'rxs2'.")

    args = parser.parse_args()

    if (args.anal_type == 'all') and (not args.q_ranges_filename):
        print ("If you require all analyses to be run a file containig the relevant Q ranges must be provided.\n")
        sys.exit(1)
    elif (not args.q_ranges_filename) and ((not args.fit_range) or (args.anal_type == 'wide')):
        print ("For all analyses except wide angle plotting ('wide') a fit range is needed in addition to the plot range.\n")
        sys.exit(1)

    return args

def process_range_file(filename):
    """Load file containing the ranges for all analyses."""

    f = file(filename,'r')
    q_ranges = yaml.load(f)

    return q_ranges

def create_header(analysis_type):
    """Creates a string containing column headings for the data produced by analyses
    selected by analType"""

    header = ''

    if analysis_type in ('rg','all'):
        header += 'Rg\tRgQmin\tRgQmax\tIo\t'

    if analysis_type in ('rxs1','all'):
        header += 'Rxs1\tRxs1Qmin\tRxs1Qmax\t'

    if analysis_type in ('rxs2','all'):
        header += 'Rxs2\tRxs2Qmin\tRxs2Qmax\t'

    return header

def create_figure(filename, x, y, title_text, x_lab, y_lab,
                  x_min, x_max, y_min, y_max, **kwargs):
    """ Outputs graph of x and y to a pdf file
    Values computed from graph (outputs) and range of R? * Qfit written on graph
    """

    fit_coeffs = kwargs.get('fitcoeffs', None)
    outputs = kwargs.get('outputs', None)
    rq_range = kwargs.get('rqrange', None)
    mask = kwargs.get('mask', None)

    # Plot the input x, y values and if provided a fit line
    # Also output data calculated from fit and R*Q over fit valies

    figure(figsize=(8, 6), dpi=300)

    ax = subplot(111, xlabel=x_lab, ylabel=y_lab, title=title_text,
                xlim=(x_min, x_max), ylim=(y_min,y_max))

    scatter(x, y, s=5)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(15)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Plot linear fit and highlight points used in its construction
    if fit_coeffs != None:
        fitLine = poly1d(fit_coeffs)
        x_points = linspace(x_min, x_max, 300)
        plot(x_points, fitLine(x_points))
        if mask == None:
            scatter(x, y,s=30)
        else:
            scatter(x[mask], y[mask],s=30)
        plt.annotate(outputs + '\n' + rq_range, xy=(0.45, 0.85),
                    xycoords='axes fraction')

    savefig(filename)

def create_output_name(prefix, analysis, q_min, q_max, fit_min, fit_max):

    q_range = "-".join([str(q_min),str(q_max)])
    fit_range = "-".join([str(fit_min),str(fit_max)])

    filename = "{0:s}_q{1:s}_fit{2:s}_{3:s}.pdf".format(
                prefix, q_range, fit_range, analysis)

    return filename

def guinier_fit(x, y, calc_type):

    fit_coeffs = polyfit(x, y, 1)

    result = {}

    result['fit'] = fit_coeffs

    if (calc_type == 'rxs1') or (calc_type == 'rxs2'):
        result['r'] = sqrt(2 * abs(fit_coeffs[0]))
        result['i'] = None
    else:
        result['r'] = sqrt(3 * abs(fit_coeffs[0]))
        result['i'] = exp(fit_coeffs[1])

    return result

def check_args(args):

    if args.header:
        print ('Filename\t' + create_header(args.anal_type))
        sys.exit(0)

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    if args.anal_type == 'all':
        analyses = ['wide', 'rg', 'rxs1', 'rxs2']
    else:
        analyses = [args.anal_type]

    if args.q_ranges_filename:
        q_ranges = process_range_file(args.q_ranges_filename)
    else:
        q_ranges = {}

        if analyses[0] == 'wide':
            # No fitting is performed for the wide angle plot
            q_ranges[analyses[0]] = {'qmin':args.plotrange[0],
                                    'qmax':args.plotrange[1],
                                    'fitmin':None,
                                    'fitmax':None}
        else:
            q_ranges[analyses[0]] = {'qmin':args.plotrange[0],
                                    'qmax':args.plotrange[1],
                                    'fitmin':args.fitrange[0],
                                    'fitmax':args.fitrange[0]}

    return analyses, q_ranges

def get_paths_output(input_filename):

    # The full pathname of the input is added to the output data.
    # The filename without suffix is used to title graphs
    # and name output graph files
    in_full_path = os.path.abspath(input_filename)
    in_path,filename = os.path.split(in_full_path)
    out_prefix, ext = os.path.splitext(filename)

    return in_full_path, out_prefix

def main():

    # Interpret command line arguments
    args = parse_arguments()

    analyses, q_ranges = check_args(args)

    # The full pathname of the input is added to the output data.
    # The filename without suffix is used to title graphs
    # and name output graph files
    in_full_path, out_prefix = get_paths_output(args.input_filename)

    # Initialise a string to contain computed values for output
    out_values = ''

    input_data = loadtxt(in_full_path, skiprows=args.skip)

    for cur_anal in analyses:

        # Different functions of Q I are plot/fit for each analysis:
        # Wide angle: Q vs ln(I)
        # Rg: Q^2 vs ln(I)
        # Rxs1/Rxs2: Q^2 vs ln(I*Q)

        if cur_anal == 'wide':
            xlabel = 'Q'
            x = input_data[:,0]
        else:
            xlabel = 'Q^2'
            x = input_data[:,0]**2

        if (cur_anal == 'wide') or (cur_anal == 'rg'):
            ylabel = 'ln(I(Q))'
            y = log(input_data[:,1])
        else:
            ylabel = 'ln(I(Q)*Q)'
            y = log(input_data[:,1] * input_data[:,0])

        # Range of the y-axes to plot
        if args.y_range == None:
            y_max = round(max(y) + 1, 0)
            y_min = round(min(y) - 1, 0)
        else:
            y_min, y_max = args.y_range

        q_min = q_ranges[cur_anal]['qmin']
        q_max = q_ranges[cur_anal]['qmax']

        if cur_anal == 'wide':

            # We don't perform any fitting on Wide Angle plots
            title_graph = 'Wide Angle ' + out_prefix

            fname = create_output_name(out_prefix, cur_anal,
                                       q_min, q_max, None, None)
            out_file = os.path.join(args.output_dir, fname)

            create_figure(out_file, x, y, title_graph,
                          xlabel, ylabel, q_min, q_max, y_min, y_max)

        else:

            fit_min = q_ranges[cur_anal]['fitmin']
            fit_max = q_ranges[cur_anal]['fitmax']

            fname = create_output_name(out_prefix, cur_anal,
                                       q_min, q_max, fit_min, fit_max)
            out_file = os.path.join(args.output_dir, fname)

            # create mask to select range of Q values for fitting
            mask = (input_data[:,0] > fit_min) & (input_data[:,0] < fit_max)

            # The plots here use Q^2, convert the range accordingly
            q2_min = q_min**2
            q2_max = q_max**2

            fit_text = '{0:s} (Q fit range = {1:5.4}-{2:5.4})'.format(
                                                out_prefix, fit_min, fit_max)

            # Perform linear fit over the fit range selected by mask
            results = guinier_fit(x[mask], y[mask], cur_anal)

            # calculate R? * Q for the ends of the range of Q used in fit
            rq_min = results['r'] * fit_min
            rq_max = results['r'] * fit_max

            # Create output for text and graph labels
            if cur_anal == 'rg':
                # Format data for text output
                out_values += '{0:0.2f}\t{1:0.2f}\t{2:0.2f}\t{3:0.2f}\t'.format(results['r'], rq_min, rq_max, results['i'])

                # Format data for output on graph for Rg
                title_graph = 'Rg ' + fit_text
                data_graph = 'Rg: {0:0.2f}\tIo: {1:0.2f}'.format(results['r'], results['i'])
                rq_range = 'Rg*Qmin: {0:0.2f}\tRg*Qmax: {1:0.2f}'.format(
                                                             rq_min, rq_max)
            else:
                # Format data for text output
                out_values += '{0:0.2f}\t{1:0.2f}\t{2:0.2f}\t'.format(
                                              results['r'], rq_min, rq_max)

                if cur_anal == 'rxs1':
                    # Format data for output on graph for Rxs1
                    title_graph = 'Rxs1 ' + fit_text
                    data_graph = 'Rxs1: {0:0.2f}'.format(results['r'])
                    rq_range = 'Rxs1*Qmin: {0:0.2f}\tRxs1*Qmax: {1:0.2f}'.format(rq_min,rq_max)
                else:
                    # Format data for output on graph for Rxs2
                    title_graph = 'Rxs2 ' + fit_text
                    data_graph = 'Rxs2: {0:0.2f}'.format(results['r'])
                    rq_range = 'Rxs2*Qmin: {0:0.2f}\tRxs2*Qmax: {1:0.2f}'.format(rq_min,rq_max)

            # Create graph including highlighted area for fit and key data
            create_figure(out_file, x, y, title_graph ,
                             xlabel, ylabel, q2_min, q2_max, y_min, y_max,
                             fitcoeffs = results['fit'], outputs = data_graph,
                             rqrange = rq_range, mask = mask)

    print (in_full_path + '\t'+ out_values)

if __name__ == "__main__":
    main()
