#!/bin/env python

# David W Wright - 09/08/2013
"""Script to analsye Q and I data from SAS experiments
Can produce the following plots and data:
wide - Wide Angle:  Q vs ln(I)
rg - Radius of gyration/Guinier:     Q^2 vs ln(I)     
rxs1/rx2 - Rxs1/Rsx2:   Q^2 vs ln(I*Q)
At the minute we assume input has no header and the 
first column is Q and second is I """

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
    parser = argparse.ArgumentParser(description= 'Basic processing of SAS data. \nProduces Rg, Io, Rxs1 and Rxs2 estimates and associated plots.\n',usage=use_message)
    
    group = parser.add_mutually_exclusive_group(required=True)   
    group.add_argument('-i','--infile', nargs='?', type=str, dest='inputFilename', help = 'Path to the input file')
    group.add_argument('--header', action='store_true', help = 'Output a header alongside output data')
    
    parser.add_argument('-s','--skip', nargs = '?', type = int, default = 0, help = 'Number of lines in input file to skip')
    
    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument('-q','--q-ranges', type = str, dest='qRangesFilename', help = 'Path to file containing Q plot and fit ranges')
    group2.add_argument('-r','--plotrange', nargs = 2, type = float, help = 'Range of Q to plot')

    parser.add_argument('-f','--fitrange', nargs = 2, type = float, help = 'Range of Q to use in fitting') 
    parser.add_argument('-o','--outdir', nargs='?', type=str, default = '.', dest='outputDir', help = 'Path for output files')
    parser.add_argument('-a','--analysis', choices = ['wide', 'rg', 'rxs1', 'rxs2','all'], default = 'wide', dest='analType', help = 'Type of analysis to run')
    parser.add_argument('-y','--yrange', nargs = 2, type = float, default = [-6,-1], help = "Range of the y-axis of plots; ln(I) for 'wide' and 'rg', ln(I*Q) for 'rxs1' and 'rxs2'.")

    args = parser.parse_args()
        
    if (args.analType == 'all') and (not args.qRangesFilename):
        print ("If you require all analyses to be run a file containig the relevant Q ranges must be provided.\n"
        "Use the following format/ordering and use a 6 line header for descriptions:\n"
        "wide		rg				rxs1				rxs2\n"
        "qmin	qmax	qmin	qmax	fitmin	fitmax	qmin	qmax	fitmin	fitmax	qmin	qmax	fitmin	fitmax\n")
        sys.exit(1)
    elif (not args.qRangesFilename) and ((not args.fitrange) or (args.analType == 'wide')):
        print ("For all analyses except wide angle plotting ('wide') a fit range is needed in addition to the plot range.\n")
        sys.exit(1)
    
    return args

def process_range_file(filename):
    """Load file containing the ranges for all analyses.
    """

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

def create_figure(filename, x, y, titleText, xLab, yLab, xMin, xMax, yMin, yMax, **kwargs):
    """ Outputs graph of x and y to a pdf file
    Values computed from graph (outputs) and range of R? * Qfit written as text on graph
    """
    
    fitCoeffs = kwargs.get('fitcoeffs', None)
    outputs = kwargs.get('outputs', None)
    rqRange = kwargs.get('rqrange', None)
    
    # Plot the input x, y values and if provided a fit line
    # Also output data calculated from fit and R*Q over fit valies
    
    figure(figsize=(8, 6), dpi=80)
    ax = subplot(111, xlabel=xLab, ylabel=yLab, title=titleText, xlim=(xMin, xMax), ylim=(yMin,yMax))
    scatter(x, y, s=5)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(15)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Plot linear fit and highlight points used in its construction
    if fitCoeffs != None:
        fitLine = poly1d(fitCoeffs)
        xPoints = linspace(xMin, xMax, 300)    
        plot(xPoints, fitLine(xPoints))
        scatter(x[mask], y[mask],s=30)
        plt.annotate(outputs + '\n' + rqRange, xy=(0.45, 0.85), xycoords='axes fraction')
        
    savefig(filename)

def create_output_name(prefix, analysis, q_min, q_max, fit_min, fit_max):
    
    q_range = "-".join(str(q_min),str(q_max))
    fit_range = "-".join(str(fit_min),str(fit_max))
    
    filename = "{0:s}_q{1:s}_fit{2:s}_{3:s}.pdf".format(
                prefix, q_range, fit_range, analysis)
    
    return filename

def main():

    # Interpret command line arguments
    args = parse_arguments()
    
    if args.header:
        print ('Filename\t' + create_header(args.analType))
        sys.exit(0)
    
    if not os.path.isdir(args.outputDir):
        os.makedirs(args.outputDir)
    
    if args.analType == 'all':
        analyses = ['wide', 'rg', 'rxs1', 'rxs2']
    else:
        analyses = [args.analType]
    
    if args.qRangesFilename:
        q_ranges = process_range_file(args.qRangesFilename)
    else:
        q_ranges = {}

        qRanges = [args.plotrange]
        
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
    
    # The full pathname of the input is added to the output data.
    # The filename without suffix is used to title graphs 
    # and name output graph files
    inFullPath = os.path.abspath(args.inputFilename)
    inPath,filename = os.path.split(inFullPath)
    out_prefix, ext = os.path.splitext(filename)
    
    # Range of the y-axes to plot
    y_min, y_max = args.yrange
    
    # Initialise a string to add the computed values to for output
    outValues = ''
    
    inputData = loadtxt(inFullPath, skiprows=args.skip)
    
    for cur_anal in analyses:
    
        q_min = q_ranges[cur_anal][qmin]
        q_max = q_ranges[cur_anal][qmax]
        
        fit_min = q_ranges[cur_anal][fitmin]
        fit_max = q_ranges[cur_anal][fitmax]
        
        fname = create_output_name(out_prefix, analysis, 
                                   q_min, q_max, fit_min, fit_max)
    
        out_path = os.path.join(args.outputDir, fname)
    
    
        if curAnal == 'wide':
            
            # We don't perform any fitting on Wide Angle plots 
            x = inputData[:,0]
            y = log(inputData[:,1])
            
            create_figure(out_path, x, y, 'Wide Angle ' + out_prefix , 
                          'Q', 'ln(I(Q))', q_min, q_max, y_min, y_max)
            
        else:
            
            # A mask is created to select out the data in the range of Q values selected for fitting
            mask = (inputData[:,0] > fit_min) & (inputData[:,0] < fit_max)
            
            # The plots here use Q^2 on the x-axis so convert the min and max of the range
            q2_min = q_min**2
            q2_max = q_max**2
            
            # TODO: basically everything
            
            fitText = outPrefix + ' (Q fit range = ' + strFitQMin + '-' + strFitQMax + ')'
            
            if curAnal == 'rg':
                
                x = inputData[:,0]**2
                y = log(inputData[:,1])
                
                # Calculate Rg and Io from linear fit of Q^2 vs ln(I)
                fitCoeffs = polyfit(x[mask], y[mask],1)
                Rg = sqrt(3 * abs(fitCoeffs[0]))
                Io = exp(fitCoeffs[1])
                
                rqMin = Rg * fitQMin
                rqMax = Rg * fitQMax
                
                # Format data for text output
                outValues = outValues + '{0:0.2f}\t{1:0.2f}\t{2:0.2f}\t{3:0.2f}\t'.format(Rg, rqMin, rqMax, Io)
                
                # Create header and format data for output on graph
                dataGraph = 'Rg: {0:0.2f}\tIo: {1:0.2f}'.format(Rg, Io)
                rqRange = 'Rg*Qmin: {0:0.2f}\tRg*Qmax: {1:0.2f}'.format(rqMin, rqMax)
                create_figure(outputFilename, x, y, 'Rg ' + fitText , 'Q^2', 'ln(I(Q))', q2Min, q2Max, yMin, yMax, fitcoeffs = fitCoeffs, outputs = dataGraph, rqrange = rqRange)
                
            else:
                
                x = inputData[:,0]**2
                y = log(inputData[:,1] * inputData[:,0])
                
                # Calculate Rxs1 or Rxs2 from linear fit of Q^2 vs ln(I*Q)
                fitCoeffs = polyfit(x[mask], y[mask],1)
                Rxs = sqrt(2 * abs(fitCoeffs[0]))
                rqMin = Rxs * fitQMin
                rqMax = Rxs * fitQMax
                
                # Format data for text output
                outValues = outValues + '{0:0.2f}\t{1:0.2f}\t{2:0.2f}\t'.format(Rxs, Rxs * fitQMin, Rxs * fitQMax)
                
                if curAnal == 'rxs1':
                    
                    # Format data for output on graph           
                    dataOut = 'Rxs1: {0:0.2f}'.format(Rxs)
                    rqRange = 'Rxs1*Qmin: {0:0.2f}\tRxs1*Qmax: {1:0.2f}'.format(rqMin, rqMax)
                    
                    create_figure(outputFilename, x, y, 'Rxs1 ' + fitText, 'Q^2', 'ln(I(Q)*Q)',  q2Min, q2Max, yMin, yMax, fitcoeffs = fitCoeffs, outputs = dataGraph, rqrange = rqRange)
                    
                elif curAnal == 'rxs2':
                    
                    # Format data for output on graph           
                    dataOut = 'Rxs2: {0:0.2f}'.format(Rxs)
                    rqRange = 'Rxs2*Qmin: {0:0.2f}\tRxs2*Qmax: {1:0.2f}'.format(rqMin, rqMax)
                    
                    create_figure(outputFilename, x, y, 'Rxs2 ' + fitText, 'Q^2', 'ln(I(Q)*Q)',  q2Min, q2Max, yMin, yMax, fitcoeffs = fitCoeffs, outputs = dataGraph, rqrange = rqRange)
                    
    print (inFullPath + '\t'+ outValues)

if __name__ == "__main__":
    main()