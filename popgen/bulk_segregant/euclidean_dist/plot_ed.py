#! python3
# plot_diffratio.py
# Script to create visualisations of the difference ratio statistics
# for assessing hypotheses of variant segregation along chromosomes.

import os, argparse, math, re, pickle
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.edistFile):
        print(f'I am unable to locate the euclidean distance file ({args.edistFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        print(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.width < 0:
        print("width must be a positive integer")
        quit()
    if args.height < 0:
        print("height must be a positive integer")
        quit()
    if args.linewidth < 1:
        print("linewidth must be an integer >= 1")
        quit()
    if args.power < 1:
        print("power must be an integer >= 1")
        quit()
    if args.bulkAlleles != []:
        if len(args.bulkAlleles) != 2:
            print("bulkAlleles must be a list of two integers")
            quit()
        if any([ not val >= 2 for val in args.bulkAlleles ]):
            print("bulkAlleles values must be integers >= 2")
            quit()
    if args.bulkOccurrence != None:
        if 0 < args.bulkOccurrence >= 1:
            print("bulkOccurrence must be a float value >0 and <=1")
    if args.reportAboveCutoff != None:
        if 0 < args.reportAboveCutoff >= 1:
            print("reportAboveCutoff must be a float value >0 and <=1")
            quit()
    # Check for conflicting arguments
    if args.onePlot and args.regions != []:
        print("You can't provide both --onePlot and --regions; please choose one and try again.")
        quit()
    if args.bulkAlleles != [] and args.bulkOccurrence == None:
        print("You must provide a --bulkOccurrence value if you provide --bulkAlleles; please try again.")
        quit()
    if args.bulkAlleles == [] and args.bulkOccurrence != None:
        print("You must provide --bulkAlleles if you provide --bulkOccurrence; please try again.")
        quit()
    # Handle regions
    for region in args.regions:
        if not re.match(r"^.+:\d+:\d+$", region):
            print(f"Region '{region}' is not in the expected format (contig:start:end)")
            print("Please provide regions in the format 'contig:start:end' and try again.")
            quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def get_edist_for_dotting(edistFile, bulkAlleles=[], bulkOccurrence=None):
    '''
    Parameters:
        edistFile -- a string pointing to the euclidean distance file
                         containing relevant statistics
        bulkAlleles -- OPTIONAL; a list of two integers indicating the maximum number of alleles
                       in each bulk OR an empty list OR None to indicate that bulk occurrence
                       should not be considered (default=[])
        bulkOccurrence -- OPTIONAL; a float value indicating the minimum fraction of occurrence
                          for one of the two bulks to be considered for plotting OR None to indicate
                          that bulk occurrence should not be considered (default=None)
    Returns:
        dotsX -- a dict pairing chromosome IDs (keys) to a list of integers indicating
                 the position for each SNP (values)
        dotsY -- a dict pairing chromosome IDs (keys) to a list of floats indicating the
                 euclidean distance for each SNP (values)
    '''
    HEADER_VALUES = ["CHROM", "POSI", "euclideanDist", "bulk1_alleles", "bulk2_alleles"]
    
    # Check that bulk values make sense
    if bulkAlleles != [] and bulkOccurrence != None:
        assert all([ isinstance(val, int) for val in bulkAlleles ]), "Bulk alleles must be integers!"
        assert all([ val >= 2 for val in bulkAlleles ]), "Bulk alleles must be >= 2!"
        assert isinstance(bulkOccurrence, float) or isinstance(bulkOccurrence, int), \
            "Bulk occurrence must be an int or float!"
        assert 0 < bulkOccurrence <= 1, "Bulk occurrence must be >0 and <=1!"
        shouldBulk = True
    else:
        shouldBulk = False
    
    # Iterate through file and grab distance values
    dotsX, dotsY = {}, {}
    with open(edistFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if firstLine:
                assert all([ hv in sl for hv in HEADER_VALUES ]), "Header line doesn't contain expected values!"
                
                chromIndex = sl.index(HEADER_VALUES[0])
                posIndex = sl.index(HEADER_VALUES[1])
                distIndex = sl.index(HEADER_VALUES[2])
                bulk1Index = sl.index(HEADER_VALUES[3])
                bulk2Index = sl.index(HEADER_VALUES[4])
                
                firstLine = False
            
            # Handle content lines
            else:
                # Parse out relevant details from this line
                chrom, pos, edist, bulk1, bulk2 = sl[chromIndex], sl[posIndex], \
                    sl[distIndex], sl[bulk1Index], sl[bulk2Index]
                
                # Check bulk occurrence and skip if applicable
                if shouldBulk:
                    # Parse out the bulk allele counts
                    bulk1Alleles, bulk2Alleles = bulkAlleles
                    
                    # Check if the bulk occurrence is met
                    if ((int(bulk1) / bulk1Alleles) < bulkOccurrence) and ((int(bulk2) / bulk2Alleles) < bulkOccurrence):
                        continue
                
                # Store data
                if not chrom in dotsX:
                    dotsX[chrom] = []
                    dotsY[chrom] = []
                
                dotsX[chrom].append(int(pos))
                dotsY[chrom].append(float(edist))
    
    return dotsX, dotsY

def smoothTriangle(data, degree):
    """
    See https://plotly.com/python/smoothing/
    """
    triangle=np.concatenate((np.arange(degree + 1), np.arange(degree)[::-1])) # up then down
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(np.sum(point)/np.sum(triangle))
    
    if smoothed == []:
        return None
    
    # Handle boundaries
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed

def weightedmovingaverage(Data, period):
    weighted = []
    Data = np.array(Data)
    for i in range(len(Data)):
            try:
                total = np.arange(1, period + 1, 1) # weight matrix
                matrix = Data[i - period + 1: i + 1]#, 3:4]
                matrix = np.ndarray.flatten(matrix)
                matrix = total * matrix # multiplication
                wma = (matrix.sum()) / (total.sum()) # WMA
                weighted = np.append(weighted, wma) # add to array
            except ValueError:
                pass
    return weighted

def WMA(s, period):
    import pandas as pd
    w = np.arange(period)+1
    w_s = w.sum()   
    swv = np.lib.stride_tricks.sliding_window_view(s.flatten(), window_shape=period)
    sw = (swv * w).sum(axis=1) / w_s
    
    # Need to now return it as a normal series
    sw = np.concatenate((np.full(period - 1, np.nan), sw))
    sw[0:period] = s[0:period]
    return pd.Series(sw)

def loess_smoothing(x, y, degree=1, frac=0.5):
    from loess.loess_1d import loess_1d
    xout, yout, wout = loess_1d(x, y, xnew=None, degree=degree, #frac=frac,
                                npoints=10000, rotate=False, sigy=None)
    return yout

def lowess_smoothing(x, y, frac=0.2, it=1, delta=10000):
    from statsmodels.nonparametric.smoothers_lowess import lowess
    smoothedY = lowess(y, x,
                    frac=frac, it=it,
                    delta=delta,
                    is_sorted=True, return_sorted=False)
    return smoothedY

def pygam_smoothing(x, y, npoints=1000):
    from pygam import LinearGAM, s, f
    gam = LinearGAM(s(0) ).fit(np.array(x), np.array(y))
    xout = gam.generate_X_grid(term=0, n=npoints).flatten()
    yout = gam.predict(xout)
    return xout, yout

def pygam_expectiles(x, y, npoints=1000):
    from pygam import ExpectileGAM
    x = np.array(x)
    x = x.reshape(len(x), 1)
    gam = ExpectileGAM(expectile=0.5).gridsearch(x, np.array(y))
    
    lam = gam.lam
    
    gam95 = ExpectileGAM(expectile=0.95, lam=lam).fit(x, np.array(y))
    gam75 = ExpectileGAM(expectile=0.75, lam=lam).fit(x, np.array(y))
    gam25 = ExpectileGAM(expectile=0.25, lam=lam).fit(x, np.array(y))
    gam05 = ExpectileGAM(expectile=0.05, lam=lam).fit(x, np.array(y))
    
    XX = gam.generate_X_grid(term=0, n=500)
    
    plt.scatter(x, y, c='k', alpha=0.2)
    plt.plot(XX, gam95.predict(XX), label='0.95')
    plt.plot(XX, gam75.predict(XX), label='0.75')
    plt.plot(XX, gam.predict(XX), label='0.50')
    plt.plot(XX, gam25.predict(XX), label='0.25')
    plt.plot(XX, gam05.predict(XX), label='0.05')
    plt.legend()
    
    xout = gam.generate_X_grid(term=0, n=npoints).flatten()
    yout = gam.predict(xout)
    return xout, yout

def alr(x, y):
    import sys
    sys.path.append(r"C:\bio\ALR\src")
    from ALR import Automated_Loess_Regression
    x, y = xList, yList
    
    alrResult = Automated_Loess_Regression(np.array(x), np.array(y), err_y=0, deg=2, alpha=0,
                                           outliers_det=False, n_sims=1000, average=True, verbose=False)

def plot_per_contig(dotsX, dotsY, powerY, lengthsDict, minimumContigSize, width, height, outputDirectory,
                    plotPDF=False, linewidth=1):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the euclidean distance value for each dot
        powerY -- a dictionary linking chromosome IDs (keys) to lists of floats
                  indicating the euclidean distance value for each position after
                  power transformation
        lengthsDict -- a dictionary with structure like:
                       {
                            'contig1': intLength1,
                            'contig2': intLength2,
                            ...
                       }
        minimumContigSize -- an integer value indicating the minimum size a contig must be
                             to be considered for plotting
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
    '''
    numContigsProcessed = 0
    numContigsPlotted = 0
    for contigID, length in lengthsDict.items():
        if length >= minimumContigSize:
            numContigsProcessed += 1
            
            # Derive our output file name and skip if already existing
            fileSuffix = "pdf" if plotPDF else "png"
            fileOut = os.path.join(outputDirectory, f"{contigID}.{fileSuffix}")
            if os.path.isfile(fileOut):
                print(f"WARNING: Plot for '{contigID}' already found in output directory; skipping...")
                continue
            
            # Skip if we found no SNPs on this contig
            if not contigID in dotsY:
                print(f"WARNING: '{contigID}' is in the Euclidean distance file but has no SNPs associated " +
                      "with it; skipping...")
                continue
            
            # Get plotting values
            x = np.array(dotsX[contigID]) / 1000000 # convert to Mbp
            y = np.array(dotsY[contigID])
            power = powerY[contigID]
            #smoothedList = smoothedY[contigID]
            smoothedList = smoothTriangle(power, 10)
            #smoothedList = WMA(np.array(yList), 100)
            #smoothedList = loess_smoothing(xList, np.array(yList), degree=1, frac=0.5)
            #smoothedX, smoothedList = pygam_smoothing(xList, yList, npoints=100000)
            
            # Skip plotting if smoothing fails
            "This probably means there are not enough data points to smooth"
            if smoothedList == None:
                print(f"WARNING: '{contigID}' has too few data points to smooth; skipping...")
                continue
            
            # Configure plot
            fig = plt.figure(figsize=(width, height))
            ax = plt.axes()
            
            ax.set_xlabel(f"Chromosomal position (Mbp)", fontweight="bold")
            ax.set_ylabel("Euclidean distance", fontweight="bold")
            ax.set_title(f"{contigID} euclidean distance plot", fontweight="bold")
            
            # Plot dots
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
            
            # Plot line
            ax.plot(x, smoothedList, zorder=1, linewidth=linewidth)
            #ax.plot(smoothedX, smoothedList, zorder=1, linewidth=linewidth)
            
            # Save output file
            plt.savefig(fileOut)
            numContigsPlotted += 1
    return numContigsProcessed, numContigsPlotted

def plot_once(dotsX, dotsY, lengthsDict, minimumContigSize, width, height, outputDirectory,
              plotPDF=False, linewidth=1):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the euclidean distance value for each dot
        lengthsDict -- a dictionary with structure like:
                       {
                            'contig1': intLength1,
                            'contig2': intLength2,
                            ...
                       }
        minimumContigSize -- an integer value indicating the minimum size a contig must be
                             to be considered for plotting
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
    '''
    numContigsProcessed = 0
    numContigsPlotted = 0
    
    # Derive our output file name and error if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    fileOut = os.path.join(outputDirectory, f"onePlot.{fileSuffix}")
    if os.path.isfile(fileOut):
        raise FileExistsError(f"'onePlot.{fileSuffix}' already found in output directory")
    
    # Calculate kbp values for plot axis labels
    kbpWindowSize = round(windowSize / 1000, 2)
    kbpStepSize = round(stepSize / 1000, 2)
    
    # Get the ordered contig IDs
    contigIDs = get_sorted_contig_ids(lengthsDict.keys())
    
    # Get each contigs' plot data
    plotData = []
    for contigID in contigIDs:
        length = lengthsDict[contigID]
        if length >= minimumContigSize:
            numContigsProcessed += 1
            
            # Skip if we found no SNPs on this contig
            if not contigID in densityDict:
                print(f"WARNING: '{contigID}' is in the difference ratio file but has no SNPs associated " +
                      "with it; skipping...")
                continue
            
            # Get the dots for this contig (if applicable)
            if dotsX != None and contigID in dotsX:
                dotsXList = dotsX[contigID]
                dotsYList = dotsY[contigID]
            else:
                dotsXList, dotsYList = None, None
            
            # Store density values
            plotData.append([densityDict[contigID], contigID, dotsXList, dotsYList])
            numContigsPlotted += 1
    
    # Format a single joined plot
    ## Produce the figure axes
    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    gs = fig.add_gridspec(1, len(plotData), hspace=0)
    axes = gs.subplots(sharey='row')
    
    ## Set the figure title
    fig.suptitle(f"Contigs difference ratio plot", fontweight="bold")
    fig.supxlabel(f"Chromosomal position ({kbpWindowSize} kbp windows; {kbpStepSize} kbp step size)", fontweight="bold")
    fig.supylabel("Mean difference ratio", fontweight="bold")
    
    ## Plot the data into each axis
    for ax, (densityList, contigID, dotsXList, dotsYList) in zip(axes, plotData):
        ax.set_title(contigID)
        
        ### Plot dots (if applicable)
        if dotsXList != None and dotsYList != None:
            ax.scatter(dotsXList, dotsYList, color="red", s=3, alpha=0.5, zorder=0)
        
        ### Plot line
        ax.plot(densityList, zorder=1, linewidth=linewidth)
    
    for ax in fig.get_axes():
        ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    
    return numContigsProcessed, numContigsPlotted

def plot_regions(dotsX, dotsY, regions, width, height, outputDirectory,
                 plotPDF=False, linewidth=1):
    '''
    Parameters:
        densityDict -- a dictionary with structure like:
                       {
                           'contig1': [ floatDensity1, floatDensity2, ... ],
                           'contig2': [ floatDensity1, floatDensity2, ... ],
                       }, where each list contains a number of float values equal to
                       the length of the contig divided by its step size
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        windowSize -- an integer value indicating what length of genome the difference ratio
                      values were averaged over
        stepSize -- an integer value indicating the step size used to move across the genome
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        dotsX -- OPTIONAL; a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the step position where a dot is located OR None for no dots
        dotsY -- OPTIONAL; a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the difference ratio value for each dot OR None for no dots
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
    '''
    # Parse out regions for plotting
    regions = [ region.split(":") for region in regions ]
    
    # Convert regions to step indices
    regions = [
        [ contig, math.floor(int(start) / stepSize), math.floor(int(end) / stepSize), start, end ]
        for contig, start, end in regions
    ]
    
    # Plot each region
    numContigsPlotted = 0
    numContigsProcessed = 0
    for contigID, startStep, endStep, origStart, origEnd in regions:
        numContigsProcessed += 1
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contigID}.{origStart}_to_{origEnd}.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Plot for '{contigID, origStart, origEnd}' already found in output directory; skipping...")
            continue
        
        # Skip if we found no SNPs on this contig
        if not contigID in densityDict:
            print(f"WARNING: '{contigID}' is in the difference ratio file but has no SNPs associated " +
                    "with it; skipping...")
            continue
        
        # Get density values for the region
        xList = [ i for i in range(startStep, endStep+1) ]
        yList = densityDict[contigID][startStep:endStep+1]
        
        # Get the dots for this region (if applicable)
        if dotsX != None and contigID in dotsX:
            dotsXList = dotsX[contigID]
            dotsYList = dotsY[contigID]
            
            # Subset to the region of interest
            regionDots = [
                (x, y)
                for x, y in zip(dotsXList, dotsYList)
                if x >= startStep and x <= endStep
            ]
            
            # Unpack the region dots back into X and Y lists
            dotsXList = [ x for x, y in regionDots ]
            dotsYList = [ y for x, y in regionDots ]
        
        else:
            dotsXList, dotsYList = None, None
        
        # Configure plot
        kbpWindowSize = round(windowSize / 1000, 2)
        kbpStepSize = round(stepSize / 1000, 2)
        
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows; {kbpStepSize} kbp step size)", fontweight="bold")
        ax.set_ylabel("Mean difference ratio", fontweight="bold")
        ax.set_title(f"{contigID} difference ratio plot", fontweight="bold")
        
        # Plot dots (if applicable)
        if dotsXList != None and dotsYList != None:
            ax.scatter(dotsXList, dotsYList, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        ax.plot(xList, yList, zorder=1, linewidth=linewidth)
        
        # Save output file
        plt.savefig(fileOut)
        numContigsPlotted += 1
    return numContigsProcessed, numContigsPlotted

def get_sorted_contig_ids(idsList):
    # Sort contig IDs by their numerical value (if possible)
    allHaveNumbers = all([ any([ c.isdigit() for c in contigID ]) for contigID in idsList ])
    if allHaveNumbers:
        numRegex = re.compile(r"\d+")
        return sorted(idsList, key=lambda x: int("".join(numRegex.findall(x))))
    else:
        return sorted(idsList)

def main():
    usage = """%(prog)s receives a difference ratio TSV file and creates
    difference ratio plots per chromosome. It specifically calculates the
    average difference ratio value over each window and plots that value.
    Hence, it may help to visualise where in the genome regions associated
    with bulks exist.
    
    Note that step size dictates where the centre of each window will be,
    with the window size dictating the size of the window.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="edistFile",
                   required=True,
                   help="Specify the location of the input euclidean distance file")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where plot files will be written")
    # Opts (metadata behaviour)
    p.add_argument("--bulkAlleles", dest="bulkAlleles",
                   required=False,
                   nargs="+",
                   type=int,
                   help="""Optionally, indicate the number of maximum possible alleles
                   in each bulk in order to calculate the occurrence fraction for
                   filtering""",
                   default=[])
    p.add_argument("--bulkOccurrence", dest="bulkOccurrence",
                   type=float,
                   required=False,
                   help="""Optionally, specify the minimum fraction of occurrence
                   for one of the two bulks to be considered for plotting""",
                   default=None)
    # Opts (plotting behaviour)
    p.add_argument("--width", dest="width",
                   type=int,
                   required=False,
                   help="""Optionally, specify the output plot width (default=10)""",
                   default=10)
    p.add_argument("--height", dest="height",
                   type=int,
                   required=False,
                   help="""Optionally, specify the output plot height (default=6)""",
                   default=6)
    p.add_argument("--showDots", dest="showDots",
                   required=False,
                   action="store_true",
                   help="""Optionally, provide this flag if you want plots to show dots
                   of each data point in addition to the line plot""",
                   default=False)
    p.add_argument("--linewidth", dest="linewidth",
                   type=int,
                   required=False,
                   help="""Optionally, specify the line width (default=1)""",
                   default=1)
    # Opts (statistical behaviour)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the minimum size of contigs which
                   should have plots created for (default=200000)"; this value
                   must exceed --window_size by at least 2x""",
                   default=200000)
    p.add_argument("--power", dest="power",
                   type=int,
                   required=False,
                   help="""Optionally, specify the power to raise euclidean distances to
                   reduce noise (default=4)""",
                   default=4)
    # Opts (output)
    p.add_argument("--regions", dest="regions",
                   required=False,
                   nargs="+",
                   help="""Optionally, indicate one or more regions to plot in greater detail
                   by providing the contig ID and start and end positions in bp (e.g.
                   contig1:10000:20000)""",
                   default=[])
    p.add_argument("--onePlot", dest="onePlot",
                   required=False,
                   action="store_true",
                   help="""Optionally, provide this flag if you want a single plot to be
                   produced with all chromosomes positioned horizontally""",
                   default=False)
    p.add_argument("--pdf", dest="plotPDF",
                   required=False,
                   action="store_true",
                   help="""Optionally, provide this flag if you want outputs to be
                   in PDF format instead of PNG format""",
                   default=False)
    p.add_argument("--reportAboveCutoff", dest="reportAboveCutoff",
                   required=False,
                   type=float,
                   help="""Optionally, indicate a cutoff value for which any window with a
                   difference ratio above this value will be reported to the console""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Figure out what our pickle file should be called
    pickleFile = os.path.join(
        args.outputDirectory,
        f"{os.path.basename(args.edistFile)}.al{'_'.join(map(str, args.bulkAlleles))}.oc{args.bulkOccurrence}.po{args.power}.pkl"
    )
    
    # Load pickle if it exists to skip computation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            dotsX, dotsY, powerY = pickle.load(fileIn)
    
    # Otherwise ...
    else:
        # Parse euclidean distance data
        dotsX, dotsY = get_edist_for_dotting(args.edistFile, args.bulkAlleles, args.bulkOccurrence)
        
        # Store power-transformed values
        powerY = {}
        for contigID in dotsX.keys():
            y = np.array(dotsY[contigID])**args.power
            powerY[contigID] = y
            
        # Save the smoothed data
        with open(pickleFile, "wb") as fileOut:
            pickle.dump([dotsX, dotsY, powerY], fileOut)
    
    # Create plots
    if args.onePlot:
        numContigsProcessed, numContigsPlotted = plot_once(dotsX, powerY,
                                                           lengthsDict, args.minimumContigSize,
                                                           args.width, args.height,
                                                           args.outputDirectory, args.plotPDF,
                                                           args.linewidth)
    elif args.regions != []:
        numContigsProcessed, numContigsPlotted = plot_regions(dotsX, powerY, args.regions,
                                                              args.width, args.height,
                                                              args.outputDirectory, args.plotPDF,
                                                              args.linewidth)
    else:
        numContigsProcessed, numContigsPlotted = plot_per_contig(dotsX, dotsY, powerY,
                                                                 lengthsDict, args.minimumContigSize,
                                                                 args.width, args.height,
                                                                 args.outputDirectory, args.plotPDF,
                                                                 args.linewidth)
    
    # Raise relevant warnings
    if numContigsProcessed == 0:
        print(f"WARNING: We didn't find any contigs which exceeded {args.minimumContigSize}bp in size")
        print("Hence, no output files have been generated! Maybe you should fix your --minimum_contig value?")
    elif numContigsPlotted == 0:
        print("WARNING: We ended up skipping every contig! This means the program has already run to completion previously.")
        print("Hence, no new output files have been generated! Maybe you should delete the existing files to restart?")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
