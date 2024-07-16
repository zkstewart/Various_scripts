#! python3
# plot_diffratio.py
# Script to create visualisations of the difference ratio statistics
# for assessing hypotheses of variant segregation along chromosomes.

import os, argparse, sys, re, pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.edistFile):
        eprint(f'I am unable to locate the euclidean distance file ({args.edistFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        eprint(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.width < 0:
        eprint("width must be a positive integer")
        quit()
    if args.height < 0:
        eprint("height must be a positive integer")
        quit()
    if args.power < 1:
        eprint("power must be an integer >= 1")
        quit()
    if args.bulkAlleles != []:
        if len(args.bulkAlleles) != 2:
            eprint("bulkAlleles must be a list of two integers")
            quit()
        if any([ not val >= 2 for val in args.bulkAlleles ]):
            eprint("bulkAlleles values must be integers >= 2")
            quit()
    if args.bulkOccurrence != None:
        if 0 < args.bulkOccurrence >= 1:
            eprint("bulkOccurrence must be a float value >0 and <=1")
    if args.reportAboveCutoff != None:
        if 0 < args.reportAboveCutoff >= 1:
            eprint("reportAboveCutoff must be a float value >0 and <=1")
            quit()
    # Check for conflicting arguments
    if args.onePlot and args.regions != []:
        eprint("You can't provide both --onePlot and --regions; please choose one and try again.")
        quit()
    if args.bulkAlleles != [] and args.bulkOccurrence == None:
        eprint("You must provide a --bulkOccurrence value if you provide --bulkAlleles; please try again.")
        quit()
    if args.bulkAlleles == [] and args.bulkOccurrence != None:
        eprint("You must provide --bulkAlleles if you provide --bulkOccurrence; please try again.")
        quit()
    # Handle regions
    for region in args.regions:
        if not re.match(r"^.+:\d+:\d+$", region):
            eprint(f"Region '{region}' is not in the expected format (contig:start:end)")
            eprint("Please provide regions in the format 'contig:start:end' and try again.")
            quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        eprint(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        eprint("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        eprint(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")
    
    # Handle mode-specific arguments
    if args.mode == "line":
        if args.wmaSize < 1:
            eprint("wmaSize must be an integer >= 1")
            quit()
        if args.linewidth < 1:
            eprint("linewidth must be an integer >= 1")
            quit()
    elif args.mode == "histogram":
        if args.binSize < 1:
            eprint("binSize must be an integer >= 1")
            quit()
        if 0 > args.binThreshold or args.binThreshold > 1:
            eprint("binThreshold must be a float value >=0 and <=1")
            quit()

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

def lowess_smoothing(x, y, frac=0.3, it=1, delta=0.01):
    """
    Slower and less effective than WMA() for our purposes, so I'm keeping this here
    for reference but not using it in the main program.
    """
    from statsmodels.nonparametric.smoothers_lowess import lowess
    smoothedY = lowess(y, x, frac=frac, it=it, delta=delta,
                       is_sorted=True, return_sorted=False)
    return smoothedY

def smoothTriangle(data, degree):
    """
    See https://plotly.com/python/smoothing/.
    Slower than WMA() and provides very similar output, so I've opted for that instead.
    But I am keeping this here since it's a good implementation of a triangle smoothing function
    and I might want to use it in the future.
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

def WMA(s, period):
    """
    See https://stackoverflow.com/questions/74518386/improving-weighted-moving-average-performance
    
    Parameters:
        s -- a numpy array of values to smooth
        period -- an integer value indicating the number of previous values to consider
                  during weighted moving average calculation
    Returns:
        sw -- a pandas Series of the smoothed values
    """
    w = np.arange(period)+1
    w_s = w.sum()
    
    try:
        swv = np.lib.stride_tricks.sliding_window_view(s.flatten(), window_shape=period)
    except ValueError:
        "Less data points than period size causes this error"
        return None
    sw = (swv * w).sum(axis=1) / w_s
    
    # Need to now return it as a normal series
    sw = np.concatenate((np.full(period - 1, np.nan), sw))
    try:
        sw[0:period] = sw[period] # set first n=period values to be same as first smoothed value
    except:
        "len(sw)==1 causes this error"
        return None
    return pd.Series(sw)

def lineplot_per_contig(dotsX, dotsY, wmaSize, width, height, power,
                        outputDirectory, plotPDF=False, showDots=True, linewidth=1):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the euclidean distance value for each dot; can have
                 had power transformation applied in advance
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        power -- an integer value indicating what power euclidean distances were raised to
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        showDots -- OPTIONAL; a boolean indicating whether to show dots for each data point
                    in addition to the line plot (default=True)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
    '''
    numContigsPlotted = 0
    for contigID in dotsX.keys():
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
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contigID}' has too few data points to smooth; skipping...")
            continue
        
        # Configure plot
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position (Mbp)", fontweight="bold")
        ax.set_ylabel(f"Euclidean distance (to power {power})", fontweight="bold")
        ax.set_title(f"{contigID} Euclidean distance plot", fontweight="bold")
        
        # Plot dots (if applicable)
        if showDots:
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        ax.plot(x, smoothedY, zorder=1, linewidth=linewidth)
        
        # Save output file
        plt.savefig(fileOut)
        numContigsPlotted += 1
    return numContigsPlotted

def lineplot_horizontal(dotsX, dotsY, wmaSize, width, height, power,
                        outputDirectory, plotPDF=False, showDots=True, linewidth=1):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the euclidean distance value for each dot; can have
                 had power transformation applied in advance
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        power -- an integer value indicating what power euclidean distances were raised to
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        showDots -- OPTIONAL; a boolean indicating whether to show dots for each data point
                    in addition to the line plot (default=True)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
    '''
    numContigsPlotted = 0
    
    # Derive our output file name and error if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    fileOut = os.path.join(outputDirectory, f"onePlot.{fileSuffix}")
    if os.path.isfile(fileOut):
        raise FileExistsError(f"'onePlot.{fileSuffix}' already found in output directory")
    
    # Get the ordered contig IDs
    contigIDs = get_sorted_contig_ids(dotsX.keys())
    
    # Get each contigs' plot data
    plotData = []
    for contigID in contigIDs:
        # Skip if we found no SNPs on this contig
        if not contigID in dotsY:
            print(f"WARNING: '{contigID}' is in the Euclidean distance file but has no SNPs associated " +
                    "with it; skipping...")
            continue
        
        # Get plotting values
        x = np.array(dotsX[contigID]) / 1000000 # convert to Mbp
        y = np.array(dotsY[contigID])
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contigID}' has too few data points to smooth; skipping...")
            continue
        
        # Store dot values
        plotData.append([contigID, x, y, smoothedY])
        numContigsPlotted += 1

    # Produce the figure axes
    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    gs = fig.add_gridspec(1, len(plotData), hspace=0)
    axes = gs.subplots(sharey='row')
    
    ## Set the figure title
    fig.suptitle(f"Euclidean distance plot", fontweight="bold")
    fig.supxlabel(f"Chromosomal position (Mbp)", fontweight="bold")
    fig.supylabel(f"Euclidean distance (to power {power})", fontweight="bold")
    
    # Plot the data into each axis
    for ax, (contigID, x, y, smoothedY) in zip(axes, plotData):
        # Set plot title
        ax.set_title(contigID)
        
        # Plot dots (if applicable)
        if showDots:
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        ax.plot(x, smoothedY, zorder=1, linewidth=linewidth)
    
    for ax in fig.get_axes():
        ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    
    return numContigsPlotted

def lineplot_regions(dotsX, dotsY, regions, wmaSize,
                     width, height, power, outputDirectory,
                     showDots=True, plotPDF=False, linewidth=1):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the euclidean distance value for each dot; can have
                 had power transformation applied in advance
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        power -- an integer value indicating what power euclidean distances were raised to
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        showDots -- OPTIONAL; a boolean indicating whether to show dots for each data point
                    in addition to the line plot (default=True)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
    '''
    # Parse out regions for plotting
    regions = [ region.split(":") for region in regions ]
    
    # Convert start and end values to int
    regions = [ (contigID, int(start), int(end)) for contigID, start, end in regions ]
    
    # Plot each region
    numContigsPlotted = 0
    for contigID, start, end in regions:        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contigID}.{start}_to_{end}.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Plot for '{contigID, start, end}' already found in output directory; skipping...")
            continue
        
        # Skip if we found no SNPs on this contig
        if not contigID in dotsY:
            print(f"WARNING: '{contigID}' is in the Euclidean distance file but has no SNPs associated " +
                    "with it; skipping...")
            continue
        
        # Get values within this region
        regionValues = [ [x, y] for x, y in zip(dotsX[contigID], dotsY[contigID]) if x >= start and x <= end ]
        x = np.array([ x / 1000000 for x, y in regionValues ]) # convert to Mbp
        y = np.array([ y for x, y in regionValues ])
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contigID, start, end}'  has too few data points to smooth; skipping...")
            continue
        
        # Configure plot
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position (Mbp)", fontweight="bold")
        ax.set_ylabel(f"Euclidean distance (to power {power})", fontweight="bold")
        ax.set_title(f"{contigID} Euclidean distance plot", fontweight="bold")
        
        # Plot dots (if applicable)
        if showDots:
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        ax.plot(x, smoothedY, zorder=1, linewidth=linewidth)
        
        # Save output file
        plt.savefig(fileOut)
        numContigsPlotted += 1
    return numContigsPlotted

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
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-d", dest="edistFile",
                    required=True,
                    help="Specify the location of the input euclidean distance file")
    p.add_argument("-f", dest="genomeFasta",
                    required=True,
                    help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                    required=True,
                    help="Output directory where plot files will be written")
    ## Opts (metadata behaviour)
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
    ## Opts (plotting behaviour)
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
    ## Opts (statistical behaviour)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                    type=int,
                    required=False,
                    help="""Optionally, specify the minimum size of contig to
                    create plots for (default=200000 i.e., 2Mb)""",
                    default=200000)
    p.add_argument("--power", dest="power",
                    type=int,
                    required=False,
                    help="""Optionally, specify the power to raise euclidean distances to
                    reduce noise (default=4)""",
                    default=4)
    ## Opts (output)
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
                    help="""Optionally, indicate a cutoff value for which any power-transformed
                    variant ('line' mode) or bin quantity ('histogram' mode) above this value will
                    be reported to the console""",
                    default=None)
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser()
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    lineparser = subparsers.add_parser("line",
                                       parents=[p],
                                       add_help=False,
                                       help="Create line plots of euclidean distances")
    lineparser.set_defaults(func=linemain)
    
    histoparser = subparsers.add_parser("histogram",
                                        aliases=["histo"],
                                        parents=[p],
                                        add_help=False,
                                        help="Create histograms of euclidean distances")
    histoparser.set_defaults(func=histomain)
    
    # Line-subparser arguments
    lineparser.add_argument("--wmaSize", dest="wmaSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the number of previous values to consider
                            during weighted moving average calculation (default=5)""",
                            default=5)
    lineparser.add_argument("--linewidth", dest="linewidth",
                            type=int,
                            required=False,
                            help="""Optionally, specify the line width (default=1)""",
                            default=1)
    lineparser.add_argument("--showDots", dest="showDots",
                            required=False,
                            action="store_true",
                            help="""Optionally, provide this flag if you want plots to show dots
                            of each data point in addition to the line plot""",
                            default=False)
    
    # Histogram-subparser arguments
    histoparser.add_argument("--binSize", dest="binSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the bin size to count variants
                            within (default=1000)""",
                            default=10000)
    histoparser.add_argument("--binThreshold", dest="binThreshold",
                            type=float,
                            required=False,
                            help="""Optionally, specify the bin size to count variants
                            within (default=1000)""",
                            default=10000)
    
    args = subParentParser.parse_args()
    validate_args(args)
    
    # Figure out what our pickle file should be called
    pickleFile = os.path.join(
        args.outputDirectory,
        f"{os.path.basename(args.edistFile)}.al{'_'.join(map(str, args.bulkAlleles))}.oc{args.bulkOccurrence}.pkl"
    )
    
    # Load pickle if it exists to skip computation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            dotsX, dotsY = pickle.load(fileIn)
    
    # Otherwise ...
    else:
        # Parse euclidean distance data
        dotsX, dotsY = get_edist_for_dotting(args.edistFile, args.bulkAlleles, args.bulkOccurrence)
        
        # Save data
        with open(pickleFile, "wb") as fileOut:
            pickle.dump([dotsX, dotsY], fileOut)
    
    # Power-transform values
    powerY = {}
    for contigID in dotsX.keys():
        y = np.array(dotsY[contigID])**args.power
        powerY[contigID] = y
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    # Drop any contigs which don't meet our length cutoff
    for record in genomeRecords:
        if len(record) < args.minimumContigSize:
            print(f"NOTE: '{record.id}' is below the minimum contig size and will be skipped")
            try:
                del dotsX[record.id]
                del powerY[record.id]
            except:
                raise ValueError(f"ERROR: '{record.id}' was not found in the Euclidean distance file but " +
                                 "was found in the genome FASTA file; this is unexpected and suggests " +
                                 "a mismatch between the two files")
    
    # Check that we still have contigs to plot
    if dotsX == {}:
        raise ValueError("ERROR: We didn't find any contigs which matched or exceeded the minimum size. " +
                         "Hence, no output files have been generated! Maybe you should fix your " +
                         "--minimum_contig value?")
    
    # Split into mode-specific functions
    if args.mode == "line":
        linemain(args, dotsX, powerY, lengthsDict)
    elif args.mode == "histogram":
        histomain(args, dotsX, powerY, lengthsDict)
    
def linemain(args, dotsX, powerY, lengthsDict):
    # Report any variants above the cutoff
    if args.reportAboveCutoff != None:
        print(f"# Euclidean distance >= {args.reportAboveCutoff} report:")
        for contigID in dotsX.keys():
            for x, y in zip(dotsX[contigID], powerY[contigID]):
                if y >= args.reportAboveCutoff:
                    print(f"# {contigID}:{x} = {y}")
    
    # Create plots
    if args.onePlot:
        numContigsPlotted = lineplot_horizontal(dotsX, powerY, args.wmaSize,
                                                args.width, args.height, args.power,
                                                args.outputDirectory, args.plotPDF,
                                                args.showDots, args.linewidth)
    elif args.regions != []:
        numContigsPlotted = lineplot_regions(dotsX, powerY, args.regions, args.wmaSize,
                                             args.width, args.height, args.power,
                                             args.outputDirectory, args.plotPDF,
                                             args.showDots, args.linewidth)
    else:
        numContigsPlotted = lineplot_per_contig(dotsX, powerY, args.wmaSize,
                                                args.width, args.height, args.power,
                                                args.outputDirectory, args.plotPDF,
                                                args.showDots, args.linewidth)
    
    # Raise errors if necessary
    if numContigsPlotted == 0:
        raise ValueError("ERROR: We ended up skipping every contig! This means the program has " + 
                         "already run to completion previously. Hence, no new output files have been " +
                         "generated! Maybe you should delete the existing files to restart?")
    
    print("Program completed successfully!")

def histomain(args, dotsX, powerY, lengthsDict):
    raise NotImplementedError("Histogram mode is not yet implemented")

if __name__ == "__main__":
    main()
