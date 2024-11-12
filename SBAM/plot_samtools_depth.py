#! python3
# plot_samtools_depth.py
# Script to create visualisations of the depth of coverage across a genome
# as calculated by samtools depth.

import os, argparse, sys, re, pickle, math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from contextlib import nullcontext

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.depthFile):
        eprint(f'I am unable to locate the samtools depth TSV file ({args.depthFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.width < 0:
        eprint("width must be a positive integer")
        quit()
    if args.height < 0:
        eprint("height must be a positive integer")
        quit()
    if args.ylim == 0 or args.ylim < -1:
        eprint("ylim must be a positive integer (to set a limit) or -1 (to set no limit)")
        quit()
    # Check for conflicting arguments
    if args.onePlot and args.regions != []:
        eprint("You can't provide both --onePlot and --regions; please choose one and try again.")
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
        if args.linewidth < 0:
            eprint("linewidth must be a float > 0")
            quit()
    elif args.mode == "histogram":
        if args.binSize < 1:
            eprint("binSize must be an integer >= 1")
            quit()

def get_depth_for_dotting(tsvFile):
    '''
    Parameters:
        tsvFile -- a string pointing to the TSV file with format:
                   contigID    position    depth
    Returns:
        dotsX -- a dict pairing chromosome IDs (keys) to a list of integers indicating
                 each position (bp)
        dotsY -- a dict pairing chromosome IDs (keys) to a list of floats indicating the
                 depth for each position (number of reads)
    '''
    # Iterate through file and grab statistical values
    dotsX, dotsY = {}, {}
    with open(tsvFile, "r") as fileIn:
        for line in fileIn:
            # Parse out relevant details from this line
            chrom, pos, depth = line.rstrip("\r\n ").split("\t")
            
            # Store data
            if not chrom in dotsX:
                dotsX[chrom] = []
                dotsY[chrom] = []
            
            dotsX[chrom].append(int(pos))
            dotsY[chrom].append(int(depth))

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

def lineplot_per_contig(dotsX, dotsY, wmaSize, ylim, width, height,
                        outputDirectory, plotPDF=False, showDots=True,
                        linewidth=1, statisticLabel="Q13 depth", createTSV=False):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the depth value for each dot
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        ylim -- an integer value indicating the maximum value that will be presented
                OR -1 to set no limit
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        showDots -- OPTIONAL; a boolean indicating whether to show dots for each data point
                    in addition to the line plot (default=True)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
    '''
    numContigsPlotted = 0
    for contigID in dotsX.keys():
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contigID}.line.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Line plot for '{contigID}' already found in output directory; skipping...")
            continue
        
        # Skip if we found no SNPs on this contig
        if not contigID in dotsY:
            print(f"WARNING: '{contigID}' is in the statistics file but has no SNPs associated " +
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
        
        # Enforce ylim if set
        if ylim != -1:
            smoothedY[smoothedY > ylim] = ylim
        
        # Configure plot
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position (Mbp)", fontweight="bold")
        ax.set_ylabel(f"{statisticLabel}", fontweight="bold")
        ax.set_title(f"{contigID} {statisticLabel} plot", fontweight="bold")
        ax.set_ylim(0, 50)
        
        # Plot dots (if applicable)
        if showDots:
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        ax.plot(x, smoothedY, zorder=1, linewidth=linewidth)
        
        # Save output file
        plt.savefig(fileOut)
        numContigsPlotted += 1
        
        # Write TSV file if requested
        if createTSV:
            with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") as fileOutTSV:
                fileOutTSV.write("contigID\tposition\tdepth\tsmoothed_depth\n")
                for xVal, yVal, smoothedYVal in zip(x*1000000, y, smoothedY): # convert back to bp
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\t{smoothedYVal}\n")
    
    return numContigsPlotted

def lineplot_horizontal(dotsX, dotsY, wmaSize, ylim, width, height,
                        outputDirectory, plotPDF=False, showDots=True,
                        linewidth=1, statisticLabel="Q13 depth", createTSV=False):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the depth value for each dot
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        ylim -- an integer value indicating the maximum value that will be presented
                OR -1 to set no limit
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        showDots -- OPTIONAL; a boolean indicating whether to show dots for each data point
                    in addition to the line plot (default=True)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
    '''
    numContigsPlotted = 0
    
    # Derive our output file name and error if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    fileOut = os.path.join(outputDirectory, f"onePlot.line.{fileSuffix}")
    if os.path.isfile(fileOut):
        raise FileExistsError(f"'onePlot.line.{fileSuffix}' already found in output directory")
    
    # Get the ordered contig IDs
    contigIDs = get_sorted_contig_ids(dotsX.keys())
    
    # Get each contigs' plot data
    plotData = []
    for contigID in contigIDs:
        # Get plotting values
        x = np.array(dotsX[contigID]) / 1000000 # convert to Mbp
        y = np.array(dotsY[contigID])
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contigID}' has too few data points to smooth; skipping...")
            continue
        
        # Enforce ylim if set
        if ylim != -1:
            smoothedY[smoothedY > ylim] = ylim
        
        # Store dot values
        plotData.append([contigID, x, y, smoothedY])
        numContigsPlotted += 1
    
    # Produce the figure axes
    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    gs = fig.add_gridspec(1, len(plotData), hspace=0)
    axes = gs.subplots(sharey='row')
    
    ## Set the figure title
    fig.suptitle(f"{statisticLabel} plot", fontweight="bold")
    fig.supxlabel(f"Chromosomal position (Mbp)", fontweight="bold")
    fig.supylabel(f"{statisticLabel}", fontweight="bold")
    
    # Plot the data into each axis
    with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
        # Write TSV header if applicable
        if createTSV:
            fileOutTSV.write("contigID\tposition\tdepth\tsmoothed_depth\n")
        
        for ax, (contigID, x, y, smoothedY) in zip(axes, plotData):
            # Set plot title
            ax.set_title(contigID)
            
            # Plot dots (if applicable)
            if showDots:
                ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
            
            # Plot line
            ax.plot(x, smoothedY, zorder=1, linewidth=linewidth)
            
            # Write TSV file if requested
            if createTSV:
                for xVal, yVal, smoothedYVal in zip(x*1000000, y, smoothedY): # convert back to bp
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\t{smoothedYVal}\n")
        
        for ax in fig.get_axes():
            ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    
    return numContigsPlotted

def lineplot_regions(dotsX, dotsY, regions, wmaSize, ylim,
                     width, height, outputDirectory,
                     plotPDF=False, showDots=True,
                     linewidth=1, statisticLabel="Q13 depth", createTSV=False):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the statistical value for each dot
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        ylim -- an integer value indicating the maximum value that will be presented
                OR -1 to set no limit
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        showDots -- OPTIONAL; a boolean indicating whether to show dots for each data point
                    in addition to the line plot (default=True)
        linewidth -- OPTIONAL; an integer value indicating the width of the line plot (default=1)
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
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
        fileOut = os.path.join(outputDirectory, f"{contigID}.{start}_to_{end}.line.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Line plot for '{contigID, start, end}' already found in output directory; skipping...")
            continue
        
        # Skip if we found no SNPs on this contig
        if not contigID in dotsY:
            print(f"WARNING: '{contigID}' is in the statistics file but has no SNPs associated " +
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
        
        # Enforce ylim if set
        if ylim != -1:
            smoothedY[smoothedY > ylim] = ylim
        
        # Configure plot
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position (Mbp)", fontweight="bold")
        ax.set_ylabel(f"{statisticLabel}", fontweight="bold")
        ax.set_title(f"{contigID} {statisticLabel} plot", fontweight="bold")
        
        # Plot dots (if applicable)
        if showDots:
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        ax.plot(x, smoothedY, zorder=1, linewidth=linewidth)
        
        # Save output file
        plt.savefig(fileOut)
        numContigsPlotted += 1
        
        # Write TSV file if requested
        if createTSV:
            with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") as fileOutTSV:
                fileOutTSV.write("contigID\tposition\tdepth\tsmoothed_depth\n")
                for xVal, yVal, smoothedYVal in zip(x*1000000, y, smoothedY): # convert back to bp
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\t{smoothedYVal}\n")
    
    return numContigsPlotted

def histo_per_contig(histoDict, width, height, ylim,
                     outputDirectory, binSize, plotPDF=False,
                     statisticLabel="Q13 depth", createTSV=False):
    '''
    Parameters:
        histoDict -- dictionary linking chromosome IDs (keys) to numpy arrays of integers
                     containing the number of SNPs in each bin that exceeded the threshold
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        ylim -- an integer value indicating the maximum value that will be presented
                OR -1 to set no limit
        outputDirectory -- a string indicating the directory where output files will be written
        binSize -- an integer value indicating the size of the bins that were used for counting
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
    '''
    numContigsPlotted = 0
    for contigID in histoDict.keys():
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contigID}.histo.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Histogram plot for '{contigID}' already found in output directory; skipping...")
            continue
        
        # Get plotting values
        x = np.arange(0, len(histoDict[contigID]))
        y = histoDict[contigID]
        stepSize = binSize / 1000 # convert to Kbp
        
        # Enforce ylim if set
        if ylim != -1:
            y[y > ylim] = ylim
        
        # Configure plot
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        ax.set_xlabel(f"Bin number ({stepSize} Kbp step and width)", fontweight="bold")
        ax.set_ylabel(f"Number of read alignments with {statisticLabel}", fontweight="bold")
        ax.set_title(f"{contigID} {statisticLabel} histogram", fontweight="bold")
        
        # Plot bars
        ax.bar(x, y, zorder=0)
        
        # Save output file
        plt.savefig(fileOut)
        numContigsPlotted += 1
        
        # Write TSV file if requested
        if createTSV:
            with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") as fileOutTSV:
                fileOutTSV.write("contigID\tposition\tdepth\n")
                for xVal, yVal in zip(x*1000000, y): # convert back to bp
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\n")
        
    return numContigsPlotted

def histo_horizontal(histoDict, width, height, ylim, outputDirectory,
                     binSize, plotPDF=False,
                     statisticLabel="Q13 depth", createTSV=False):
    '''
    Parameters:
        histoDict -- dictionary linking chromosome IDs (keys) to numpy arrays of integers
                     containing the number of SNPs in each bin that exceeded the threshold
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        ylim -- an integer value indicating the maximum value that will be presented
                OR -1 to set no limit
        outputDirectory -- a string indicating the directory where output files will be written
        binSize -- an integer value indicating the size of the bins that were used for counting
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
    '''
    numContigsPlotted = 0
    
    # Derive our output file name and error if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    fileOut = os.path.join(outputDirectory, f"onePlot.histo.{fileSuffix}")
    if os.path.isfile(fileOut):
        raise FileExistsError(f"'onePlot.histo.{fileSuffix}' already found in output directory")
    
    # Get the ordered contig IDs
    contigIDs = get_sorted_contig_ids(histoDict.keys())
    
    # Get each contigs' plot data
    plotData = []
    for contigID in contigIDs:
        x = np.arange(0, len(histoDict[contigID]))
        y = histoDict[contigID]
        
        # Enforce ylim if set
        if ylim != -1:
            y[y > ylim] = ylim
        
        plotData.append([contigID, x, y])
        numContigsPlotted += 1
    
    # Produce the figure axes
    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    gs = fig.add_gridspec(1, len(plotData), hspace=0)
    axes = gs.subplots(sharey='row')
    
    ## Set the figure title
    stepSize = binSize / 1000 # convert to Kbp
    fig.suptitle(f"{statisticLabel} histogram", fontweight="bold")
    fig.supxlabel(f"Bin number ({stepSize} Kbp step and width)", fontweight="bold")
    fig.supylabel(f"Number of read alignments with {statisticLabel}", fontweight="bold")
    
    # Plot the data into each axis
    with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
        # Write TSV header if applicable
        if createTSV:
            fileOutTSV.write("contigID\tposition\tdepth\n")
        
        for ax, (contigID, x, y) in zip(axes, plotData):
            # Set plot title
            ax.set_title(contigID)
            
            # Plot bars
            ax.bar(x, y, zorder=0)
            
            # Write TSV file if requested
            if createTSV:
                for xVal, yVal in zip(x*1000000, y): # convert back to bp
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\n")
        
        for ax in fig.get_axes():
            ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    
    return numContigsPlotted

def histo_regions(histoDict, regions, width, height, ylim, outputDirectory,
                  binSize, plotPDF=False,
                  statisticLabel="Q13 depth", createTSV=False):
    '''
    Parameters:
        histoDict -- dictionary linking chromosome IDs (keys) to numpy arrays of integers
                     containing the number of SNPs in each bin that exceeded the threshold
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        ylim -- an integer value indicating the maximum value that will be presented
                OR -1 to set no limit
        outputDirectory -- a string indicating the directory where output files will be written
        binSize -- an integer value indicating the size of the bins that were used for counting
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
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
        fileOut = os.path.join(outputDirectory, f"{contigID}.{start}_to_{end}.histo.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Histogram plot for '{contigID, start, end}' already found in output directory; skipping...")
            continue
        
        # Figure out which bins fall within this region
        binStart = start // binSize
        binEnd = end // binSize
        regionBins = histoDict[contigID][binStart:binEnd+1]
        
        # Raise error if incorrect number of bins found
        if len(regionBins) != binEnd - binStart + 1:
            raise ValueError(f"Region '{contigID, start, end}' only has {len(regionBins)} bins " +
                             f"but should have {binEnd - binStart + 1} bins. This probably means " +
                             "that your region coordinates are incorrect.")
        
        # Format values within this region
        x = np.arange(binStart, binEnd+1)
        y = regionBins
        
        # Enforce ylim if set
        if ylim != -1:
            y[y > ylim] = ylim
        
        # Configure plot
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        stepSize = binSize / 1000 # convert to Kbp
        ax.set_xlabel(f"Bin number ({stepSize} Kbp step and width)", fontweight="bold")
        ax.set_ylabel(f"Number of read alignments with {statisticLabel}", fontweight="bold")
        ax.set_title(f"{contigID} {statisticLabel} histogram", fontweight="bold")
        
        # Plot bars
        ax.bar(x, y, zorder=0)
        
        # Save output file
        plt.savefig(fileOut)
        numContigsPlotted += 1
        
        # Write TSV file if requested
        if createTSV:
            with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") as fileOutTSV:
                fileOutTSV.write("contigID\tposition\tdepth\n")
                for xVal, yVal in zip(x*1000000, y): # convert back to bp
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\n")
    
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
    usage = """%(prog)s receives a samtools depth output file and creates
    line plots and histograms of the alignment depth in various manners
    (per contig, horizontally, or for specific regions). Through calculating
    a weighted moving average (for line plots), it can help to visualise where
    in the genome regions that gaps in coverage are likely to be found.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-d", dest="depthFile",
                    required=True,
                    help="Specify the location of the input depth TSV file")
    p.add_argument("-o", dest="outputDirectory",
                    required=True,
                    help="Output directory where plot files will be written")
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
    p.add_argument("--ylim", dest="ylim",
                    type=int,
                    required=False,
                    help="""Optionally, specify the ylim value to set the maximum
                    value that will be presented; default is -1, which means that
                    NO limit will be set""",
                    default=-1)
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
    p.add_argument("--tsv", dest="createTSV",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want to have the plotted
                    data output as a TSV file""",
                    default=False)
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser()
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    lineparser = subparsers.add_parser("line",
                                       parents=[p],
                                       add_help=False,
                                       help="Create line plots of alignment depth")
    lineparser.set_defaults(func=linemain)
    
    histoparser = subparsers.add_parser("histogram",
                                        aliases=["histo"],
                                        parents=[p],
                                        add_help=False,
                                        help="Create histograms of alignment depth")
    histoparser.set_defaults(func=histomain)
    
    # Line-subparser arguments
    lineparser.add_argument("--wmaSize", dest="wmaSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the number of previous values to consider
                            during weighted moving average calculation (default=5)""",
                            default=5)
    lineparser.add_argument("--linewidth", dest="linewidth",
                            type=float,
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
                            help="""Optionally, specify the bin size to count reads
                            within (default=10000)""",
                            default=10000)
    
    args = subParentParser.parse_args()
    validate_args(args)
    
    # Figure out what our pickle file should be called
    pickleFile = os.path.join(
        args.outputDirectory,
        f"{os.path.basename(args.depthFile)}.pkl"
    )
    
    # Load pickle if it exists to skip computation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            dotsX, dotsY = pickle.load(fileIn)
    
    # Otherwise ...
    else:
        # Parse Euclidean distance data
        dotsX, dotsY = get_depth_for_dotting(args.depthFile)
        
        # Save data
        with open(pickleFile, "wb") as fileOut:
            pickle.dump([dotsX, dotsY], fileOut)
    
    # Figure out the lengths of each contig
    lengthsDict = {}
    for contigID in dotsX.keys():
        lengthsDict[contigID] = dotsX[contigID][-1]
    
    # Drop any contigs which don't meet our length cutoff
    for contigID, length in lengthsDict.items():
        if length < args.minimumContigSize:
            print(f"NOTE: '{contigID}' is below the minimum contig size and will be skipped")
            del dotsX[contigID]
            del dotsY[contigID]
    
    # Check that we still have contigs to plot
    if dotsX == {}:
        raise ValueError("ERROR: We didn't find any contigs which matched or exceeded the minimum size. " +
                         "Hence, no output files have been generated! Maybe you should fix your " +
                         "--minimum_contig value?")
    
    # Split into mode-specific functions
    if args.mode == "line":
        linemain(args, dotsX, dotsY)
    elif args.mode in ["histogram", "histo"]:
        histomain(args, dotsX, dotsY, lengthsDict)

def linemain(args, dotsX, dotsY):
    # Create plots
    if args.onePlot:
        numContigsPlotted = lineplot_horizontal(dotsX, dotsY, args.wmaSize, args.ylim,
                                                args.width, args.height,
                                                args.outputDirectory, args.plotPDF,
                                                args.showDots, args.linewidth,
                                                createTSV=args.createTSV)
    elif args.regions != []:
        numContigsPlotted = lineplot_regions(dotsX, dotsY, args.regions,
                                             args.wmaSize, args.ylim,
                                             args.width, args.height,
                                             args.outputDirectory, args.plotPDF,
                                             args.showDots, args.linewidth,
                                             createTSV=args.createTSV)
    else:
        numContigsPlotted = lineplot_per_contig(dotsX, dotsY, args.wmaSize, args.ylim,
                                                args.width, args.height,
                                                args.outputDirectory, args.plotPDF,
                                                args.showDots, args.linewidth,
                                                createTSV=args.createTSV)
    
    # Raise errors if necessary
    if numContigsPlotted == 0:
        raise ValueError("ERROR: We ended up skipping every contig! This means the program has " + 
                         "already run to completion previously. Hence, no new output files have been " +
                         "generated! Maybe you should delete the existing files to restart?")
    
    print("Program completed successfully!")

def histomain(args, dotsX, dotsY, lengthsDict):
    # Bin data into histograms
    histoDict = {}
    for contigID in get_sorted_contig_ids(dotsX.keys()):
        histoDict[contigID] = np.array([ 0
                for windowChunk in range(math.ceil(lengthsDict[contigID] / args.binSize))
            ])
        for x, y in zip(dotsX[contigID], dotsY[contigID]):
            binIndex = x // args.binSize
            histoDict[contigID][binIndex] += y
    
    # Create plots
    if args.onePlot:
        numContigsPlotted = histo_horizontal(histoDict,
                                             args.width, args.height, args.ylim,
                                             args.outputDirectory, args.binSize,
                                             args.plotPDF, createTSV=args.createTSV)
    elif args.regions != []:
        numContigsPlotted = histo_regions(histoDict, args.regions,
                                          args.width, args.height, args.ylim,
                                          args.outputDirectory, args.binSize,
                                          args.plotPDF, createTSV=args.createTSV)
    else:
        numContigsPlotted = histo_per_contig(histoDict,
                                             args.width, args.height, args.ylim,
                                             args.outputDirectory, args.binSize,
                                             args.plotPDF, createTSV=args.createTSV)
    
    # Raise errors if necessary
    if numContigsPlotted == 0:
        raise ValueError("ERROR: We ended up skipping every contig! This means the program has " + 
                         "already run to completion previously. Hence, no new output files have been " +
                         "generated! Maybe you should delete the existing files to restart?")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
