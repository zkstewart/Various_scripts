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
    if not os.path.isfile(args.diffratioFile):
        print(f'I am unable to locate the difference ratio file ({args.diffratioFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        print(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.windowSize < 1:
        print("windowSize must be a positive integer")
        quit()
    if args.stepSize < 1:
        print("stepSize must be a positive integer")
        quit()
    if args.stepSize > args.windowSize:
        print("stepSize must be less than or equal to windowSize; otherwise, gaps will occur")
        quit()
    if args.minimumContigSize < (args.windowSize * 2):
        print("minimumContigSize must be at least 2x the window size")
        quit()
    if args.width < 0:
        print("width must be a positive integer")
        quit()
    if args.height < 0:
        print("height must be a positive integer")
        quit()
    if args.reportAboveCutoff != None:
        if 0 < args.reportAboveCutoff > 1:
            print("reportAboveCutoff must be a float value between 0 and 1")
            quit()
    # Check for conflicting arguments
    if args.onePlot and args.regions != []:
        print("You can't provide both --onePlot and --regions; please choose one and try again.")
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

def get_diffratio_density(diffratioFile, lengthsDict, windowSize=100000, stepSize=100000):
    '''
    Parameters:
        diffratioFile -- a string pointing to the difference ratio file
                         containing relevant statistics
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        windowSize -- OPTIONAL; an integer value indicating what size (diameter in bp) to average
                      difference ratio values over (default=100000)
        stepSize -- OPTIONAL; an integer value indicating the step size to move across the genome
                    (default=100000)
    '''
    HEADER_VALUES = ["CHROM", "POSI", "differenceRatio"]
    ratioDict = {}
    with open(diffratioFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if firstLine:
                assert all([ hv in sl for hv in HEADER_VALUES ]), "Header line doesn't contain expected values!"
                
                chromIndex = sl.index(HEADER_VALUES[0])
                posIndex = sl.index(HEADER_VALUES[1])
                ratioIndex = sl.index(HEADER_VALUES[2])
                
                firstLine = False
            
            # Handle content lines
            else:
                # Parse out relevant details from this line
                chrom, pos, ratio = sl[chromIndex], sl[posIndex], sl[ratioIndex]
                
                # Set up data storage bins for this chromosome if it doesn't exist
                if chrom not in ratioDict:
                    ratioDict[chrom] = {
                        "indices": np.array([ 0.0
                            for windowChunk in range(math.ceil(lengthsDict[chrom] / stepSize))
                        ]),
                        "counts": np.array([ 0
                            for windowChunk in range(math.ceil(lengthsDict[chrom] / stepSize))
                        ]),
                    }
                
                # Add difference ratio into any bins it belongs within
                "As determined by step size and window size; windowSize > stepSize will lead to bin overlap"
                
                ## Add to centre bin
                centreIndex = math.floor(int(pos) / stepSize)
                ratioDict[chrom]["indices"][centreIndex] += float(ratio)
                ratioDict[chrom]["counts"][centreIndex] += 1
                
                ## Look left for overlapping bins
                if windowSize > stepSize: # If the window size is larger than the step size, it may overlap with other bins
                    for leftIndex in range(centreIndex - 1, -1, -1):
                        # Calculate the window start and end positions for this step index
                        leftCentrePos = (leftIndex * stepSize)
                        leftStartPos = leftCentrePos - (windowSize / 2)
                        leftEndPos = leftCentrePos + (windowSize / 2)
                        
                        # If the SNP is within this window, add it to the bin
                        if int(pos) >= leftStartPos and int(pos) <= leftEndPos:
                            ratioDict[chrom]["indices"][leftIndex] += float(ratio)
                            ratioDict[chrom]["counts"][leftIndex] += 1
                        else:
                            break
                
                ## Look right
                if windowSize > stepSize:
                    for rightIndex in range(centreIndex + 1, len(ratioDict[chrom]["indices"])):
                        # Calculate the window start and end positions for this step index
                        rightCentrePos = (rightIndex * stepSize)
                        rightStartPos = rightCentrePos - (windowSize / 2)
                        rightEndPos = rightCentrePos + (windowSize / 2)
                        
                        # If the SNP is within this window, add it to the bin
                        if int(pos) >= rightStartPos and int(pos) <= rightEndPos:
                            ratioDict[chrom]["indices"][rightIndex] += float(ratio)
                            ratioDict[chrom]["counts"][rightIndex] += 1
                        else:
                            break
    
    # Average the difference ratio per window
    ## See https://stackoverflow.com/questions/26248654/how-to-return-0-with-divide-by-zero/37977222#37977222
    densityDict = {}
    for chrom in ratioDict.keys():
        densityDict.setdefault(chrom, [])
        a = ratioDict[chrom]["indices"]
        b = ratioDict[chrom]["counts"]
        average = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
        densityDict[chrom] = list(average)
    
    return densityDict

def plot_per_contig(lengthsDict, densityDict, minimumContigSize, windowSize, stepSize, width, height, outputDirectory, plotPDF=False):
    '''
    Parameters:
        lengthsDict -- a dictionary with structure like:
                       {
                            'contig1': intLength1,
                            'contig2': intLength2,
                            ...
                       }
        densityDict -- a dictionary with structure like:
                       {
                           'contig1': [ floatDensity1, floatDensity2, ... ],
                           'contig2': [ floatDensity1, floatDensity2, ... ],
                       }, where each list contains a number of float values equal to
                       the length of the contig divided by its step size
        minimumContigSize -- an integer value indicating the minimum size a contig must be
                             to be considered for plotting
        windowSize -- an integer value indicating what length of genome the difference ratio
                      values were averaged over
        stepSize -- an integer value indicating the step size used to move across the genome
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
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
            if not contigID in densityDict:
                print(f"WARNING: '{contigID}' is in the difference ratio file but has no SNPs associated " +
                      "with it; skipping...")
                continue
            
            # Get density values
            densityList = densityDict[contigID]
            
            # Smooth the curve for better visualisation
            "Smoothing the plot is probably a mistake, it's better to change the window size to emulate this effect"
            smoothedDensityList = densityList
            
            # Configure plot
            kbpWindowSize = round(windowSize / 1000, 2)
            kbpStepSize = round(stepSize / 1000, 2)
            
            fig = plt.figure(figsize=(width, height))
            ax = plt.axes()
            
            ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows; {kbpStepSize} kbp step size)", fontweight="bold")
            ax.set_ylabel("Mean difference ratio", fontweight="bold")
            ax.set_title(f"{contigID} difference ratio plot", fontweight="bold")
            
            ax.plot(smoothedDensityList)
            
            # Save output file
            plt.savefig(fileOut)
            numContigsPlotted += 1
    return numContigsProcessed, numContigsPlotted

def plot_once(lengthsDict, densityDict, minimumContigSize, windowSize, stepSize, width, height, outputDirectory, plotPDF=False):
    '''
    Parameters:
        lengthsDict -- a dictionary with structure like:
                       {
                            'contig1': intLength1,
                            'contig2': intLength2,
                            ...
                       }
        densityDict -- a dictionary with structure like:
                       {
                           'contig1': [ floatDensity1, floatDensity2, ... ],
                           'contig2': [ floatDensity1, floatDensity2, ... ],
                       }, where each list contains a number of float values equal to
                       the length of the contig divided by its step size
        minimumContigSize -- an integer value indicating the minimum size a contig must be
                             to be considered for plotting
        windowSize -- an integer value indicating what length of genome the difference ratio
                      values were averaged over
        stepSize -- an integer value indicating the step size used to move across the genome
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory where output files will be written
        plotPDF -- OPTIONAL; a boolean indicating whether to save output files as
                   PDFs (True) or PNGs (False)
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
            
            # Store density values
            plotData.append([densityDict[contigID], contigID])
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
    for ax, (densityList, contigID) in zip(axes, plotData):
        ax.plot(densityList)
        ax.set_title(contigID)
    for ax in fig.get_axes():
        ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    
    return numContigsProcessed, numContigsPlotted

def plot_regions(densityDict, regions, windowSize, stepSize, width, height, outputDirectory, plotPDF=False):
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
        
        # Configure plot
        kbpWindowSize = round(windowSize / 1000, 2)
        kbpStepSize = round(stepSize / 1000, 2)
        
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows; {kbpStepSize} kbp step size)", fontweight="bold")
        ax.set_ylabel("Mean difference ratio", fontweight="bold")
        ax.set_title(f"{contigID} difference ratio plot", fontweight="bold")
        
        ax.plot(xList, yList)
        
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
    p.add_argument("-d", dest="diffratioFile",
                   required=True,
                   help="Specify the location of the input VCF file")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where plot files will be written")
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
    # Opts (statistical behaviour)
    p.add_argument("--window_size", dest="windowSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the size of the window to sum
                   SNPs within (default=50000)""",
                   default=50000)
    p.add_argument("--step_size", dest="stepSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the size of the step to move across
                   a contig (default=50000)""",
                   default=50000)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the minimum size of contigs which
                   should have plots created for (default=200000)"; this value
                   must exceed --window_size by at least 2x""",
                   default=200000)
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
        f"{os.path.basename(args.diffratioFile)}.window{args.windowSize}.step{args.stepSize}.pkl"
    )
    
    # Load pickle if it exists
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            densityDict = pickle.load(fileIn)
    
    # Otherwise, tally SNPs over windows per contig
    else:
        densityDict = get_diffratio_density(args.diffratioFile, lengthsDict,
                                            args.windowSize, args.stepSize)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump(densityDict, fileOut)
    
    # Report any windows above the cutoff (if applicable)
    if args.reportAboveCutoff != None:
        print("# Difference ratio report:")
        for chrom, diffratioList in densityDict.items():
            for windowIndex, diffratio in enumerate(diffratioList):
                if diffratio >= args.reportAboveCutoff:
                    windowStart = (windowIndex * args.stepSize) - (args.windowSize / 2)
                    windowEnd = (windowIndex * args.stepSize) + (args.windowSize / 2)
                    print(f"# > Window index {windowIndex} from {windowStart} to {windowEnd} on '{chrom}' has a difference ratio of {diffratio}")
    
    # Create plots
    if args.onePlot:
        numContigsProcessed, numContigsPlotted = plot_once(lengthsDict, densityDict, args.minimumContigSize,
                                                           args.windowSize, args.stepSize,
                                                           args.width, args.height,
                                                           args.outputDirectory, args.plotPDF)
    elif args.regions != []:
        numContigsProcessed, numContigsPlotted = plot_regions(densityDict, args.regions,
                                                              args.windowSize, args.stepSize,
                                                              args.width, args.height,
                                                              args.outputDirectory, args.plotPDF)
    else:
        numContigsProcessed, numContigsPlotted = plot_per_contig(lengthsDict, densityDict, args.minimumContigSize,
                                                                 args.windowSize, args.stepSize,
                                                                 args.width, args.height,
                                                                 args.outputDirectory, args.plotPDF)
    
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
