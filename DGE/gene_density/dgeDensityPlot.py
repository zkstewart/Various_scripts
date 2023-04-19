#! python3
# dgeDensityPlot.py
# Script to create visualisations of DE gene density
# for assessing hypotheses of DE gene distribution
# across chromosomes.

import os, argparse, math, sys
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.ndimage.filters import gaussian_filter1d

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_GFF3IO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.dgeTable):
        print(f'I am unable to locate the input DGE file ({args.dgeTable})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        print(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3File):
        print(f'I am unable to locate the input GFF3 file ({args.gff3File})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.rowIndex < 1:
        print("rowIndex must be a positive integer")
        quit()
    if args.columnIndex < 1:
        print("columnIndex must be a positive integer")
        quit()
    if args.windowSize < 1:
        print("windowSize must be a positive integer")
        quit()
    if args.minimumContigSize < (args.windowSize * 2):
        print("minimumContigSize must be at least 2x the window size")
        quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
        print('This program will not allowing overwriting. If you want to re-do an analysis, stop now and fix this issue.')

def get_deg_density(dgeTable, lengthsDict, gff3Obj, rowIndex, columnIndex, windowSize=100000):
    '''
    Parameters:
        dgeTable -- a string pointing to the file where DGE results have been tabulated
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        gff3Obj -- a ZS_GFF3IO.GFF3 object of the genome annotation for this species
        rowIndex -- a 1-based integer indicating which row in the table to start at
        columnIndex -- a 1-based integer which column in the table DGE results start at
        windowSize -- an integer value indicating what size to bin genes in
    '''
    densityDict = {}
    ongoingCount = 0
    with open(dgeTable, "r") as fileIn:
        for line in fileIn:
            ongoingCount += 1
            # Skip if we're not at the starting row
            if rowIndex > ongoingCount:
                continue
            
            # Parse out DE results
            sl = line.rstrip("\r\n ").split("\t")
            geneID = sl[0]
            deColumns = sl[columnIndex-1:]
            
            # Get gene details
            geneFeature = gff3Obj[geneID]
            contigID = geneFeature.contig
            start, end = geneFeature.coords
            
            # Figure out which bin this gene overlaps with
            "start -> end gives us the range of window chunks the gene model overlaps"
            startWindowIndex = math.floor(start / windowSize)
            endWindowIndex = math.floor(end / windowSize)
            
            # Store DE results
            densityDict.setdefault(contigID, {
                i : [ 0
                        for windowChunk in range(math.ceil(lengthsDict[contigID] / windowSize))
                    ]
                for i in range(len(deColumns))
            })
            
            for windowChunkIndex in range(startWindowIndex, endWindowIndex+1):
                for i in range(len(deColumns)):
                    if deColumns[i] != ".":
                        densityDict[contigID][i][windowChunkIndex] += 1
    
    return densityDict

def normalise_density(densityList, minValue=None, maxValue=None):
    '''
    Using min-max normalisation over a density list. Can compute minimum
    and maximum values in function if the graph is specific to just one line,
    or outside of function and feed that in through method params.
    
    Parameters:
        densityList -- a list containing integers of any length
        minValue -- None, OR an integer for the min value to normalise by
        maxValue -- None, OR an integer for the max value to normalise by
    '''
    if minValue == None or maxValue == None:
        minValue, maxValue = min(densityList), max(densityList)
    
    normalisedDensityList = [ (x - minValue) / (maxValue - minValue) for x in densityList ]
    return normalisedDensityList

def main():
    usage = """%(prog)s receives a DGE results table and creates density plots per
    chromosome representing the number of DE genes in each chromosomal window.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="dgeTable",
                   required=True,
                   help="Specify the location of the input DGE results table file")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Specify the location of the GFF3 annotations file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where plot files will be written")
    p.add_argument("-r", dest="rowIndex",
                   type=int,
                   required=True,
                   help="Specify the row where DGE results start from; use 1-based counting")
    p.add_argument("-c", dest="columnIndex",
                   type=int,
                   required=True,
                   help="""Specify the column where DGE results are indicated;
                   anything not "." will be considered a DE result; use 1-based counting""")
    # Opts
    p.add_argument("--window_size", dest="windowSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the size of the window to sum
                   DEGs within (default=100000)""",
                   default=100000)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the minimum size of contigs which
                   should have plots created for (default=1000000)"; this value
                   must exceed window_size by at least 2x""",
                   default=1000000)
    
    args = p.parse_args()
    validate_args(args)
    os.makedirs(args.outputDirectory, exist_ok=True)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }    
    
    # Load in GFF3
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    
    # Tally DEGs over windows per contig
    densityDict = get_deg_density(args.dgeTable, lengthsDict, gff3Obj, args.rowIndex, args.columnIndex, args.windowSize)
    
    # Create plots per contig
    numContigsProcessed = 0
    numContigsPlotted = 0
    for contigID, length in lengthsDict.items():
        if length >= args.minimumContigSize:
            numContigsProcessed += 1
            
            # Derive our output file name and skip if already existing
            fileOut = os.path.join(args.outputDirectory, f"{contigID}.png")
            if os.path.isfile(fileOut):
                print(f"WARNING: Plot for '{contigID}' already found in output directory; skipping...")
                continue
            
            # Skip if we found no DEGs on this contig
            if not contigID in densityDict:
                print(f"WARNING: '{contigID}' is in the genome FASTA but has no DEGs associated " +
                    "with it; skipping...")
                continue
            
            # Configure plot
            kbpWindowSize = round(args.windowSize / 1000, 2)
            
            fig, ax = plt.subplots(figsize=(10,6), nrows=1, ncols=1)
            ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows)", fontweight="bold")
            ax.set_ylabel("Min-max normalised DEG number per window", fontweight="bold")
            ax.set_title(f"{contigID} DEG density plot", fontweight="bold")
            
            # Loop through each DE comparison's density values and plot them
            maxValue = max([ v for i, values in densityDict[contigID].items() for v in values ])
            minValue = min([ v for i, values in densityDict[contigID].items() for v in values ])
            for i in densityDict[contigID].keys():
                densityList = densityDict[contigID][i]
                try:
                    normalisedDensityList = normalise_density(densityList, minValue, maxValue)
                except:
                    continue
                
                # Smooth the curve for better visualisation
                smoothedDensityList = gaussian_filter1d(normalisedDensityList, sigma=1)
                
                # Plot the line
                ax.plot(smoothedDensityList, label = f"DE column {i+1}")
            
            # Save output file
            fig.legend()
            fig.savefig(fileOut)
            plt.close(fig)
            
            numContigsPlotted += 1
        
    # Raise relevant warnings
    if numContigsProcessed == 0:
        print(f"WARNING: We didn't find any contigs which exceeded {args.minimumContigSize}bp in size")
        print("Hence, no output files have been generated! Maybe you should fix your --minimum_contig value?")
    elif numContigsPlotted == 0:
        print("WARNING: We ended up skipping every contig! This means the program has already run to completion previously.")
        print("Hence, no new output files have been generated! Maybe you should delete the existing folder to restart?")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
