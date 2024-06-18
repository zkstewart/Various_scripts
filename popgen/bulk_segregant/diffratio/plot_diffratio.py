#! python3
# plot_diffratio.py
# Script to create visualisations of the difference ratio statistics
# for assessing hypotheses of variant segregation along chromosomes.

import os, argparse, math
import matplotlib.pyplot as plt
from Bio import SeqIO
#from scipy.ndimage import gaussian_filter1d

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
    if args.minimumContigSize < (args.windowSize * 2):
        print("minimumContigSize must be at least 2x the window size")
        quit()
    # if args.smoothingSigma < 0:
    #     print("smoothingSigma must be a float greater than or equal to zero")
    #     quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def get_diffratio_density(diffratioFile, lengthsDict, windowSize=100000):
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
        windowSize -- an integer value indicating what size to bin genes in
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
                
                # Add difference ratios into window-sized bins
                ratioDict.setdefault(chrom, {
                    "indices": [ 0
                        for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                    ],
                    "counts": [ 0
                        for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                    ],
                    }
                )
                
                windowChunkIndex = math.floor(int(pos) / windowSize)
                ratioDict[chrom]["indices"][windowChunkIndex] += float(ratio)
                ratioDict[chrom]["counts"][windowChunkIndex] += 1
    
    # Average the difference ratio per window
    densityDict = {}
    for chrom in ratioDict.keys():
        densityDict.setdefault(chrom, [])
        for windowChunkIndex in range(len(ratioDict[chrom]["indices"])):
            try:
                average = ratioDict[chrom]["indices"][windowChunkIndex] / ratioDict[chrom]["counts"][windowChunkIndex]
            except:
                average = 0.0
            densityDict[chrom].append(average)
      
    return densityDict

def main():
    usage = """%(prog)s receives a difference ratio TSV file and creates
    difference ratio plots per chromosome. It specifically calculates the
    average difference ratio value over each window and plots that value.
    Hence, it may help to visualise where in the genome regions associated
    with bulks exist.
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
    # Opts
    p.add_argument("--window_size", dest="windowSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the size of the window to sum
                   SNPs within (default=50000)""",
                   default=50000)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the minimum size of contigs which
                   should have plots created for (default=200000)"; this value
                   must exceed --window_size by at least 2x""",
                   default=200000)
    # p.add_argument("--smoothing", dest="smoothingSigma",
    #                type=float,
    #                required=False,
    #                help="""Optionally, specify how much smoothing should be applied
    #                to the plot; should be float value >= 0.0 (default=0.0)""",
    #                default=0.0)
    p.add_argument("--pdf", dest="plotPDF",
                   required=False,
                   action="store_true",
                   help="""Optionally, provide this flag if you want outputs to be
                   in PDF format instead of PNG format""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Tally SNPs over windows per contig
    densityDict = get_diffratio_density(args.diffratioFile, lengthsDict, args.windowSize)
    
    # Create plot per contig
    numContigsProcessed = 0
    numContigsPlotted = 0
    for contigID, length in lengthsDict.items():
        if length >= args.minimumContigSize:
            numContigsProcessed += 1
            
            # Derive our output file name and skip if already existing
            fileSuffix = "pdf" if args.plotPDF else "png"
            fileOut = os.path.join(args.outputDirectory, f"{contigID}.{fileSuffix}")
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
            kbpWindowSize = round(args.windowSize / 1000, 2)
            fig = plt.figure(figsize=(10,6))
            ax = plt.axes()
            
            ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows)", fontweight="bold")
            ax.set_ylabel("Mean difference ratio per window", fontweight="bold")
            ax.set_title(f"{contigID} difference ratio plot", fontweight="bold")
            
            ax.plot(smoothedDensityList)
            
            # Save output file
            plt.savefig(fileOut)
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
