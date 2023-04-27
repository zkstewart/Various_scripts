#! python3
# bsa_snpIndexDensityPlot.py
# Script to create visualisations of SNP-index statistics
# generated by QTLseq for assessing hypotheses that
# segregating SNPs congregate in genomic locations

import os, argparse, math
import matplotlib.pyplot as plt
from Bio import SeqIO
#from scipy.ndimage.filters import gaussian_filter1d

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.indexFile):
        print(f'I am unable to locate the input SNP-index file ({args.indexFile})')
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
    if args.smoothingSigma < 0:
        print("smoothingSigma must be a float greater than or equal to zero")
        quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
        print('This program will not allowing overwriting. If you want to re-do an analysis, stop now and fix this issue.')

def get_snpindex_density(indexFile, lengthsDict, windowSize=100000):
    '''
    Parameters:
        vcfFile -- a string pointing to the SNP-index file created by QTLseq;
                   it should probably be the normal file i.e., not the p95 or
                   p99 files
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        windowSize -- an integer value indicating what size to bin genes in
    '''
    indexDict = {}
    with open(indexFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines (if relevant)
            if line.startswith("CHROM"):
                continue
            
            # Handle content lines
            else:
                # Parse out relevant details from this line
                chrom, pos, delta = *sl[0:2], float(sl[9])
                
                # Calculate the SNP-index-like value for this SNP
                windowChunkIndex = math.floor(int(pos) / windowSize)
                
                # Add SNP-index
                indexDict.setdefault(chrom, {
                    "indices": [ 0
                        for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                    ],
                    "counts": [ 0
                        for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                    ],
                    }
                )
                indexDict[chrom]["indices"][windowChunkIndex] += delta # may want to use abs() of this...?
                indexDict[chrom]["counts"][windowChunkIndex] += 1
    
    # Average the SNP-index per window
    densityDict = {}
    for chrom in indexDict.keys():
        densityDict.setdefault(chrom, [])
        for windowChunkIndex in range(len(indexDict[chrom]["indices"])):
            try:
                average = indexDict[chrom]["indices"][windowChunkIndex] / indexDict[chrom]["counts"][windowChunkIndex]
            except:
                average = 0.0
            densityDict[chrom].append(average)
      
    return densityDict

def main():
    usage = """%(prog)s receives a SNP-index file generated by QTLseq and
    creates SNP-index density plots per chromosome. It averages the SNP-index
    per window and plots that value. Hence, it may help to visualise
    where in the genome regions associated with bulks exist.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="indexFile",
                   required=True,
                   help="Specify the location of the input SNP-index file")
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
    p.add_argument("--smoothing", dest="smoothingSigma",
                   type=float,
                   required=False,
                   help="""TURNED OFF: Optionally, specify how much smoothing should be applied
                   to the plot; should be float value >= 0.0 (default=1.0)""",
                   default=1.0)
    
    args = p.parse_args()
    validate_args(args)
    os.makedirs(args.outputDirectory, exist_ok=True)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Tally SNPs over windows per contig
    densityDict = get_snpindex_density(args.indexFile, lengthsDict, args.windowSize)
    
    # Create plot per contig
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
            
            # Skip if we found no SNPs on this contig
            if not contigID in densityDict:
                continue # warning is unnecessary here
            
            # Get density values
            densityList = densityDict[contigID]
            
            # Smooth the curve for better visualisation
            "Smoothing is a mistake for this kind of analysis"
            # smoothedDensityList = gaussian_filter1d(densityList, sigma=args.smoothingSigma)
            smoothedDensityList = densityList
            
            # Configure plot
            kbpWindowSize = round(args.windowSize / 1000, 2)
            fig = plt.figure(figsize=(10,6))
            ax = plt.axes()
            
            ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows)", fontweight="bold")
            ax.set_ylabel("Mean SNP-index-like value per window", fontweight="bold")
            ax.set_title(f"{contigID} SNP-index plot", fontweight="bold")
            
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
