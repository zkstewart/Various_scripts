#! python3
# snpDensityPlot.py
# Script to create visualisations of SNP density
# for assessing hypotheses of SNP distribution
# across chromosomes.

import os, argparse, math, gzip
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.ndimage.filters import gaussian_filter1d
from contextlib import contextmanager

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.vcfFile):
        print(f'I am unable to locate the input VCF file ({args.vcfFile})')
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
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
        print('This program will not allowing overwriting. If you want to re-do an analysis, stop now and fix this issue.')

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def get_vcf_density(vcfFile, lengthsDict, windowSize=100000):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        windowSize -- an integer value indicating what size to bin genes in
    '''
    densityDict = {}
    with open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            # Parse out contig lengths
            if line.startswith("#"):
                if "contig=" in line:
                    contigID = line.split("ID=")[1].split(",length=")[0].rstrip(">")
            # Handle content lines
            elif len(sl) >= 10:
                contigID, position = sl[0:2]
                densityDict.setdefault(contigID,
                    [ 0
                    for windowChunk in range(math.ceil(lengthsDict[contigID] / windowSize))
                    ]
                )
                windowChunkIndex = math.floor(int(position) / windowSize)
                densityDict[contigID][windowChunkIndex] += 1
    return densityDict

def normalise_density(densityList):
    '''
    Using min-max normalisation
    
    Parameters:
        densityList -- a list containing integers of any length
    '''
    minValue, maxValue = min(densityList), max(densityList)
    normalisedDensityList = [ (x - minValue) / (maxValue - minValue) for x in densityList ]
    return normalisedDensityList

def main():
    usage = """%(prog)s receives a VCF and creates SNP density plots per
    chromosome.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcfFile",
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
                   SNPs within (default=100000)""",
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
    
    # Tally SNPs over windows per contig
    densityDict = get_vcf_density(args.vcfFile, lengthsDict, args.windowSize)
    
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
                print(f"WARNING: '{contigID}' is in the VCF header but has no SNPs associated " +
                      "with it; skipping...")
                continue
            
            # Get density values
            densityList = densityDict[contigID]
            normalisedDensityList = normalise_density(densityList)
            
            # Smooth the curve for better visualisation
            smoothedDensityList = gaussian_filter1d(normalisedDensityList, sigma=1)
            
            # Configure plot
            fig = plt.figure(figsize=(10,6))
            ax = plt.axes()
            
            ax.set_xlabel("Chromosomal position (100 kbp windows)", fontweight="bold")
            ax.set_ylabel("Min-max normalised SNP number per window", fontweight="bold")
            ax.set_title(f"{contigID} SNP density plot", fontweight="bold")
            
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
