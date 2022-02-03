#! python3
# pileup_to_snps.py
# Script to parse a samtools pileup produced by mpileup and
# identify SNPs using custom filtering rules

import os, argparse
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.pileupFile):
        print('I am unable to locate the samtools mpileup file (' + args.pileupFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()
    # Validate numeric inputs
    if args.floorCount <= 0:
        print("floorCount should be greater than 0")
        quit()
    if args.coverageCutoff < 0:
        print("coverageCutoff should be greater than or equal to 0")
        quit()
    if args.maf <= 0.00:
        print("maf should be greater than 0")
        quit()
    if args.maf >= 1.00:
        print("maf should be less than 1")
        quit()

## Data filtering
def mpileup_to_snpPiles(mpileupFile):
    snpPiles = []
    with open(mpileupFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            # Extract relevant information
            chrom = l[0]
            pos = l[1]
            ref = l[2]
            piles = [l[i:i+3] for i in range(3, len(l), 3)]
            # Store results
            snpPiles.append([chrom, pos, ref, piles])
    return snpPiles

def filter_snpPiles(snpPiles, floorCount, coverageCutoff):
    filteredPiles = []
    for pile in snpPiles:
        # Filter 1: Floor count
        if int(pile[0]) < floorCount:
            continue
        # Filter 2: Total coverage
        totalCoverage = sum([int(p[0]) for p in pile[3]]) # pile[3] is the list of triplets [coverage, description, ASCII score]
        if totalCoverage > coverageCutoff:
            continue
        # Store results if above filters pass
        filteredPiles.append(pile)
    return filteredPiles

## Data exploration
def plot_pile_statistics(snpPiles, boxplotName, histogramName):
    # Calculate position coverage
    covs = []
    for value in snpPiles:
        piles = value[3] # pile[3] is the list of triplets [coverage, description, ASCII score]
        totalCoverage = sum([int(pile[0]) for pile in piles])
        covs.append(totalCoverage)
    covs = np.array(covs)
    
    # Plot 1: Boxplot of position coverage
    fig = plt.figure(figsize = (5, 12))
    plt.boxplot(covs)
    #plt.show()
    plt.savefig(boxplotName)
    
    # Plot 2: Histogram of position coverage
    fig = plt.figure(figsize = (12, 6))
    plt.hist(x=covs, bins=50, color='#0504aa',
             alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    #plt.show()
    plt.savefig(histogramName)

## Main
def main():
    # User input
    usage = """%(prog)s reads in a samtools pileup (mpileup) output file and,
    according to rules identified from Haenel et al.'s 2021 paper, identifies
    SNPs.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="pileupFile", required=True,
        help="Input samtools mpileup output file")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the filtered SNPs")
    ## Optional
    p.add_argument("--floor", dest="floorCount", type=int,
        help="""This number is the minimum read depth required
        from all samples to be considered; default=0""",
        default=0)
    p.add_argument("--repeatCutoff", dest="coverageCutoff", type=int,
        help="""This number is the maximum cumulative read coverage
        allowed for a site before being discarded as a repeat region; default=0""",
        default=0)
    p.add_argument("--maf", dest="maf", type=float,
        help="""This number is the minor allele frequency (MAF) required
        for a SNP to be called; default=0.04""",
        default=0.04)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get SNP piles
    snpPiles = mpileup_to_snpPiles(args.pileupFile, args.floorCount)
    
    # Generate exploratory plots
    plot_pile_statistics(snpPiles, 'piles_boxplot.png', 'piles_histogram.png')
    
    # Filter SNP piles
    filteredPiles = filter_snpPiles(snpPiles, args.floorCount, args.coverageCutoff)
    
    # Convert SNP piles to output format
    ## TBD
    
    ## TESTING
    pickle.dump(snpPiles, open("snpPiles.pickle", "wb"))


if __name__ == "__main__":
    main()


## Testing zone
highCov = []
covs = []
for spile in snpPiles:
    piles = spile[2]
    totalCoverage = sum([int(pile[0]) for pile in piles])
    if totalCoverage > 300:
        highCov.append(spile)
        covs.append(totalCoverage)