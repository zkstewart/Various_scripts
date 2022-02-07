#! python3
# pileup_to_snps.py
# Script to parse a samtools pileup produced by mpileup and
# identify SNPs using custom filtering rules

import os, argparse, copy, re
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.pileupFile):
        print('I am unable to locate the pileup file (' + args.pileupFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the bcf file (' + args.vcfFile + ')')
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
    if args.mafCutoff <= 0.00:
        print("mafCutoff should be greater than 0")
        quit()
    if args.mafCutoff >= 1.00:
        print("mafCutoff should be less than 1")
        quit()

## Data filtering
def mpileup_to_snpPiles(pileupFile):
    '''
    This is going to produce a snpPiles data structure with format like:
    {chromosome:
        {pos:
            [coverage_1],
            ...
            [coverage_N]
        }
    }
    
    The sub-dictionary indexed by 'pos' will contain a list for each sample
    including the coverage at the chromosomal position. This list can be
    augmented by later functions.
    '''
    snpPiles = {}
    with open(pileupFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            # Extract relevant information
            chrom = l[0]
            pos = l[1]
            ref = l[2]
            coverages = [int(l[i]) for i in range(3, len(l), 3)] # This ignores the type of alignment (match, mismatch) and ASCII scores
            #piles = [l[i:i+3] for i in range(3, len(l), 3)]
            # Establish dictionary structure
            if chrom not in snpPiles:
                snpPiles[chrom] = {}
            # Store results
            snpPiles[chrom][pos] = coverages
    return snpPiles

def augment_snpPiles_with_GT_from_vcf(snpPiles, vcfFile):
    '''
    This assumes a dictionary-based snpPiles data structure is being
    received as from mpileup_to_snpPiles(). It adds the genotype information
    for each sample to the per-samples list which should already contain, at
    a minimum, the coverage at the position.
    
    It also assumes the order of the pileup file and VCF file are the same
    in terms of their sample listing.
    '''
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"): continue
            
            l = line.rstrip("\r\n").split("\t")
            # Extract relevant information
            chrom = l[0]
            pos = l[1]
            ref = l[3]
            alt = l[4]
            fieldsDescription = l[8]
            # Determine which field position we're extracting
            if ":" not in fieldsDescription:
                pos = -1
            else:
                pos = fieldsDescription.split(":").index("GT")
            # Parse genotype per sample
            ongoingCount = 0 # This gives us the index for the sample in order
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                if pos != -1:
                    gtField = sampleResult.split(":")[pos]
                else:
                    gtField = sampleResult
                genotype = gtField.replace("0", ref).replace("1", alt)
                # Store results
                snpPiles[chrom][pos][ongoingCount].append(genotype)
                ongoingCount += 1
    return snpPiles

def filter_snpPiles(snpPiles, floorCount, coverageCutoff, mafCutoff):
    '''
    TBD: Update to dictionary-based structure
    '''
    filteredPiles = []
    for pileGroup in snpPiles:
        # Filter 1: Floor count
        if int(pileGroup[0]) < floorCount:
            continue
        # Filter 2: Total coverage
        totalCoverage = sum([int(pile[0]) for pile in pileGroup[3]]) # pileGroup[3] is the list of triplets [coverage, description, ASCII score]
        if totalCoverage > coverageCutoff:
            continue
        # Filter 3: MAF cutoff
        alleleCount = []
        for pile in pileGroup[3]:
            # alleleFreqs = allele_frequency_from_subpile(pileGroup[2], pile) # pileGroup[2] is the reference base
            # diploidAlleles = call_diploid_alleles_from_alleleFreqs(alleleFreqs)
            pass
        # Store results if above filters pass
        filteredPiles.append(pile)
    return filteredPiles

## Data exploration
def plot_pile_statistics(snpPiles, boxplotName, histogramName):
    # Calculate position coverage
    covs = []
    for chrom in snpPiles.keys():
        for _, value in snpPiles[chrom].items():
            totalCoverage = sum([int(v[0]) for v in value]) # value gives us a list with [[coverage, ...], ...]
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
    p.add_argument("-ip", dest="pileupFile", required=True,
        help="Input samtools PILEUP format file")
    p.add_argument("-ib", dest="vcfFile", required=True,
        help="Input bcftools view VCF output file")
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
    p.add_argument("--mafCutoff", dest="mafCutoff", type=float,
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
    # filteredPiles = filter_snpPiles(snpPiles, args.floorCount, args.coverageCutoff)
    
    # Convert SNP piles to output format
    ## TBD
    
    ## TESTING
    pickle.dump(snpPiles, open("snpPiles.pickle", "wb"))


if __name__ == "__main__":
    main()


## Testing zone
# highCov = []
# covs = []
# for spile in snpPiles:
#     piles = spile[2]
#     totalCoverage = sum([int(pile[0]) for pile in piles])
#     if totalCoverage > 300:
#         highCov.append(spile)
#         covs.append(totalCoverage)

# bcfFile = r"F:\flies\chapa_2022\pileup\btrys06_mpileup.bcf"
# with open(bcfFile, "rb") as fileIn:
#     for line in fileIn:
#         stophere