#! python3
# pileup_to_snps.py
# Script to parse a samtools pileup produced by mpileup and
# identify SNPs using custom filtering rules

import os, argparse
from collections import Counter

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
    
    if args.mafCutoff <= 0.00:
        print("mafCutoff should be greater than 0")
        quit()
    if args.mafCutoff >= 1.00:
        print("mafCutoff should be less than 1")
        quit()
    
    if args.coverageCutoff == -1:
        pass # This is an accepted option
    elif args.coverageCutoff < 0:
        print("coverageCutoff should be greater than or equal to 0, or -1")
        quit()

## Data parsing and filtering
def mpileup_to_snpPiles(pileupFile, floorCount, coverageCutoff):
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
    
    Basic filtration occurs here to reduce the amount of data we need to store
    in memory.
    '''
    snpPiles = {}
    with open(pileupFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            # Extract relevant information
            chrom = l[0]
            pos = l[1]
            ref = l[2]
            coverages = [[int(l[i])] for i in range(3, len(l), 3)] # This ignores the type of alignment (match, mismatch) and ASCII scores
            # Perform basic filtration to reduce memory burden
            ## Filter 1: Floor count
            success = all([c[0] >= floorCount for c in coverages])
            if not success: continue
            ## Filter 2: Total coverage
            if coverageCutoff != -1:
                success = sum([c[0] for c in coverages]) <= coverageCutoff
                if not success: continue
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
    
    Lastly, the input snpPiles may have been filtered already, so we check to
    see if a position is present in snpPiles. If it's not, we just skip it.
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
            # Skip if we've filtered this position already
            if chrom not in snpPiles:
                continue
            elif pos not in snpPiles[chrom]:
                continue
            # Determine which field position we're extracting
            if ":" not in fieldsDescription:
                fieldIndex = -1
            else:
                fieldIndex = fieldsDescription.split(":").index("GT")
            # Parse genotype per sample
            ongoingCount = 0 # This gives us the index for the sample in order
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                if fieldIndex != -1:
                    gtField = sampleResult.split(":")[fieldIndex]
                else:
                    gtField = sampleResult
                genotype = gtField.replace("0", ref).replace("1", alt)
                # Store results
                snpPiles[chrom][pos][ongoingCount].append(genotype)
                ongoingCount += 1
    return snpPiles

def snpPiles_to_geno(snpPiles, mafCutoff):
    '''
    This assumes a snpPiles dictionary structure has been created through at least
    a two-step process of 1) mpileup_to_snpPiles, and 2) augment_snpPiles_with_GT_from_vcf.
    
    It also performs the MAF cutoff filtration
    '''
    geno = []
    for chrom in snpPiles.keys():
        for pos, value in snpPiles[chrom].items():
            # Filter 3: MAF cutoff
            alleles = [allele for v in value for allele in v[1].split("/")] # v[1] will be the genotype like "A/T" or "G/."
            alleleCounter = Counter(alleles)
            
            majorAlleleCount = alleleCounter.most_common(1)[0][1] # works like [('G', 3)][0][1]
            minorAlleleCount = len(alleles) - majorAlleleCount
            
            maf = minorAlleleCount / len(alleles)
            success = maf >= mafCutoff
            if not success: continue
            
            # Store results if above filters pass
            if geno == []: # add header line
                geno.append("#CHROM\tPOS\t{0}\n".format(
                    "\t".join(["ind" + str(i+1) for i in range(0, len(value))])
                ))
            geno.append("{0}\n".format(
                "\t".join([chrom, pos, "\t".join(
                    [v[1] for v in value]
                )])
            ))
    return geno

## File out
def write_geno_file(geno, outputFileName):
    with open(outputFileName, "w") as fileOut:
        fileOut.write("".join(geno))

## Data exploration
def plot_pile_statistics(snpPiles, boxplotName, histogramName):
    # Calculate position coverage
    covs = []
    for chrom in snpPiles.keys():
        for _, value in snpPiles[chrom].items():
            totalCoverage = sum([v[0] for v in value]) # value gives us a list with [[coverage, ...], ...]
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
    p.add_argument("-iv", dest="vcfFile", required=True,
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
        allowed for a site before being discarded as a repeat region;
        default=-1 which means no repeat filtering occurs""",
        default=-1)
    p.add_argument("--mafCutoff", dest="mafCutoff", type=float,
        help="""This number is the minor allele frequency (MAF) required
        for a SNP to be called; default=0.04""",
        default=0.04)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get SNP piles
    snpPiles = mpileup_to_snpPiles(args.pileupFile, args.floorCount, args.coverageCutoff)
    
    # Generate exploratory plots
    plot_pile_statistics(snpPiles, 'piles_boxplot.png', 'piles_histogram.png')
    
    # Augment SNP piles with genotype
    snpPiles = augment_snpPiles_with_GT_from_vcf(snpPiles, args.vcfFile)
    pickle.dump(snpPiles, open("snpPiles_augment.pickle", "wb")) ## TESTING
    
    # Filter SNP piles
    geno = snpPiles_to_geno(snpPiles, args.mafCutoff)
    
    # Write output .geno file
    write_geno_file(geno, args.outputFileName)

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