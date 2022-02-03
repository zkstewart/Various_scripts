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
    if args.mafCutoff <= 0.00:
        print("mafCutoff should be greater than 0")
        quit()
    if args.mafCutoff >= 1.00:
        print("mafCutoff should be less than 1")
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
            snpPiles.append([chrom, pos, ref, piles]) # This forms a pileGroup
    return snpPiles

def filter_snpPiles(snpPiles, floorCount, coverageCutoff, mafCutoff):
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
            alleleFreqs = allele_frequency_from_subpile(pileGroup[2], pile) # pileGroup[2] is the reference base
            diploidAlleles = call_diploid_alleles_from_alleleFreqs(alleleFreqs)
            
        # Store results if above filters pass
        filteredPiles.append(pile)
    return filteredPiles

def allele_frequency_from_subpile(ref, pile):
    '''
    This function assumes it is receiving a triplet-formatted pile, where
    the first value is the read count, the second value are the samtools read
    bases encoding, and the third is the ASCII quality score for each read.
    Refer to http://www.htslib.org/doc/samtools-mpileup.html for details on the
    samtools pileup format.
    
    It also needs to know the reference base. Provide that please.
    '''
    newpile = copy.deepcopy(pile) # prevent changes to the original pile
    bases = newpile[1]
    
    # Remove symbols irrelevant to our calculations
    bases = re.sub("\^.", "", bases) # ^ is followed by the ASCII quality of the read alignment; we don't need it
    bases = bases.replace(">", "").replace("<", "") # I don't know how to handle reference skips
    
    # Get the allele for each position
    ## Handle insertions and deletions AFTER the read base
    insertAfter = re.findall(r"\+[0-9]+[ACGTNacgtn*#]+", bases) # count...
    bases = re.sub(r"\+[0-9]+[ACGTNacgtn*#]+", "", bases) # then remove them
    deleteAfter = re.findall(r"-[0-9]+[ACGTNacgtn]+", bases) # count...
    bases = re.sub(r"-[0-9]+[ACGTNacgtn]+", "", bases) # then remove them
    
    ## Handle simple matches, mismatches, and deletions
    matches = re.findall(r"[\.,]", bases)
    mismatches = re.findall(r"[ACGTNacgtn]", bases)
    deletions = re.findall(r"[\*#]", bases)
    
    # Tally the alleles
    alleles = {}
    for entry in [insertAfter, deleteAfter]:
        for value in entry:
            if value not in alleles:
                alleles[value] = 1
            else:
                alleles[value] += 1
                
    alleles[ref] = len(matches)
    
    for mm in mismatches:
        mm = mm.upper()
        if mm not in alleles:
            alleles[mm] = 1
        else:
            alleles[mm] += 1
    
    if len(deletions) > 0:
        alleles["deletion"] = len(deletions)
    
    # Calculate frequencies
    total = sum(alleles.values())
    alleleFreqs = {}
    for key, value in alleles.items():
        freq = value / total
        alleleFreqs[key] = freq
    
    return alleleFreqs

def call_diploid_alleles_from_alleleFreqs(alleleFreqs):
    '''
    alleleFreqs is assumed to be a dictionary with key:value pairs
    corresponding to the putative allele, and its frequency
    as a float ratio which should sum to approximately 1.0.
    
    A simple assumption is made in this function that all samples are
    diploid. This means we're going to see a monoallelic or biallelic
    position. If it's biallelic, our ratio should be close to 0.5 : 0.5
    for the two variants, otherwise we should see something close to
    1.0 : 0.0. Deviations outside of this range will be assumed to be noise.
    '''
    MONOALLELIC_CUTOFF = 0.7
    BIALLELIC_CUTOFF = 0.3
    
    # Get an ordered list of frequencies from alleleFreqs
    orderedFreqs = []
    for key, value in alleleFreqs.items():
        orderedFreqs.append([key, value])
    orderedFreqs.sort(key = lambda x: -x[1])
    
    # Handle scenario 1: monoallelic
    if orderedFreqs[0][1] >= MONOALLELIC_CUTOFF:
        return [orderedFreqs[0][0], orderedFreqs[0][0]]
    # Handle scenario 2: biallelic
    if len(orderedFreqs) > 1 and orderedFreqs[1][1] >= BIALLELIC_CUTOFF:
        return [orderedFreqs[0][0], orderedFreqs[1][0]]
    # Handle uncertainty (default to null hypothesis - there is no SNP here.)
    else:
        return [orderedFreqs[0][0], orderedFreqs[0][0]]

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