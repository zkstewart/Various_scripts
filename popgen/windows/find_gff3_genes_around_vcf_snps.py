#! python3
# find_gff3_genes_around_vcf_snps.py
# Script to take in a GFF3 file and a VCF to perform
# several steps and extract genes surrounding SNPs of
# interest (usually outlier SNPs). It's being built for
# the Chapa 2022 project according to Pete's advice.

import os, argparse, math

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3):
        print('I am unable to locate the GFF3 file (' + args.gff3 + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()
    # Validate numeric inputs
    if args.numNearest < 0:
        print('numNearest value should be greater than or equal to 0')

def get_gff3_overlap_dicts(gff3File):
    '''
    This function assumes the GFF3 only includes protein-coding genes that follow
    a canonical structure of gene -> mRNA -> exon/CDS. It should also be organised
    such that a gene line ALWAYS precedes its relevant mRNA lines, which ALWAYS
    precede any exon/CDS lines.
    '''
    chromosomeGeneRegions = {}
    geneCodingRegions = {}
    mrnaToGeneMap = {} # This will let us map back our mRNA IDs to their parent gene ID
    
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            l = line.rstrip("\r\n").split("\t")
            chrom = l[0]
            annotType = l[2]
            start = int(l[3])
            end = int(l[4])
            details = l[8].strip("\"").split(';')
            detail_dict = {}
            for i in range(len(details)):
                if details[i] == '':
                    continue
                split_details = details[i].split('=', maxsplit=1)
                detail_dict[split_details[0]] = split_details[1]
            
            # Set up storage structure at the start of each chromosome
            if chrom not in chromosomeGeneRegions:
                chromosomeGeneRegions[chrom] = []
            
            # Store gene details into their relevant dictionary structures
            if annotType == "gene":
                geneID = detail_dict["ID"]
                chromosomeGeneRegions[chrom].append([start, end, geneID]) # Store overall gene boundaries
                geneCodingRegions[geneID] = [] # Set up storage structure for exon/CDS regions
            elif annotType == "mRNA" or annotType == "transcript":
                mrnaID = detail_dict["ID"]
                geneID = detail_dict["Parent"]
                mrnaToGeneMap[mrnaID] = geneID # Allows us to track the mRNA back to its gene ID when handling exons and CDS
            elif annotType == "exon" or annotType == "CDS":
                mrnaID = detail_dict["Parent"]
                geneID = mrnaToGeneMap[mrnaID]
                geneCodingRegions[geneID].append([start, end, annotType])
            
    return chromosomeGeneRegions, geneCodingRegions

def get_snp_positions_from_vcf(vcfFile):
    snpPositions = {}
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            l = line.rstrip("\r\n").split("\t")
            chrom = l[0]
            pos = int(l[1])
            
            # Store in dictionary
            if chrom not in snpPositions:
                snpPositions[chrom] = []
            snpPositions[chrom].append(pos)
            
    return snpPositions

def find_genes_containing_snp(chromosomeGeneRegions, geneCodingRegions, snpPositions):
    genes = {} # we use a dictionary to make it easy to "overwrite" with more relevant annotations; also keeps non-redundant
    usedSNPs = {} # this is so we don't reuse the SNP for any other metrics
    for chrom, positions in snpPositions.items():
        for pos in positions:
            if chrom not in chromosomeGeneRegions:
                continue
            for gStart, gEnd, geneID in chromosomeGeneRegions[chrom]:
                if gStart <= pos and pos <= gEnd: # if it overlaps...
                    # Find out which part of the gene it overlaps
                    exon, cds = False, False
                    for rStart, rEnd, annotType in geneCodingRegions[geneID]:
                        if rStart <= pos and pos <= rEnd: # if it overlaps...
                            if annotType == "exon":
                                exon = True
                            else:
                                cds = True
                    usedSNPs.setdefault(chrom, [])
                    usedSNPs[chrom].append(pos)
                    # Derive the position of the SNP relative to the gene
                    if cds:
                        genes[geneID] = "CDS" # always overwrite with most important
                    elif exon:
                        if geneID in genes and genes[geneID] != "CDS": # overwrite only if it's not more important
                            genes[geneID] = "UTR"
                    else:
                        if geneID not in genes: # only write if there's nothing better
                            genes[geneID] = "intron"
    return genes, usedSNPs

def update_snpPositions_with_usedSNPs(snpPositions, usedSNPs):
    newSnpPositions = {}
    for chrom, positions in snpPositions.items():
        if chrom in usedSNPs:
            usedPositions = usedSNPs[chrom]
            freePositions = list(set(positions).symmetric_difference(usedPositions))
            if freePositions != []:
                newSnpPositions.setdefault(chrom, [])
                newSnpPositions[chrom] = freePositions
        else:
            newSnpPositions.setdefault(chrom, [])
            newSnpPositions[chrom] = positions
    return newSnpPositions

def find_genes_surrounding_snp(chromosomeGeneRegions, snpPositions):
    '''
    This function assumes SNPs have been filtered so that SNPs overlapping a gene are
    absent. It won't cause any bugs, but the notion of finding a gene to the left
    and right of a SNP that's inside of a gene loses meaning since we'll find the
    same gene in both directions.
    '''
    genes = []
    for chrom, positions in snpPositions.items():
        for pos in positions:
            if chrom not in chromosomeGeneRegions:
                continue
            
            nearestLeft = [math.inf, None]
            nearestRight = [math.inf, None]
            for gStart, gEnd, geneID in chromosomeGeneRegions[chrom]:
                isLeft = gEnd < pos
                isRight = gStart > pos
                if isLeft:
                    dist = pos - gEnd
                    if dist < nearestLeft[0]:
                        nearestLeft = [dist, geneID]
                elif isRight:
                    dist = gStart - pos
                    if dist < nearestRight[0]:
                        nearestRight = [dist, geneID]
                
            if nearestLeft[1] != None:
                genes.append(nearestLeft[1])
            if nearestRight[1] != None:
                genes.append(nearestRight[1])
            
    return list(set(genes))

def find_genes_nearest_snp(chromosomeGeneRegions, snpPositions, numNearest):
    genes = []
    for chrom, positions in snpPositions.items():
        for pos in positions:
            if chrom not in chromosomeGeneRegions:
                continue

            nNearest = []
            for gStart, gEnd, geneID in chromosomeGeneRegions[chrom]:
                isLeft = gEnd < pos
                isRight = gStart > pos
                if isLeft:
                    dist = pos - gEnd
                    nNearest.append([dist, geneID])
                elif isRight:
                    dist = gStart - pos
                    nNearest.append([dist, geneID])
            nNearest.sort(key = lambda x: x[0])
            
            for n in range(numNearest):
                if n < len(nNearest): # only append as many entries as exist if numNearest is large
                    genes.append(nNearest[n][1])
            
    return list(set(genes))

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3 and VCF files to obtain genes that are surrounding
    SNPs of interest e.g., outlier SNPs. The output file will be a list of genes including
    the reason for their identification i.e., which SNP they're near to, and what metric was
    used to identify it as a key gene of interest.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcf", required=True,
        help="Input VCF containing only SNPs of interest (filter it beforehand)")
    p.add_argument("-g", dest="gff3", required=True,
        help="Input GFF3 file")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output text file name")
    ## Optional
    p.add_argument("--N", dest="numNearest", required=False, type=int,
        help="Optionally, specify how many nearest genes should be found for a SNP (default==2)",
        default=2)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3 to get dictionarys that will let us locate gene overlaps / nearest genes
    chromosomeGeneRegions, geneCodingRegions = get_gff3_overlap_dicts(args.gff3)
    
    # Parse VCF as genotypes dictionary
    snpPositions = get_snp_positions_from_vcf(args.vcf)
    
    # Interesting gene metric 1: SNP overlapping a gene
    metric1GenesDict, usedSNPs = find_genes_containing_snp(chromosomeGeneRegions, geneCodingRegions, snpPositions)
    freePositions = update_snpPositions_with_usedSNPs(snpPositions, usedSNPs)
    
    # Interesting gene metric 2: genes to the left and right of a SNP
    metric2GenesList = find_genes_surrounding_snp(chromosomeGeneRegions, freePositions)
    
    # Interesting gene metric 3: nearest N genes to a SNP
    metric3GenesList = find_genes_nearest_snp(chromosomeGeneRegions, freePositions, args.numNearest)
    
    # Merge all metric results
    for geneID in metric2GenesList:
        if geneID not in metric1GenesDict:
            metric1GenesDict[geneID] = "near"
    for geneID in metric3GenesList:
        if geneID not in metric1GenesDict:
            metric1GenesDict[geneID] = "near"
    
    # Produce output file0
    with open(args.outputFileName, "w") as fileOut:
        for key, value in metric1GenesDict.items():
            fileOut.write("{0}\t{1}\n".format(key, value))
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
