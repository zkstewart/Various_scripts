#! python3
# beagle_gene_haplotypes.py
# Script to take in a GFF3 file and a phased VCF
# as output by Beagle and predicts haplotypes for
# gene sequences (CDS only)

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find windows
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find GFF3IO
from windows import snp_proximity_report
from Function_packages import ZS_GFF3IO, ZS_SeqIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.beagleVCF):
        print('I am unable to locate the VCF file (' + args.beagleVCF + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3):
        print('I am unable to locate the GFF3 file (' + args.gff3 + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if 0 > args.minSnps:
        print("minSnps should be 0 or greater")
        quit()
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print('Specified output directory already exists and contains files.')
        print('This script would rather not mess with the possibility of overwriting stuff.')
        print('Make sure you specify a non-existing directory or delete the files in the directory and try again.')
        quit()

def get_phased_genotypes_from_vcf(vcfFile):
    '''
    This function expects the provided VCF file to contain phased genotypes
    as per Beagle's output. It doesn't necessarily have to be a VCF produced
    by Beagle, but it must have a GT field we can parse. Any ambiguous genotypes
    in any samples will result in the variant being completely ignored
    and a warning being printed.
    
    Returns:
        snpGenotypes -- a dictionary with structure like:
                        {
                            'contig1': {
                                "ref_alt": [refAllele, [altAllele1, ...]],
                                pos1: {
                                    'sample1': genotype,
                                    'sample2': genotype,
                                    ...
                                },
                                ...
                            },
                            'contig2': {
                                "ref_alt": [...],
                                pos1: { ... },
                                ...
                            }
                            ...,
                        }
    '''
    #snpPositions = {}
    snpGenotypes = {}
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = l[9:] # This gives us the ordered sample IDs
            if line.startswith("#"):
                continue
            
            # Determine which field position we're extracting to get our GT value
            fieldsDescription = l[8]
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Extract relevant details of the SNP
            chrom = l[0]
            pos = int(l[1])
            ref = l[3]
            alt = l[4].split(",")
            
            # Format a dictionary to store sample genotypes for this position
            posGenotypeDict = {}
            imputeWarningSkip = False
            ongoingCount = 0 # This gives us the index for our samples header list 
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                # Grab our genotype
                if gtIndex != -1:
                    genotype = sampleResult.split(":")[gtIndex]
                else:
                    genotype = sampleResult
                
                # Validate genotype is formatted correctly and not missing data
                if "." in genotype or "/" in genotype:
                    print(f"Warning: unphased/non-imputed genotype being skipped at {chrom}:{pos}")
                    imputeWarningSkip = True
                    break
                
                # Parse and store genotype
                samplePopulation = samples[ongoingCount]
                posGenotypeDict[samplePopulation] = list(map(int, genotype.split("|")))
                
                ongoingCount += 1
            if imputeWarningSkip is True:
                continue
            
            # Store in our overall dictionaries
            # snpPositions.setdefault(chrom, [])
            # snpPositions[chrom].append(pos)
            
            snpGenotypes.setdefault(chrom, {})
            snpGenotypes[chrom]["ref_alt"] = [ref, alt]
            snpGenotypes[chrom][pos] = posGenotypeDict
    
    #return snpPositions, snpGenotypes
    return snpGenotypes

def get_genotyped_snps_for_genes(gff3Obj, snpGenotypes):
    '''
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object with NCLS indexing of the gene features
        snpGenotypes -- a dictionary with structure like:
                        {
                            'contig1': {
                                "ref_alt": [refAllele, [altAllele1, ...]],
                                pos1: {
                                    'sample1': genotype,
                                    'sample2': genotype,
                                    ...
                                },
                                ...
                            },
                            'contig2': {
                                "ref_alt": [...],
                                pos1: { ... },
                                ...
                            }
                            ...,
                        }
    Returns:
        geneSnpDict -- a dictionary with structure like:
                       {
                           'geneID1': {
                               "ref_alt": [refAllele, [altAllele1, ...]],
                               pos1: {
                                   'sample1': genotype,
                                   'sample2': genotype,
                                   ...
                               }
                           },
                           ...
                       }
    '''
    geneSnpDict = {}
    
    for contig, positionsDict in snpGenotypes.items():
        for pos, genotypesDict in positionsDict.items():
            matches = gff3Obj.ncls_finder(pos, pos, "contig", contig)
            # If the SNP is located within a gene
            if matches != []:
                # Get the longest/representative mRNA for each gene
                mrnaFeatures = [ZS_GFF3IO.GFF3.longest_isoform(m) for m in matches]
                
                # For each mRNA, locate the SNP in the gene
                for feature in mrnaFeatures:
                    snpLocation, geneID = snp_proximity_report.locate_snp_within_gene([feature], pos) # hacky use of existing func
                    
                    # If we found the SNP in the CDS, index it
                    if snpLocation == "CDS":
                        geneSnpDict.setdefault(geneID, {})
                        geneSnpDict[geneID][pos] = genotypesDict
    
    return geneSnpDict

def get_haplotype_sequences(gff3Obj, genomeFASTA_obj, geneSnpDict, minSnps=0):
    '''
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object containing all the genes
                   referred to in geneSnpDict
        genomeFASTA_obj -- a ZS_SeqIO.FASTA object containing the contigs
                   that are referred to in our gff3Obj
        geneSnpDict -- a dictionary with structure like:
                       {
                           'geneID1': {
                               "ref_alt": [refAllele, [altAllele1, ...]],
                               pos1: {
                                   'sample1': genotype,
                                   'sample2': genotype,
                                   ...
                               }
                           },
                           ...
                       }
        minSnps -- an integer indicating how many SNPs must be
                   present in a gene for us to make a haplotype of it.
                   Genes without this many SNPs will not have any
                   sequences reported.
    Returns:
        haploSeqDict -- TBD
    '''
    assert isinstance(minSnps, int), \
        "minSnps must be an integer value"
    
    haploSeqDict = {}
    for geneID, snpDict in geneSnpDict.items():
        # Skip genes which fail minSnps threshold
        if len(snpDict) < minSnps:
            continue
        
        # Get the representative mRNA feature and sequence
        mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(gff3Obj[geneID])
        cds_FastASeq_obj, cds_featureType, cds_startingFrame = gff3Obj.retrieve_sequence_from_FASTA(genomeFASTA_obj, mrnaFeature.ID, "CDS")
        
        ## TBD: Continue from here

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3 and a phased VCF file as output by Beagle
    and generates haplotype predictions for genes containing a variant within its
    CDS. Outputs will be written to the provided directory.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="beagleVCF",
                   required=True,
                   help="Input phased VCF")
    p.add_argument("-g", dest="gff3",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-f", dest="fasta",
                   required=True,
                   help="Input genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    ## Optional
    p.add_argument("--minSnps", dest="minSnps",
                   type=int,
                   required=False,
                   help="""Optionally, specify how many SNPs are needed
                   to be considered a haplotype (default==2)""",
                   default=2)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Parse VCF to get SNP genotypes
    snpGenotypes = get_phased_genotypes_from_vcf(args.beagleVCF)
    
    # Parse genome FASTA
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.fasta)
    
    # Get the SNPs that are present within gene model CDS regions
    geneSnpDict = get_genotyped_snps_for_genes(gff3Obj, snpGenotypes)
    
    # Get haplotypes for each gene and sample
    haploSeqDict = get_haplotype_sequences(gff3Obj, genomeFASTA_obj, geneSnpDict)
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
