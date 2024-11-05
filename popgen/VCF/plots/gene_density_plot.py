#! python3
# gene_density_plot.py
# Script to create visualisations of SNP density by gene
# across a genome, using a VCF file as input.

import os, argparse, math, gzip, sys
import matplotlib.pyplot as plt
from Bio import SeqIO
from contextlib import contextmanager

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_GFF3IO, ZS_VCFIO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.gff3File):
        print(f'I am unable to locate the input GFF3 file ({args.gff3File})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.vcfFile):
        print(f'I am unable to locate the input VCF file ({args.vcfFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        print(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Check for logical consistency of optional argument combination
    if args.filter != []:
        if not args.weightByOccurrence and not args.skipMonoallelic:
            raise ValueError("You have provided a list of samples to filter on, " + 
                             "but you have not enabled either the --weightByOccurrence " + 
                             "or --skipMonoallelic flags. Please enable at least one of these " + 
                             "flags to proceed, or omit --filter values.")
    
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
        print('This program will not allowing overwriting. If you want to re-do an analysis, stop now and fix this issue.')
    elif os.path.exists(args.outputDirectory):
        raise ValueError(f"The specified output directory ({args.outputDirectory}) is not a directory!")
    else:
        os.makedirs(args.outputDirectory, exist_ok=True)
        print(f"Output directory created at '{args.outputDirectory}' as part of argument validation.")

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def tally_variants_within_features(tallyDict, vcfFile, gff3, onlyCDS=False, weightByOccurrence=False,
                                   skipMonoallelic=False, filterSamples=[]):
    '''
    Modifies the tallyDict in place to count the number of SNPs within each gene feature.
    
    Parameters:
        tallyDict -- a dictionary with structure like:
                     {
                         "geneID1": [snpCount, geneLength],
                         "geneID2": [snpCount, geneLength],
                         ...
                     }
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        gff3 -- a ZS_GFF3IO.GFF3 object
        onlyCDS -- OPTIONAL; a boolean flag to indicate whether to only count
                   variants within CDS regions (default == False)
        weightByOccurrence -- OPTIONAL; a boolean flag to indicate whether to
                              count a variant each time it occurs in a sample
                              (default == False)
        skipMonoallelic -- OPTIONAL; a boolean flag to indicate whether to skip
                           counting a variant when ALL samples have a non-variant
                           allele (default == False)
        filterSamples -- OPTIONAL; a list of strings indicating the sample IDs
                         to filter on. Only relevant if you are using the
                         --weightByOccurrence flag or the --skipMonoallelic flag
                         (default == [])
    '''
    with ZS_VCFIO.open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # replace quotations may help with Excel files
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
                # for sample in filterSamples:
                #     if not sample in samples:
                #         raise ValueError(f"Sample indicated for filtration ({sample}) not found in VCF header!")
                filterIndices = set([ samples.index(sample) for sample in filterSamples ])
                continue
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            chrom = sl[0]
            pos = int(sl[1])
            sampleDetails = sl[9:]
            
            fieldsDescription = sl[8]
            if ":" not in fieldsDescription:
                gtIndex = 0
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Iterate through any gene features this variant overlaps
            geneFeatures = gff3.ncls_finder(pos, pos, "contig", chrom)
            for geneFeature in geneFeatures:
                geneID = geneFeature.ID
                
                # Narrow down to variants within the CDS region if onlyCDS
                if onlyCDS:
                    mrnaFeatures = [
                        childFeature
                        for childFeature in geneFeature.children
                        if childFeature.type == "mRNA"
                    ]
                    
                    cdsFeatures = [
                        childFeature
                        for mrnaFeature in mrnaFeatures
                        for childFeature in mrnaFeature.children
                        if childFeature.type == "CDS"
                    ]
                    
                    isWithinCDS = any([ cdsFeature.start <= pos <= cdsFeature.end for cdsFeature in cdsFeatures ])
                    if not isWithinCDS:
                        continue
                
                # Get sample genotypes if relevant
                if skipMonoallelic or weightByOccurrence:
                    genotypes = [
                        sampleDetails[i].split(":")[gtIndex].replace("|", "/").split("/")
                        for i in range(len(sampleDetails))
                        if i not in filterIndices
                    ]
                
                # Skip if we're skipping monoallelic and all samples are the same
                if skipMonoallelic:
                    genotypeSet = set(*genotypes).difference(".") # remove missing data
                    if len(genotypeSet) == 1:
                        continue
                
                # See how many times we should count this variant if relevant
                if weightByOccurrence:
                    snpCount = sum([ 1 for genotype in genotypes if set(genotype).difference(".") != "0" ])
                else:
                    snpCount = 1
                
                # Add to the tally
                tallyDict[geneID][0] += snpCount

def main():
    usage = """%(prog)s receives a VCF, a genome FASTA, and the GFF3 file for the genome.
    It then creates a SNP density plot for each contig in the genome, with the SNP density
    calculated per gene. The output is a series of PNG files, one per contig, showing the
    SNP density across the contig; unless you specify --onePlot, in which case a single plot
    will be created with all contigs positioned vertically.
    
    Note that the onlyCDS flag affects the geneLength value that is computed.
    If onlyCDS is True, the geneLength will be the sum of the lengths of all CDS
    regions within the gene. If onlyCDS is False, the geneLength will just be the
    length from the start to end of the entire region.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Specify the location of the input VCF file")
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Specify the location of the input GFF3 file")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where plot files will be written")
    # Opts
    p.add_argument("--onePlot", dest="onePlot",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want a single plot to be
                    produced with all chromosomes positioned horizontally""",
                    default=False)
    p.add_argument("--onlyCDS", dest="onlyCDS",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want to only count
                    variants within CDS regions""",
                    default=False)
    p.add_argument("--weightByOccurrence", dest="weightByOccurrence",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want to count a variant
                    each time it occurs in a sample""",
                    default=False)
    p.add_argument("--skipMonoallelic", dest="skipMonoallelic",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want to skip counting a
                    variant when ALL samples have a non-variant allele""",
                    default=False)
    p.add_argument("--filter", dest="filter",
                    required=False,
                    nargs="+",
                    help="""Optionally, indicate one or more samples to filter on. Only
                    relevant if you are using the --weightByOccurrence flag or the
                    --skipMonoallelic flag""",
                    default=[])
    
    args = p.parse_args()
    validate_args(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Parse GFF3 with NCLS indexing
    gff3 = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    gff3.create_ncls_index(typeToIndex="gene")
    
    # Establish tally dictionary
    tallyDict = {}
    for geneFeature in gff3.types["gene"]:
        geneID = geneFeature.ID
        if args.onlyCDS:
            longestMrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(geneFeature)
            geneLength = sum([ cdsFeature.end - cdsFeature.start + 1 for cdsFeature in longestMrnaFeature.CDS ])
        else:
            geneLength = geneFeature.end - geneFeature.start + 1
        tallyDict[geneID] = [0, geneLength]
    
    # Count SNPs within features
    tally_variants_within_features(tallyDict, args.vcfFile, gff3,
                                   args.onlyCDS, args.weightByOccurrence,
                                   args.skipMonoallelic, args.filter)
    
    ## TBD...
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
