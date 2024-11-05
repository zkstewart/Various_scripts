#! python3
# calculate_persample_ed.py
# Script to perform a process akin to QTL-Seq, but for per-sample SNP data.
# It does not compute SNP-indices but instead borrows from Hill et. al. 2013
# and calculates the Euclidean distance between the two bulks based on allele
# frequencies.

# Load normal/pip packages
import os, argparse, sys
from math import sqrt

# Load functions from persample script
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find diffratio
from diffratio.calculate_persample_diffratio import parse_pops_metadata

# Load functions from ZS functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 4 dirs up is where we find Function_packages
from Function_packages import ZS_VCFIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print(f'I am unable to locate the VCF file ({args.vcfFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print(f'I am unable to locate the metadata file ({args.metadataFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def calculate_snp_ed(b1Gt, b2Gt):
    '''
    Parameters:
        b1Gt / b2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                       [
                           [0, 1],
                           [0, 0],
                           [1, 1],
                           ...
                       ]
    '''
    # Get all the unique alleles
    alleles = list(set([ allele for gt in b1Gt + b2Gt for allele in gt ]))
    if 0 not in alleles:
        alleles.append(0)
    alleles.sort()
    
    # Tally for bulk 1
    b1Count = { allele: 0 for allele in alleles }
    for allele1, allele2 in b1Gt:
        b1Count[allele1] += 1
        b1Count[allele2] += 1
    
    # Tally for bulk 2
    b2Count = { allele: 0 for allele in alleles }
    for allele1, allele2 in b2Gt:
        b2Count[allele1] += 1
        b2Count[allele2] += 1
    
    # Sum the number of genotyped alleles for each bulk
    b1Sum = sum(b1Count.values()) # bulk1_alleles
    b2Sum = sum(b2Count.values()) # bulk2_alleles
    
    # Calculate the Euclidean distance between the two bulks if possible
    if b1Sum == 0 and b2Sum == 0:
        return b1Sum, b2Sum, "." # euclidean distance is null
    elif b1Sum == 0 or b2Sum == 0:
        return b1Sum, b2Sum, 1 # euclidean distance is 1
    else:
        # Derive our euclidean distance value
        """Refer to "Euclidean distance calculation" in Hill et al. 2013"""
        edist = sqrt(sum([
            ((b1Count[allele] / b1Sum) - (b2Count[allele] / b2Sum))**2
            for allele in alleles
        ]))
        
        # Return the values
        return b1Sum, b2Sum, edist

def main():
    # User input
    usage = """%(prog)s receives a VCF produced from per-sample variant calling
    and a TSV metadata file containing the sample names and what segregating bulk
    they belong to; the TSV should contain no header, and have two columns, the first
    indicating the sample ID and the second indicating whether it's part of "bulk1"
    or "bulk2". It produces an output file similar to what's seen with QTL-Seq, but
    with the difference being that Euclidean distance is calculated between bulks
    rather than SNP-index.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-vcf", dest="vcfFile",
                   required=True,
                   help="Specify the VCF file name")
    p.add_argument("-meta", dest="metadataFile",
                   required=True,
                   help="Specify the metadata file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    # Opts
    p.add_argument("--ignoreIdentical", dest="ignoreIdentical",
                   required=False,
                   action="store_true",
                   help="""Optionally, provide this flag if you'd like variants where
                   both bulks are identical to be ignored; this can occur when both bulks
                   have the same variant with respect to the reference genome""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse metadata file
    metadataDict = parse_pops_metadata(args.metadataFile)
    metadataSamples = set([
        sample
        for samples in metadataDict.values()
        for sample in samples
    ])
    
    # Iteratively build output table
    with open(args.outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "CHROM", "POSI", "variant",
            "bulk1_alleles", "bulk2_alleles", 
            "euclideanDist"
        ])))
        
        # Parse through VCF file
        vcfIterator = ZS_VCFIO.SimpleGenotypeIterator(args.vcfFile)
        firstLine = True
        for values in vcfIterator:
            # Grab header line containing sample IDs
            if firstLine:
                vcfSamples = set(values)
                firstLine = False
                
                # Check that the metadata file matches the VCF file
                if vcfSamples != metadataSamples:
                    vcfDiff = vcfSamples.difference(metadataSamples)
                    metadataDiff = metadataSamples.difference(vcfSamples)
                    
                    print("ERROR: Metadata file does not match the VCF file!")
                    
                    if len(vcfDiff) > 0:
                        print("In your VCF, the following samples exist which are " + 
                            "absent from the metadata: ", ", ".join(vcfDiff))
                    if len(metadataDiff) > 0:
                        print("In your metadata, the following samples exist which are " + 
                            "absent from the VCF: ", ", ".join(metadataDiff))
                    quit()
            # Handle content lines
            else:
                contig, pos, ref, alt, snpDict = values
                ref_alt = [ref, *alt]
                
                # Figure out what type of variant this is
                if any([ x == "." for x in ref_alt ]):
                    variant = "indel"
                elif any([ len(ref_alt[0]) != len(ref_alt[x]) for x in range(1, len(ref_alt))]):
                    variant = "indel"
                else:
                    variant = "snp"
                
                # Split sample genotypes into bulk1 and bulk2
                bulk1 = [ snpDict[sample] for sample in metadataDict["bulk1"] if sample in snpDict ]
                bulk2 = [ snpDict[sample] for sample in metadataDict["bulk2"] if sample in snpDict ]
                
                # Calculate difference ratio
                bulk1_alleles, bulk2_alleles, \
                    euclideanDist = calculate_snp_ed(bulk1, bulk2)
                
                # Skip if both bulks are identical
                if args.ignoreIdentical and euclideanDist == 0:
                    bulk1Dedup = set(( tuple(x) for x in bulk1 ))
                    bulk2Dedup = set(( tuple(x) for x in bulk2 ))
                    "if both have set len==1, are the same, and have 0 Euclidean distance, they are identical non-reference alleles"
                    if len(bulk1Dedup) == 1 and bulk1Dedup == bulk2Dedup:
                        continue
                
                # Format and write output line
                outputLine = f"{contig}\t{pos}\t{variant}\t{bulk1_alleles}\t" + \
                    f"{bulk2_alleles}\t{euclideanDist}\n"
                
                fileOut.write(outputLine)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
