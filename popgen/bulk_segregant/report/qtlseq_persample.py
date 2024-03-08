#! python3
# qtlseq_persample.py
# Script to perform a process akin to QTL-Seq, but for per-sample SNP data.
# It does not compute SNP-indices in the true sense, but (what I believe is)
# a more relevant metric when looking at per-sample calling.

# Load normal/pip packages
import os, argparse, sys

# Load functions from other scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find windows
import haplotypes

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

def parse_pops_metadata(metadataFile):
    # Parse file
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            if "\t" in line:
                sl = line.rstrip("\r\n ").split("\t")
            elif "," in line:
                sl = line.rstrip("\r\n ").split(",")
            else:
                print("ERROR: Metadata file should be tab-delimited or comma-delimited!")
                quit()
            
            if sl == []:
                continue
            else:
                sample, pop = sl
                metadataDict.setdefault(pop, set())
                metadataDict[pop].add(sample)
    
    # Make sure it's valid for this script
    if len(metadataDict) == 2:
        assert all([ pop in ['bulk1', 'bulk2'] for pop in metadataDict.keys() ]), \
            "ERROR: Metadata file should have 2 or 3 populations: bulk1, bulk2, and (optionally) parent!"
    elif len(metadataDict) == 3:
        assert all([ pop in ['bulk1', 'bulk2', 'parent'] for pop in metadataDict.keys() ]), \
            "ERROR: Metadata file should have 2 or 3 populations: bulk1, bulk2, and (optionally) parent!"
    else:
        print("ERROR: Metadata file should have 2 or 3 populations: bulk1, bulk2, and (optionally) parent!")
        quit()
    
    return metadataDict

def calculate_snp_indices(b1Gt, b2Gt):
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
    
    # Calculate reference index for each bulk if possible
    if b1Sum == 0:
        b1RefIndex = "."
    else:
        b1RefIndex = b1Count[0] / b1Sum
    
    if b2Sum == 0:
        b2RefIndex = "."
    else:
        b2RefIndex = b1Count[0] / b2Sum
    
    # Calculate the delta reference index if possible
    if b1RefIndex == "." or b2RefIndex == ".":
        deltaRefIndex = "."
    else:
        deltaRefIndex = abs(b1RefIndex - b2RefIndex)
    
    # Calculate the difference ratio if possible
    if b1Sum == 0 and b2Sum == 0:
        return b1Sum, b2Sum, b1RefIndex, b2RefIndex, deltaRefIndex, "." # difference ratio is null
    elif b1Sum == 0 or b2Sum == 0:
        return b1Sum, b2Sum, b1RefIndex, b2RefIndex, deltaRefIndex, 1 # difference ratio is 1
    else:
        proportions = {
            "b1": [
                alleleCount / b1Sum
                for alleleCount in b1Count.values()
            ],
            "b2": [
                alleleCount / b2Sum
                for alleleCount in b2Count.values()
            ]
        }
        
        # Derive our difference ratio value
        proportionCommon = sum([
            min(proportions["b1"][x], proportions["b2"][x])
            for x in range(len(proportions["b1"]))
        ])
        differenceRatio = 1 - proportionCommon
        
        # Return the values
        return b1Sum, b2Sum, b1RefIndex, b2RefIndex, deltaRefIndex, differenceRatio

def main():
    # User input
    usage = """%(prog)s receives several files associated with a per-sample SNP
    prediction project dealing with two segregating populations. It produces
    an output file similar to what's seen with QTL-Seq.
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
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse metadata file
    metadataDict = parse_pops_metadata(args.metadataFile)
    
    # Parse VCF data for outlier SNP genotypes
    snpGenotypes = haplotypes.get_genotypes_from_vcf(args.vcfFile, snpPositions=None, imputeMissing=False)
    
    # Format output table
    with open(args.outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "CHROM", "POSI", "variant", "bulk1_alleles",
            "bulk2_alleles", "bulk1_refIndex",
            "bulk2_refIndex", "delta_refIndex",
            "differenceRatio"
        ])))
        
        # Write content lines
        for contig, posDict in snpGenotypes.items():
            for pos, snpDict in posDict.items():
                # Figure out what type of variant this is
                if any([ x == "." for x in snpDict["ref_alt"] ]):
                    variant = "indel"
                elif len(snpDict["ref_alt"][0]) == len(snpDict["ref_alt"][1]):
                    variant = "snp"
                else:
                    variant = "indel"
                
                # Split sample genotypes into bulk1 and bulk2
                bulk1 = [ snpDict[sample] for sample in metadataDict["bulk1"] if sample in snpDict ]
                bulk2 = [ snpDict[sample] for sample in metadataDict["bulk2"] if sample in snpDict ]
                
                # Calculate difference ratio
                bulk1_alleles, bulk2_alleles, bulk1_refIndex, \
                    bulk2_refIndex, delta_refIndex, \
                    differenceRatio = calculate_snp_indices(bulk1, bulk2)
                
                # Format and write output line
                outputLine = f"{contig}\t{pos}\t{variant}\t{bulk1_alleles}\t" + \
                    f"{bulk2_alleles}\t{bulk1_refIndex}\t{bulk2_refIndex}\t" + \
                    f"{delta_refIndex}\t{differenceRatio}\n"
                
                fileOut.write(outputLine)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
