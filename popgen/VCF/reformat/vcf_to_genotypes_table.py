#! python3
# vcf_to_genotypes_table.py
# A script to receive a VCF file and reformat it into a
# table with ref/alt genotypes rather than digit/digit genotypes.
# Useful for manual inspection.

import os, argparse, gzip
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.inputVCF):
        raise FileNotFoundError(f"I am unable to locate the variant calls VCF file ({args.inputVCF})")
    if (args.sampleOrderFile != None) and (not os.path.isfile(args.sampleOrderFile)):
        raise FileNotFoundError(f"I am unable to locate the sample order file ({args.sampleOrderFile})")
    # Validate output file location
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"The output file ({args.outputFileName}) already exists. " +
                              "Please specify a unique file name.")

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def parse_sample_order(filename):
    sampleOrder = []
    with open(filename, "r") as fileIn:
        for line in fileIn:
            sampleOrder.append(line.strip())
    return sampleOrder

def recode_genotypes_as_refalt(sampleData, refAlt, gtIndex, sampleNames):
    '''
    Parameters:
        sampleData -- a list of strings, each representing a sample's genotype data
        refAlt -- a list of reference and alternate allele(s) as strings
        gtIndex -- the index of the GT field in the sample data format
        sampleNames -- a list of sample names corresponding to the sampleData
    Returns:
        sampleData -- a dictionary where keys are sample names and values are their
                      recoded genotypes in ref/alt format (as opposed to digit/digit)
    '''
    reformattedSampleData = []
    for sample in sampleData:
        # Extract the genotype field
        fields = sample.split(":")
        gt = fields[gtIndex]
        
        # Recode the genotype
        if "." in gt:
            reformattedSampleData.append(".")
        else:
            alleles = gt.split("/")
            reformattedGT = "/".join([ refAlt[int(allele)] for allele in alleles ])
            reformattedSampleData.append(reformattedGT)
    
    sampleData = {sampleNames[i]: reformattedSampleData[i] for i in range(len(sampleNames))}
    return sampleData

def parse_vcf(filename, sampleOrder=None):
    '''
    Parameters:
        filename -- a string containing the path to the VCF file
        sampleOrder -- a list of sample names in the desired order, or None to use VCF header order
    Returns:
        vcfDict -- a dictionary where keys are variant identifiers (chrom:pos:ref:alt)
                   and values are dictionaries with keys "chrom", "pos", "ref", "alt",
                   and "sampleData" (a dictionary of sample names to their genotypes)
    '''
    vcfDict = {}
    with open_vcf_file(filename) as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#CHROM"):
                sampleNames = line.strip().split("\t")[9:]
                if sampleOrder != None and sampleOrder != []:
                    soSet = set(sampleOrder)
                    snSet = set(sampleNames)
                    
                    if soSet != snSet:
                        errorMessage = []
                        soDiff = soSet.difference(snSet)
                        if len(soDiff) > 0:
                            errorMessage.append(f"Sample order contains samples not in VCF: \"{', '.join(soDiff)}\"")
                        
                        snDiff = snSet.difference(soSet)
                        if len(snDiff) > 0:
                            errorMessage.append(f"VCF contains samples not in sample order: \"{', '.join(snDiff)}\"")
                        
                        raise ValueError("; and ".join(errorMessage))
            if line.startswith("#"):
                continue
            
            # Process variant lines
            sl = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, formatField = sl[0:9]
            refAlt = [ref, *alt.split(",")]
            gtIndex = formatField.split(":").index("GT")
            sampleData = sl[9:]
            
            # Format sample data into a dict for this variant
            sampleData = recode_genotypes_as_refalt(sampleData, refAlt, gtIndex, sampleNames)
            
            # Store the sample data in the dictionary
            key = f"{chrom}:{pos}:{ref}:{alt}" # ensure unique key
            assert key not in vcfDict, f"Duplicate key found: {key}"
            vcfDict[key] = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "sampleData": sampleData
            }
    return vcfDict

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF file and reformats it into a table with ref/alt genotypes
    in human-readable format, rather than digit/digit genotypes.
    """
    p = argparse.ArgumentParser(description=usage)
    # Req
    p.add_argument("-i", dest="inputVCF",
                   required=True,
                   help="Input the VCF file to convert")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for genotypes table")
    # Opts
    p.add_argument("--sampleOrder", dest="sampleOrderFile",
                   required=False,
                   help="""Optionally, specify a text file listing the order of samples
                   (as columns) in the output table; default is to use the VCF header order""",
                   default=None)
    args = p.parse_args()
    validate_args(args)
    
    # Parse sample order
    if args.sampleOrderFile is None:
        sampleOrder = None
    else:
        sampleOrder = parse_sample_order(args.sampleOrderFile)
    
    # Parse VCF file
    vcfDict = parse_vcf(args.inputVCF, sampleOrder)
    
    # Write output file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#CHROM\tPOS\tREF\tALT\t" + "\t".join(sampleOrder) + "\n")
        
        # Write each variant
        for key in vcfDict:
            variant = vcfDict[key]
            sampleData = [ variant["sampleData"][sample] for sample in sampleOrder ] # reorder sample data
            line = f"{variant['chrom']}\t{variant['pos']}\t{variant['ref']}\t{variant['alt']}\t" + "\t".join(sampleData)
            fileOut.write(line + "\n")
    
    print("Program finished successfully!")

if __name__ == "__main__":
    main()
