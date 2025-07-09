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
    # Validate numeric arguments
    if args.ploidy < 1:
        raise ValueError("Ploidy must be a positive integer greater than 0.")
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

def recode_genotypes_as_msa(sampleData, refAlt, gtIndex, sampleNames, ploidy=2):
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
    # Get the length of the longest allele
    alleleLength = max([ len(allele) for allele in refAlt ])
    
    reformattedSampleData = []
    for sample in sampleData:
        # Extract the genotype field
        fields = sample.split(":")
        gt = fields[gtIndex]
        
        # Recode the genotype
        if "." in gt:
            reformattedSampleData.append("-"*alleleLength*ploidy)  # Missing genotype
        else:
            alleles = gt.split("/")
            msaGT = "".join([ refAlt[int(allele)].ljust(alleleLength, "-") for allele in alleles ])
            reformattedSampleData.append(msaGT)
    
    sampleData = {sampleNames[i]: reformattedSampleData[i] for i in range(len(sampleNames))}
    return sampleData

def parse_vcf_to_msaDict(filename, ploidy=2):
    '''
    Parameters:
        filename -- a string containing the path to the VCF file
    Returns:
        msaDict -- a dictionary where keys are sample names and values are lists of
                   genotypes in MSA format (ref/alt with '-' padding for length equivalence)
    '''
    msaDict = {}
    with open_vcf_file(filename) as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#CHROM"):
                sampleNames = line.strip().split("\t")[9:]
                msaDict = {sample: [] for sample in sampleNames}
            if line.startswith("#"):
                continue
            
            # Process variant lines
            sl = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, formatField = sl[0:9]
            refAlt = [ref, *alt.split(",")]
            gtIndex = formatField.split(":").index("GT")
            sampleData = sl[9:]
            
            # Format sample data into a dict for this variant
            sampleData = recode_genotypes_as_msa(sampleData, refAlt, gtIndex, sampleNames, ploidy=ploidy)
            
            # Store the sample data in the dictionary
            for sampleID, msaGT in sampleData.items():
                msaDict[sampleID].append(msaGT)
    
    return msaDict

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF file and reformats it into a FASTA with genotypes encoded
    in MSA-style format amenable to phylogenetic analysis.
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
    p.add_argument("--ploidy", dest="ploidy",
                   required=False,
                   type=int,
                   help="Optionally, indicate the ploidy of the samples (default: 2)",
                   default=2)
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF file
    msaDict = parse_vcf_to_msaDict(args.inputVCF, ploidy=args.ploidy)
    
    # Write output file
    with open(args.outputFileName, "w") as fileOut:
        for sampleID, msaGTs in msaDict.items():
            seq = "".join(msaGTs)
            fileOut.write(f">{sampleID}\n{seq}\n")
    
    print("Program finished successfully!")

if __name__ == "__main__":
    main()
