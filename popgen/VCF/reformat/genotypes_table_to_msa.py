#! python3
# genotypes_table_to_msa.py
# Converts a genotypes table (as produced by 'vcf_to_genotypes_table.py')
# and converts it to a multiple sequence alignment (MSA) format which can
# be used for phylogenetic analysis.

import os, argparse, gzip
import pandas as pd
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.inputGenotypesTable):
        raise FileNotFoundError(f"I am unable to locate the genotypes table file ({args.inputGenotypesTable})")
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

def parse_genotypes_table(filename, ploidy=2):
    '''
    Parameters:
        filename -- a string representing the path to the genotypes table file.
    Returns:
        samplesDict -- a dictionary with structure like:
                       {
                           "sample1": [["A", "A"], ["C", "T"], ...],
                           "sample2": [["G", "G"], ["-", "-"], ...],
                           ...
                       }
    '''
    df = pd.read_csv(filename, sep="\t", comment=None)
    
    # Check for ref/alt length
    if not all(df["REF"].apply(lambda x: len(x)) == df["ALT"].apply(lambda x: len(x))):
        raise ValueError("The reference and alternate alleles must have the same length for all variants; " +
                         "you will need to use vcf_to_msa.py to handle this scenario.")
    
    # Store genotypes in a dictionary
    samplesDict = {}
    for col in df.columns[4:]: # Skip first four columns (chromosome, position, ref, alt)
        genotypes = df[col].tolist()
        alleles = [ gt.split("/") if not "." in gt else ["-"]*ploidy for gt in genotypes ]
        samplesDict[col] = alleles
    return samplesDict

def convert_to_fasta(samplesDict, outputFileName):
    '''
    Parameters:
        samplesDict -- a dictionary with structure like:
                       {
                           "sample1": [["A", "A"], ["C", "T"], ...],
                           "sample2": [["G", "G"], ["-", "-"], ...],
                           ...
                       }
        outputFileName -- a string representing the path to the output FASTA file.
    '''
    with open(outputFileName, "w") as fileOut:
        for sample, alleles in samplesDict.items():
            sequence = "".join([ "".join(allele) for allele in alleles ])
            fileOut.write(f">{sample}\n{sequence}\n")

## Main
def main():
    # User input
    usage = """%(prog)s reads in a genotypes table file (as produced by 'vcf_to_genotypes_table.py')
    and converts it to a multiple sequence alignment (MSA) format which can be used for phylogenetic analysis.
    The output is a FASTA file with each sample as a separate sequence; assuming diploid data,
    each genotype is represented as two alleles per position. If variants include indels, you should
    make sure to align the output FASTA prior to phylogenetic analysis.
    """
    p = argparse.ArgumentParser(description=usage)
    # Req
    p.add_argument("-i", dest="inputGenotypesTable",
                   required=True,
                   help="Input genotypes table file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name MSA FASTA file")
    # Opts
    p.add_argument("--ploidy", dest="ploidy",
                   required=False,
                   type=int,
                   help="Optionally, indicate the ploidy of the samples (default: 2)",
                   default=2)
    args = p.parse_args()
    validate_args(args)
    
    # Parse genotypes table
    gtSamplesDict = parse_genotypes_table(args.inputGenotypesTable, ploidy=args.ploidy)
    
    # Write to MSA FASTA file
    convert_to_fasta(gtSamplesDict, args.outputFileName)
    
    print("Program finished successfully!")

if __name__ == "__main__":
    main()
