#! python3
# filter_vcf_to_match_snp_index.py
# Script to enable easy filtering of the qtlseq.vcf file so as to
# only contain the SNPs identified in the snp_index.p##.tsv file.

import os, argparse, gzip
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.snpIndexFile):
        print('I am unable to locate the SNP index p## file (' + args.snpIndexFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the QTLseq VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def parse_qtlseq_snp_index_file_to_set(snpIndexFile):
    '''
    This function will parse a a snp_index.p##.tsv file produced by QTL-seq
    into an OrderedDict containing all relevant values.
    
    Parameters:
        snpIndexFile -- a string indicating the file location of the QTL-seq
                        snp_index.p##.tsv file.
    Returns:
        snpIndexSet -- a set containing "-" concatenated values of the CHROM
                       and POSI values from the snp_index.p##.tsv file e.g.,
                       "chr1-10191" or "utg00001-1923402". It will ALSO
                       contain just the CHROM values to check to see if that
                       contig shows up (used for filtering the comment lines)
    '''
    snpIndexSet = set()
    with open(snpIndexFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            # Handle header line
            if firstLine is True:
                firstLine = False
                continue
            # Handle content lines
            else:
                chrom, pos = line.rstrip("\r\n ").split("\t")[0:2]
                snpIndexSet.add(chrom)
                snpIndexSet.add(f"{chrom}-{pos}")
    return snpIndexSet

def main():
    # User input
    usage = """%(prog)s filters a QTL-seq VCF file down to only those SNPs
    included in the snp_index.p##.tsv file produced by QTL-seq.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="snpIndexFile",
                   required=True,
                   help="Input SNP index p## file produced by QTL-seq")
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input VCF file produced by QTL-seq")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    args = p.parse_args()
    validate_args(args)
    
    # Parse SNP index file
    snpIndexSet = parse_qtlseq_snp_index_file_to_set(args.snpIndexFile)
    
    # Write new VCF file
    with open_vcf_file(args.vcf) as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            # Handle comment lines
            if line.startswith("##"):
                if "contig=" not in line:
                    fileOut.write(line)
                else:
                    contigID = line.split(",length")[0].split("ID=")[1]
                    if contigID in snpIndexSet:
                        fileOut.write(line)
            elif line.startswith("#"):
                fileOut.write(line)
            # Handle content lines
            if "-".join(line.split("\t")[0:2]) in snpIndexSet:
                fileOut.write(line)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
