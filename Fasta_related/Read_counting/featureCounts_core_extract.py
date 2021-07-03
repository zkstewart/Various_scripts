#! python3
# featureCounts_core_extract.py
# Script to extract the mapped reads associated with each feature
# from the results of featureCounts when run with output to CORE format.

import os, argparse
import gzip
from Bio import SeqIO

# Define functions
def validate_args(args):
    # Validate that all arguments have been provided
    for key, value in vars(args).items():
        if value == None:
            print(key + ' argument was not specified. Fix this and try again.')
            quit()
    # Validate input file locations
    if not os.path.isfile(args.coreFile):
        print('I am unable to locate the CORE file (' + args.coreFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastqFile[0]):
        print('I am unable to locate the FASTQ file (' + args.fastqFile[0] + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastqFile[1]):
        print('I am unable to locate the FASTQ file (' + args.fastqFile[1] + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()

def core_to_dict(coreFile):
    coreDict = {}
    with open(coreFile, "r") as fileIn:
        for line in fileIn:
            # Parse line and skip irrelevant ones
            sl = line.rstrip("\r\n").split("\t")
            if sl[1] != "Assigned":
                continue
            # Extract details
            read = sl[0]
            feature = sl[3]
            # Add to dict
            coreDict[read] = feature
            # if feature not in coreDict:
            #     coreDict[feature] = [read]
            # else:
            #     coreDict[feature].append(read)
    return coreDict

def fastq_maybe_gz_handle(fastqFile):
    if fastqFile.endswith(".gz"):
        return gzip.open(fastqFile, "rt")
    else:
        return open(fastqFile, "r")

def coredict_and_fastq_to_feature_fastqs(coreDict, fastqFile1, fastqFile2):
    with fastq_maybe_gz_handle(fastqFile1) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if (record.id in coreDict):
                with open("{0}_R1.fastq".format(coreDict[record.id]), "a") as fileOut:
                    SeqIO.write(record, fileOut, "fastq")
    
    with fastq_maybe_gz_handle(fastqFile2) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if (record.id in coreDict):
                with open("{0}_R2.fastq".format(coreDict[record.id]), "a") as fileOut:
                    SeqIO.write(record, fileOut, "fastq")

def main():
    # User input
    usage = """%(prog)s reads in CORE formatted file produced by featureCounts (-R CORE) alongside
    the FASTQ used for mapping and, for each feature, outputs a FASTQ file of the reads that mapped
    to said feature. The input FASTQ can be in gzip format.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-c", dest="coreFile",
        help="Input CORE file name")
    p.add_argument("-fq", dest="fastqFile", nargs=2,
        help="Input FASTQ file name")
    args = p.parse_args()
    validate_args(args)

    # Parse CORE file
    coreDict = core_to_dict(args.coreFile)

    # Parse FASTQ file and produce each output
    coredict_and_fastq_to_feature_fastqs(coreDict, args.fastqFile[0], args.fastqFile[1])

if __name__ == "__main__":
    main()
