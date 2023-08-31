#! python3
# fasta_extractSubset.py
# Extract a subset of sequences from an input fasta or fastq file
# Behaviour of subsetting can be specified, and paired fastq file
# subsetting is supported.

import os, argparse, random, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Define functions for later use
def validate_args(args):
    # Ensure inputFastas numbers are sensible
    if len(args.inputFastas) == 0:
        print('Detected no input fasta files. Specify these files in the command-line.')
        quit()
    elif len(args.inputFastas) > 2:
        print('Detected more than 2 input fasta files. This script does not allow that to happen. Fix your inputs and try again.')
        quit()
    # Validate file locations
    for entry in args.inputFastas:
        if not os.path.isfile(entry):
            print('I am unable to locate the input FASTA file (' + entry + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Handle file type processing
    if len(args.inputFastas) == 2:
        print('You\'ve specified two files, which means I\'m assuming it\'s paired fastq files. Validating this...')
        lineNum1 = 0    # We'll use this value to count the number of lines in the file; if they differ, they're unlikely to be properly paired fastq files
        firstChar1 = None
        for line in file_iterator(args.inputFastas[0]):
            if firstChar1 == None:
                firstChar1 = line[0]
                if firstChar1 != '@':
                    break
            lineNum1 += 1
        lineNum2 = 0
        firstChar2 = None
        for line in file_iterator(args.inputFastas[1]):
            if firstChar2 == None:
                firstChar2 = line[0]
                if firstChar2 != '@':
                    break
            lineNum2 += 1
        # Ensure file type is correct
        if firstChar1 == '@' and firstChar2 == '@':
            print('Tentatively believe these are fastq files (both files start with "@"). If this isn\'t true, expect an error.')
            seqType = 'fastq'
        else:
            print('I don\'t think these are fastq files! One or both files do not start with "@" as I would expect. Fix your inputs to continue.')
            quit()
        # Ensure file length is correct
        if lineNum1 != lineNum2:
            print('You\'ve provided two fastq files, but they don\'t seem to have the same number of lines! This is bad, and could cause errors or it could result in nonsense results where pairing is disrupted. Fix your inputs and try again.')
            quit()
    else:
        print('You\'ve specified one file. Will check to find its sequence type...')
        for line in file_iterator(args.inputFastas[0]):
            firstChar1 = line[0]
            break
        
        if firstChar1 == '@':
            print('Tentatively believe this is a fastq file (file starts with "@"). If this isn\'t true, expect an error.')
            seqType = 'fastq'
        elif firstChar1 == '>':
            print('Tentatively believe this is a fasta file (file starts with ">"). If this isn\'t true, expect an error.')
            seqType = 'fasta'
        else:
            print('I don\'t recognise this file! It should start with "@" (fastq) or ">" (fasta). Fix your inputs to continue.')
            quit()
    # Ensure input number is sensible
    if args.sampleNum < 1:
        print('Sample number argument must be >= 1')
        quit()
    # Handle output file names
    if len(args.outputFileNames) != len(args.inputFastas):
        print("{0} input FASTAs given, but a non-equal number of output files given ({1})".format(len(args.inputFastas), len(args.outputFileNames)))
        print("Correct this issue and try again.")
        quit()
    for outputFileName in args.outputFileNames:
        if os.path.isfile(outputFileName):
            print('Specified file name already exists (' + outputFileName + ')')
            print('This program won\'t overwrite an existing file; provide a new name and try again.')
            quit()
    # Ensure argument combinations are sensible
    if args.sampleType == "longest" and seqType == "fastq" and len(args.inputFastas) == 2:
        print("Sampling the longest FASTQ sequences from paired files usually doesn't make sense...")
        print("For that reason, we don't support it here. Sorry about that.")
        quit()
    return seqType

def is_gz(file):
    with gzip.open(file, 'r') as fh:
        try:
            fh.read(1)
            return True
        except OSError:
            return False

def file_iterator(fileName):
    if is_gz(fileName):
        with gzip.open(fileName, "rt") as handle:
            yield from handle
    else:
        with open(fileName, "r") as handle:
            yield from handle

def fasta_iterator(fileName):
    if is_gz(fileName):
        with gzip.open(fileName, "rt") as handle:
            yield from SimpleFastaParser(handle)
    else:
        with open(fileName, "r") as handle:
            yield from SimpleFastaParser(handle)

def fastq_iterator(fileName):
    if is_gz(fileName):
        with gzip.open(fileName, "rt") as handle:
            yield from FastqGeneralIterator(handle)
    else:
        with open(fileName, "r") as handle:
            yield from FastqGeneralIterator(handle)

def subsample_output(inputFileName, outputFileName, seqType, sampleIndices, sampleNum):
    ongoingCount = 0
    outputCount = 0 # Use this to stop program as soon as we got what we were sampling for
    with open(outputFileName, "w") as fileOut:
        # Iterate through FASTA file
        if seqType == "fasta":
            for id, seq in fasta_iterator(inputFileName):
                # If this sequence index is one we intend to sample
                if ongoingCount in sampleIndices:
                    fileOut.write(f">{id}\n{seq}\n")
                    outputCount += 1
                ongoingCount += 1
                # End iteration if we've sampled as many as desired
                if outputCount == sampleNum:
                    break
        # Iterate through FASTQ file
        else:
            for id, seq, qual in fastq_iterator(inputFileName):
                # If this sequence index is one we intend to sample
                if ongoingCount in sampleIndices:
                    fileOut.write(f"@{id}\n{seq}\n+\n{qual}\n")
                    outputCount += 1
                ongoingCount += 1
                # End iteration if we've sampled as many as desired
                if outputCount == sampleNum:
                    break

def count_num_seqs(inputFileName, seqType):
    numSeqs = 0
    if seqType == "fasta":
        for _ in fasta_iterator(inputFileName):
            numSeqs += 1
    else:
        for _ in fastq_iterator(inputFileName):
            numSeqs += 1
    return numSeqs

def get_n_longest_seq_indices(inputFileName, seqType, sampleNum):
    seqLens = [] # stores [seqID, seqLen] pairs to enable sorting by seqLen with seqID remaining paired
    ongoingCount = 0
    if seqType == "fasta":
        for id, seq in fasta_iterator(inputFileName):
            seqLens.append([ongoingCount, len(seq)])
            ongoingCount += 1
    else:
        for id, seq, qual in fastq_iterator(inputFileName):
            seqLens.append([ongoingCount, len(seq)])
            ongoingCount += 1
        
    # Return top n sequence indices
    seqLens.sort(key = lambda x: -x[1])
    return [s[0] for s in seqLens[0: sampleNum]]

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s reads in a provided individual fast(a/q) or paired fastq files
    and subsets the file(s) into a specified number of sequences. Accepts .gz files,
    but writes uncompressed files only.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", "--input", dest="inputFastas", nargs="+",
            help="""Accepts 1-2 fast(a/q) file name(s). If subsetting paired fastq files, enter both file
            names to ensure they use the same random seed. If using regular fasta, just one file at a time
            please. Will also accept a single fastq file.""")
    p.add_argument("-n", "--num", dest="sampleNum", type=int,
            required=True,
            help="number of samples to obtain")
    p.add_argument("-t", "--type", dest="sampleType",
            choices=['random', 'first', 'last', 'longest'],
            required=True,
            help="""type of sampling mechanism [first means we sample n sequences from the start of the file,
            last means we sample n sequences from the end of the file, random means random selection from all
            sequences in the file, longest means we select the n longest sequences from the file", default = 'random'""")
    p.add_argument("-o", "--output", dest="outputFileNames", nargs="+",
            required=True,
            help="""Optionally, specify the output file names; must have same number of arguments
            as the number of input files; if not provided, names will be generated for you""")
    
    # Parse arguments
    args = p.parse_args()
    seqType = validate_args(args)     # totalCount corresponds to the length of the first value in inputFastas
    
    # Derive details from the input file(s)
    numSeqs = count_num_seqs(args.inputFastas[0], seqType)
    
    # Figure out which sequences we sample
    if args.sampleNum >= numSeqs:
        print("You specified a sampling number greater than or equal to the number of sequences " +
              "in the input file. No point performing this script unless you sample less than " +
              f"{str(numSeqs)} sequences.""")
        quit()
    if args.sampleType == 'random':
        sampleIndices = random.sample(range(0, numSeqs), args.sampleNum)
    elif args.sampleType == 'first':
        sampleIndices = list(range(0, args.sampleNum))
    elif args.sampleType == 'last':
        sampleIndices = list(range(numSeqs-args.sampleNum, numSeqs))
    elif args.sampleType == 'longest':
        sampleIndices = get_n_longest_seq_indices(args.inputFastas[0], seqType, args.sampleNum)
    
    # Convert sampleIndices into a set (saves time when working with very large files)
    sampleIndices = set(sampleIndices)
    
    # Perform sampling
    for i in range(len(args.inputFastas)):
        outputFileName = args.outputFileNames[i]
        subsample_output(args.inputFastas[i], outputFileName, seqType, sampleIndices, args.sampleNum)
    
    # Done!
    print('Program exited successfully!')

if __name__ == "__main__":
    main()
