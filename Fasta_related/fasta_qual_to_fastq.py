#! python3

import os, argparse
from Bio import SeqIO

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.fastaFile):
                print('I am unable to locate the fasta file (' + args.fastaFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.fastqFile):
                print('I am unable to locate the fastq file (' + args.fastqFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate the index value
        if args.index < 0:
                print('Index value must be a number greater than 0.')
                quit()
        with open(args.fastaFile, 'r') as fileIn:
                for line in fileIn:
                        firstFasta = line.strip('>\r\n')
                        break
        if args.index == 0:
                qualCode = firstFasta
        else:
                splitLine = firstFasta.split(' ')
                if len(splitLine) < args.index:
                        print('The index value you specified seems to be incorrect. The split line looks like this.')
                        print(splitLine)
                        print('You specified an index of ' + str(args.index) + ' which exceeds the length of this list.')
                        print('Make sure your index is correct and try again.')
                        quit()
                qualCode = splitLine[args.index-1]
        found = False
        with open(args.fastqFile, 'r') as fileIn:
                for line in fileIn:
                        if not line.startswith('@'):
                                continue
                        qualID = line.strip('@\r\n')
                        if qualID == qualCode:
                                found = True
                                break
        if found == False:
                print('The code in your first fasta entry (' + qualCode + ') couldn\'t be found in your fastq file!')
                print('Check that your index value is correct, and see if this code is actually in the fastq file. If it\'s missing, something is wrong.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

## Output related
def demultiplexed_fasta_fix(fastaFile, origFastq, idIndex, outputFileName):
        fastaParse = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        fastqParse = SeqIO.to_dict(SeqIO.parse(open(origFastq, 'r'), 'fastq'))  # oof, memory hit here... will only run on a HPC now probably
        with open(outputFileName, 'w') as fileOut:
                for record in fastaParse:
                        # Parse the qual identifying string
                        seqid = record.description
                        if idIndex == 0:
                                qualCode = seqid
                        else:
                                splitID = seqid.split(' ')
                                qualCode = splitID[idIndex-1]   # Assume the input idIndex is 1-based
                        # Process the qual line to make it match the fasta sequence if it was altered
                        qualLine = fasta_fastq_qual_match(str(record.seq), fastqParse[qualCode])
                        # Make sure it matches
                        assert len(record) == len(qualLine)
                        # Produce output
                        fileOut.write('>' + seqid + '\n' + str(record.seq) + '\n+\n' + qualLine + '\n')

def fasta_fastq_qual_match(fastaSeq, origFastqRecord):
        # Process the fastq record to retrieve lines
        fqFormat = origFastqRecord.format("fastq")
        fqLines = fqFormat.split('\n')
        fqSeq = fqLines[1]
        fqQual = fqLines[3]
        # Find the new fasta sequence in the original
        newStart = fqSeq.find(fastaSeq)
        newEnd = newStart + len(fastaSeq)
        # Extract the fastq line which matches this sequence
        newQual = fqQual[newStart:newEnd]
        # Make sure its OK
        newSeq = fqSeq[newStart:newEnd]
        assert newSeq == fastaSeq       # If this is True, then we know that the qual region should match the sequence region properly
        # Return
        return newQual

##### USER INPUT SECTION

usage = """%(prog)s reads in a demultiplexed .fasta and original .fastq format file pair and associates
the correct quality scores to the demultiplexed fasta file, producing a new .fastq output file. The two
file do not need to be sorted, however, the sequence ID of the original fastq sequence must be within the sequence ID
of the demultiplexed .fasta file. If the .fasta file's sequence ID has space characters in it, you
can specify the location of the identifying string by specifying its index in 1-based
notation (e.g., for "STP10C_1 9QUDO:00017:00096 orig_bc=AACCATCCGCT", index 1 would correspond
to STP10C_1 and index 2 would correspond to 9QUDO:00017:00096 which would be the identifier
to pair the .fasta to the .fastq file). Note that because this script works with unsorted data, it needs
to load the full .fastq file into memory which is time and memory consuming.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-f", dest="fastaFile",
		  help="Specify new fasta file")
p.add_argument("-q", dest="fastqFile",
		  help="Specify original fastq file")
p.add_argument("-i", dest="index", type=int,
		  help="""Specify the index of the identifying string. If there are no space characters, then you can leave this at its
          default value of 1. If there are space characters and you want to keep them in for the matching process, enter 0 here instead.""", default=1)
p.add_argument("-o", "-output", dest="outputFileName",
	       help="Output file name")

args = p.parse_args()
validate_args(args)

# Parse fasta file and produce output
demultiplexed_fasta_fix(args.fastaFile, args.fastqFile, args.index, args.outputFileName)
