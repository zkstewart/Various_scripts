#! python3
# htseq_tabulate.py
# Script to combine multiple HTSeq outputs into a single table file for use with 
# edgeR / DESeq2

# Import necessary packages
import os, argparse

# Define functions
def validate_args(args):
    # Validate that all arguments are provided
    for key, value in vars(args).items():
        if value == None:
            print(key + " argument was not specified. Fix this and try again.")
            quit()
    
    # Validate input file locations
    for file in args.inputFiles:
        if not os.path.isfile(file):
            print("I am unable to locate the HTSeq count file ({0})".format(file))
            print("Make sure you\'ve typed the file name or location correctly and try again.")
            quit()

#### USER INPUT SECTION
usage = """%(prog)s will produce a single tabulated file from multiple HTSeq output files. This
table is formatted for DGE analysis with edgeR. The header of this file will have columns labelled
according to each filename up to the first '.' (e.g., 'SITE1.htseq.counts' will have its column
labelled 'SITE1').
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="inputFiles", nargs="+",
    help="Specify each input file to be tabulated")
p.add_argument("-o", dest="output",
    help="output file name")
args = p.parse_args()

# Tabulate in list data structure
for i in range(len(args.inputFiles)):
    # Handle the first file
    if i == 0:
        # Turn the entire file into a list
        with open(args.inputFiles[i], 'r') as fileIn:
            htseqFile = fileIn.read()
        combinedList = htseqFile.split('\n')
        # Clean up potential blank lines in the list
        while combinedList[-1] == '' or combinedList[-1].replace('\t', '') == '':
            del combinedList[-1]
        # Separate the list into gene IDs on the left, and counts on the right
        for i in range(len(combinedList)):
            combinedList[i] = combinedList[i].split('\t')
    # Handle all subsequent files
    else:
        with open(args.inputFiles[i], 'r') as fileIn:
            for line in fileIn:
                if line == "" or line.replace("\t", "") == "":
                    continue
                row = line.split("\t")
                combinedList[i].append(row[1])
    # Alert user to script progress
    print('Finished file (' + args.inputFiles[i] + ')')

# Format output file header
sampleNames = ['gene_id']
for name in args.inputFiles:
    sampleNames.append(name.split('.')[0])
header = '\t'.join(sampleNames) + '\n'

# Format each line for output to text
for i in range(len(combinedList)):
    combinedList[i] = '\t'.join(combinedList[i])
output = '\n'.join(combinedList)

# Write output to file
with open(args.output, 'w') as fileOut:
    output = header + output
    fileOut.write(output)

print('All done')
