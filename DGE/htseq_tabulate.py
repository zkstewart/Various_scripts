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
    
    # Validate output file
    if os.path.exists(args.output):
        print("Specifed output file already exists.")
        print("This program does not allow file overwriting. \
            Specify a new name or move the existing file and try again.")
        quit()

def tabulate_htseq(htseqFiles, outputFile):
    # Extraneous HTSeq lines
    junkLines = ["__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"]
    # Format tabular data
    for i in range(len(htseqFiles)):
        # Handle the first file
        if i == 0:
            # Turn the entire file into a list
            with open(htseqFiles[i], 'r') as fileIn:
                htseqFile = fileIn.read()
            table = htseqFile.split('\n')
            # Clean up potential blank lines in the list
            while table[-1] == '' or table[-1].replace('\t', '') == '':
                del table[-1]
            # Separate the list into gene IDs on the left, and counts on the right
            for i in range(len(table)):
                table[i] = table[i].split('\t')
            # Delete extraneous HTSeq lines
            while table[-1][0] in junkLines:
                del table[-1]
        # Handle all subsequent files
        else:
            ongoingCount = 0
            with open(htseqFiles[i], 'r') as fileIn:
                for line in fileIn:
                    # Skip potential blank lines
                    if line == "" or line.replace("\t", "") == "":
                        continue
                    # Skip extraneous HTSeq lines
                    row = line.rstrip("\r\n").split("\t")
                    if row[0] in junkLines:
                        continue
                    # Format line for addition to table
                    table[ongoingCount].append(row[1])
                    ongoingCount += 1
    
    # Format output file header
    sampleNames = ['gene_id']
    for name in htseqFiles:
        baseName = os.path.basename(name)
        sampleNames.append(baseName.split('.')[0])
    table.insert(0, sampleNames)
    
    # Format each line for output to text
    for i in range(len(table)):
        table[i] = '\t'.join(table[i])

    # Write output to file
    with open(outputFile, 'w') as fileOut:
        fileOut.write("\n".join(table))

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

# Run program
tabulate_htseq(args.inputFiles, args.output)

# Alert user to successful program completion
print('Program completed successfully!')
