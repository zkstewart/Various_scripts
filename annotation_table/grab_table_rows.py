#! python3
# grab_table_rows
# Script to obtain rows from a table in a TSV file format based on
# row names directly obtained from a GFF3 file or from a regular text
# file.

import os, argparse, re
from Bio import SeqIO

# Define functions
def validate_args(args):
    # Validate input file arguments
    if not os.path.isfile(args.tableFileName):
        print('I am unable to locate the input TSV file (' + args.tableFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate sequence ID arguments
    if (args.idString == None or args.idString == []) and args.textFileName == None and args.gff3FileName == None and args.fastaFileName == None:
        print('You must specify at least one argument for -t, -s, -f, or -g fields; this will give us a list of sequence IDs to retrieve')
        print('Fix your inputs and try again.')
        quit()
    if args.textFileName != None:
        if not os.path.isfile(args.textFileName):
            print('I am unable to locate the input ID list file (' + args.textFileName + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.gff3FileName != None:
        if not os.path.isfile(args.gff3FileName):
            print('I am unable to locate the input GFF3 file (' + args.gff3FileName + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate numeric arguments
    if args.numHeaderLines < 0:
        print('The number of header lines must be a non-negative integer.')
        print('Fix your inputs and try again.')
        quit()
    # Handle file overwrites
    if os.path.exists(args.outputFileName):
        print(f"'{args.outputFileName}' already exists; Specify a different name or delete, move, or rename the " +
              "existing file and run the program again.")
        quit()

def text_file_to_dict(textFile):
    '''
    Parameters:
        textFile -- a string representing the path to a text file
    Returns:
        outDict -- a dictionary with the lines as keys and None as values
    '''
    outDict = {}
    with open(textFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ").lstrip(">")
            if l != "":
                assert l not in outDict, f"Duplicate entry found in text file: {l}"
                outDict.setdefault(l)
    return outDict

def gff3_file_to_dict(gff3File, gff3Types):
    '''
    Parameters:
        gff3File -- a string representing the path to a GFF3 file
        gff3Types -- a list of strings representing the feature types to obtain IDs from
                     e.g. ["gene", "mRNA"]
    Returns:
        outDict -- a dictionary with the feature IDs as keys and None as values
    '''
    idsRegex = re.compile(r"(^|;)ID=([^;]+)")
    
    outDict = {}
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            contig, source, featureType, start, end, \
                score, strand, frame, attributes \
                = line.rstrip("\t\n").split("\t")
            if featureType in gff3Types:
                try:
                    seqID = idsRegex.search(attributes).groups()[1]
                except:
                    raise ValueError(f"Unable to find feature ID in GFF3 line: {line}")
                
                assert seqID not in outDict, f"Duplicate entry found in GFF3 file: {seqID}"
                outDict.setdefault(seqID)
    return outDict

def fasta_to_dict(fastaFile):
    '''
    Parameters:
        fastaFile -- a string representing the path to a GFF3 file
    Returns:
        outDict -- a dictionary with the sequence IDs as keys and None as values
    '''
    outDict = {}    
    records = SeqIO.parse(fastaFile, "fasta")
    for record in records:
        outDict.setdefault(record.id)
    return outDict

def main():
    # User input
    usage = """%(prog)s reads in sequence IDs from one or more sources (text, GFF3, FASTA,
    command-line) and retrieves the indicated rows from a TSV table file. The row
    identifier is expected to be in the first column of the table. The output file
    will contain the same header as the input file, followed by the rows that match
    the sequence IDs.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="tableFileName",
                   required=True,
                   help="Input TSV file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name")
    # Opts (sequence IDs)
    p.add_argument("-t", dest="textFileName",
                   required=False,
                   help="Optionally specify the name of a text file listing sequence IDs")
    p.add_argument("-g", dest="gff3FileName",
                   required=False,
                   help="Input GFF3 file name")
    p.add_argument("-f", dest="fastaFileName",
                   required=False,
                   help="Input fasta file name")
    p.add_argument("-s", dest="idString",
                   required=False,
                   nargs = "+",
                   help="Optionally input sequence IDs as text here; entries will be separated at space characters",
                   default=[])
    # Opts (behavioural)
    p.add_argument("--gff3Types", dest="gff3Types",
                   required=False,
                   nargs="+",
                   help="""Optionally, if providing a -g file, specify one or more feature types to obtain IDs from;
                   default is 'gene' and 'mRNA'; this is case sensitive!""",
                   default=["gene", "mRNA"])
    p.add_argument("--numHeaderLines", dest="numHeaderLines",
                   required=False,
                   type=int,
                   help="""Optionally specify the number of header lines in the TSV file; default is 1""",
                   default=1)
    
    args = p.parse_args()
    validate_args(args)
    
    # Obtain our sequence IDs
    seqIDs = {} # use a dict to prevent duplication and keep insertion order
    if args.textFileName != None:
        seqIDs = text_file_to_dict(args.textFileName)
    
    if args.gff3FileName != None:
        gff3IDs = gff3_file_to_dict(args.gff3FileName, args.gff3Types)
        for idValue in gff3IDs:
            assert idValue not in seqIDs, f"Duplicate entry found in -g file: {idValue}"
            seqIDs.setdefault(idValue)
    
    if args.fastaFileName != None:
        fastaIDs = fasta_to_dict(args.fastaFileName)
        for idValue in fastaIDs:
            assert idValue not in seqIDs, f"Duplicate entry found in -f file: {idValue}"
            seqIDs.setdefault(idValue)
    
    if args.idString != None and args.idString != []:
        for idValue in args.idString:
            idValue = idValue.lstrip(">")
            if idValue != "":
                assert idValue not in seqIDs, f"Duplicate entry found in -s input: {idValue}"
                seqIDs.setdefault(idValue)
    
    # Read in the table file, writing out filtered rows
    with open(args.tableFileName, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        ongoingCount = 0
        for line in fileIn:
            # Handle header lines
            if ongoingCount < args.numHeaderLines:
                # Alert user of the line being skipped over as a header
                abbreviatedLine = line.rstrip("\r\n")[:20] + "..." \
                    if len(line) > 20 else line.rstrip("\r\n")
                print(f"# Writing header line {ongoingCount+1} of {args.numHeaderLines}: " +
                      abbreviatedLine)
                
                # Write the line and iterate header counter
                fileOut.write(line)
                ongoingCount += 1
            
            # Handle body lines
            sl = line.rstrip("\r\n ").split("\t")
            if sl[0] in seqIDs:
                fileOut.write(line)
    
    print("Program completed successfully!")

if __name__ == '__main__':
    main()
