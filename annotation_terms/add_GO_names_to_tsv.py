#! python3
# add_GO_names_to_tsv
# Script to take in a TSV file containing GO codes and insert a column with the GO names
# associated with those codes. This script is format agnostic and intended to be flexible.

# Load packages
import os, argparse, sys
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_GO

# Define functions
def validate_args(args):
    # Validate input locations
    if not os.path.isfile(args.inputTSV):
        print('I am unable to locate the input TSV file (' + args.inputTSV + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for oboFile in args.goOboFiles:
        if not os.path.isfile(oboFile):
            print('I am unable to locate the input GO.obo file (' + oboFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate numeric arguments
    if args.locationIndex < 1:
        print("locationIndex must be 1 or greater")
        quit()
    if args.insertionIndex < 1:
        print("insertionIndex must be 1 or greater")
        quit()
    # Handle file overwrites
    if os.path.exists(args.outputFileName):
        print(f"{args.outputDirectory} already exists, and I will not overwrite it.")
        print("Delete/move whatever exists here, or specify a different output name when you try again.")
        quit()

def parse_annotation_file(annotationTSV, annotationDelimiter="; ", blankAnnotation="0", hasHeader=False):
    '''
    Parameters:
        annotationTSV -- a string representing the path to the annotation TSV file
        annotationDelimiter -- OPTIONAL; a string representing the delimiter used to separate
                               annotation terms in the annotation TSV file; default == "; "
        blankAnnotation -- OPTIONAL; a string representing the value to expect for sequences without
                           annotations; default == "0"
        hasHeader -- OPTIONAL; a boolean indicating if the annotation TSV file has a header;
                     default == False
    Returns:
        annotationDict -- a dictionary with sequence IDs as keys and annotation terms as
                          a set of annotation strings
    '''
    annotationDict = {}
    with open(annotationTSV, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            # Handle potential header
            if firstLine:
                firstLine = False
                if hasHeader:
                    continue
            # Handle content lines
            seqid, annotations = line.strip("\r\n ").split("\t")
            if annotations == blankAnnotation:
                annotations = set()
            else:
                annotations = set(annotations.split(annotationDelimiter))
            annotationDict[seqid] = annotations
    return annotationDict

def parse_codes_from_df(df, goDelimiter="; ", blankCharacter="0", locationIndex=0):
    ''''
    Parameters:
        df -- a pandas.DataFrame object of a file containing GO codes delimited by goDelimiter
        goDelimiter -- OPTIONAL; a string representing the delimiter used to separate GO codes
                       in the DataFrame; default == "; "
        blankCharacter -- OPTIONAL; a string representing the value to expect for sequences without
                          annotations; default == "0"
        locationIndex -- OPTIONAL; an integer representing the 0-indexed column where the GO codes
                         are located; default == 0
    '''
    # Obtain GO codes as a list
    goCodesList = df.iloc[:, locationIndex].to_list()
    
    # Separate GO code by delimiter
    goCodesList = [
        x.split(goDelimiter) if x != blankCharacter else [blankCharacter]
        for x in goCodesList
    ]
    
    return goCodesList

def main():
    # User input
    usage = """%(prog)s will receive a file containing GO codes in the indicated column and will insert
    an additional column containing the GO names associated with those codes. To make this work in a format
    agnostic way, you must 1) indicate whether your input has a header, 2) specify the column number where
    the GO codes are located, 3) specify the delimiter of the GO codes, and 4) specify the index where you'd
    like the new column to be inserted (0-indexed). The new header will be 'GO_names'.
    NOTE 1: this script can handle header-ed files, but only if they contain a single line (and if you
    specify --hasHeader). NOTE 2: this script can accept multiple GO.obo files to parse, and you are
    encouraged to use the GO.obo file that you used to generate the annotations initially, as well as the
    most recent version of the GO.obo file.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputTSV",
                   required=True,
                   help="Specify input TSV file")
    p.add_argument("-g", dest="goOboFiles",
                   required=True,
                   nargs="+",
                   help="Specify one or more go.obo files to parse")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the modified TSV with GO names column")
    # Opts (file format)
    p.add_argument("--location", dest="locationIndex",
                   required=False,
                   type=int,
                   help="""Preferably, specify the 1-based location index of the column containing
                   the GO codes; default == 1 (1st column)""",
                   default=1)
    p.add_argument("--insertion", dest="insertionIndex",
                   required=False,
                   type=int,
                   help="""Preferably, specify the 1-based index for where you'd like the GO names
                   column to be inserted; default == 2 (as the 2nd column)""",
                   default=2)
    p.add_argument("--delimiter", dest="annotationDelimiter",
                   required=False,
                   help="""Optionally, specify the delimiter used in the 
                   annotation TSV file separating terms; default == '; '""",
                   default="; ")
    p.add_argument("--blank", dest="blankCharacter",
                   required=False,
                   help="""Optionally, specify the character used in the 
                   input TSV file to indicate no annotations; default == '0'""",
                   default="0")
    p.add_argument("--hasHeader", dest="hasHeader",
                   required=False,
                   action="store_true",
                   help="Optionally indicate if the annotation TSV file has a SINGLE header line",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GO obo
    go = ZS_GO.GO(args.goOboFiles)
    
    # Parse the input TSV
    df = pd.read_csv(args.inputTSV, sep="\t", comment=None,
                     header=0 if args.hasHeader else None)
    
    # Parse the GO codes out of the DataFrame
    goCodesList = parse_codes_from_df(df, goDelimiter="; ", blankCharacter="0",
                                      locationIndex=args.locationIndex-1) # make it 0-based
    
    # Convert GO codes to GO names
    goNames = go.convert_codes_to_names(goCodesList)
    
    # Flatten the list of lists and drop Nones and blanks
    goNames = [
        args.annotationDelimiter.join([ x for x in sublist if x != None and x != args.blankCharacter ])
        for sublist in goNames
    ]
    
    # Replace empty strings with blankCharacter
    goNames = [ args.blankCharacter if x == "" else x for x in goNames ]
    
    # Add column to DataFrame
    df.insert(args.insertionIndex-1, "GO_names", goNames)
    
    # Write output file
    df.to_csv(args.outputFileName,
          sep="\t", index=False, lineterminator="\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
