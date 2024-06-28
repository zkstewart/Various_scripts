#! python3
# merge_annotation_tsvs
# Script to take in a multiple annotation files linking shared IDs (left column)
# to annotation terms and combines these annotations to the shared ID level.

# Load packages
import os, argparse

def validate_args(args):
    # Validate input locations
    for annotationTSV in args.annotationTSVs:
        if not os.path.isfile(annotationTSV):
            print('I am unable to locate the annotation file (' + annotationTSV + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
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

def combine_annotations(annotDictList):
    '''
    Combines annotations from each annotDict value in annotDictList to their
    shared key level.
    
    Parameters:
        annotDictList -- a list containing dictionaries with common keys and values as sets
    Returns:
        combinedAnnots -- a dictionary with structure like:
                          {
                              "key1": set(["term1", "term2", ...]),
                              "key2": set(["term3", "term4", ...]),
                              ...
                          }
    '''
    combinedAnnots = {}
    for index, annotDict in enumerate(annotDictList):
        # Validate shared keys in 2nd and subsequent annotation files
        if index != 0:
            if set(combinedAnnots.keys()) != set(annotDict.keys()):
                raise ValueError("Not all keys are shared between annotation files!")
        
        # Iterate through annotations
        for key, annotationSet in annotDict.items():
            # Establish the shared keys that all files must contain
            if index == 0:
                combinedAnnots[key] = annotationSet
            
            # Combine annotations
            else:
                combinedAnnots[key].update(annotationSet)
    return combinedAnnots

def main():
    # User input
    usage = """%(prog)s will receive 2 or more annotation TSV files linking identifiers to annotation terms
    and combine these annotations to the shared identifier level. The output will be a new TSV file
    with the combined annotations.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="annotationTSVs",
                   required=True,
                   nargs="+",
                   help="Specify multiple annotation TSV files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for combined annotations")
    # Opts (behavioural)
    
    # Opts (file format)
    p.add_argument("--delimiter", dest="annotationDelimiter",
                   required=False,
                   help="""Optionally, specify the delimiter used in the 
                   annotation TSV file separating terms; default == '; '""",
                   default="; ")
    p.add_argument("--blank", dest="blankCharacter",
                   required=False,
                   help="""Optionally, specify the character used in the 
                   annotation TSV file to indicate no annotations; default == '0'""",
                   default="0")
    p.add_argument("--hasHeader", dest="hasHeader",
                   required=False,
                   action="store_true",
                   help="Optionally indicate if the annotation TSV file has a header",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse each annotation file
    annotDictList = []
    for annotationTSV in args.annotationTSVs:
        annotDict = parse_annotation_file(annotationTSV, args.annotationDelimiter,
                                          args.blankCharacter, args.hasHeader)
        annotDictList.append(annotDict)
    
    # Combine annotations to orthogroups
    combinedAnnots = combine_annotations(annotDictList)
    
    # Write output file
    with open(args.outputFileName, "w") as outFile:
        for key, annotations in combinedAnnots.items():
            if annotations == set():
                annotations = args.blankCharacter
            else:
                annotations = args.annotationDelimiter.join(annotations)
            outFile.write(f"{key}\t{annotations}\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
