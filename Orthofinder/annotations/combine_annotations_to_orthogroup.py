#! python3
# combine_annotations_to_orthogroup
# Script to take in an annotation file linking gene IDs to annotation terms
# and combines these annotations to the orthogroup level that they were
# clustered within.

# Load packages
import os, argparse

def validate_args(args):
    # Validate input locations
    if not os.path.isfile(args.orthogroupsTSV):
        print('I am unable to locate the Orthogroups file (' + args.orthogroupsTSV + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.annotationTSV):
        print('I am unable to locate the annotation file (' + args.annotationTSV + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.exists(args.outputFileName):
        print(f"{args.outputFileName} already exists, and I will not overwrite it.")
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

def parse_orthogroups_tsv(orthogroupsFile, speciesID):
    '''
    Simple function to parse a Orthogroups.tsv file from OrthoFinder into
    a dictionary structure from which all relevant data can be obtained for
    this script.
    
    Parameters:
        orthogroupsFile -- a string representing the path to the Orthogroups.tsv file
        speciesID -- a string representing the species ID to be used to extract gene IDs
    Returns:
        orthoDict -- a dictionary with structure like:
                     {
                         "OG000001": set(["gene1", "gene2", ...]), # all genes in orthogroup
                         "OG000002": set(["gene3", "gene4", ...]), # regardless of species
                         ...
                     }
    '''
    orthoDict = {}
    firstLine = True
    with open(orthogroupsFile , "r") as orthoFile:
        for line in orthoFile:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # this won't support sequence IDs containing quotation marks
            # Handle the header line
            if firstLine == True:
                assert speciesID in sl, \
                    f"Species ID {speciesID} not found in Orthogroups header: {sl}"
                speciesIndex = sl.index(speciesID)
                firstLine = False
            # Handle all other lines
            else:
                orthogroupName, genes = sl[0], sl[speciesIndex]
                if genes != "":
                    genes = set(genes.split(", "))
                else:
                    genes = set()
                orthoDict[orthogroupName] = genes
    return orthoDict

def aggregate_annotations(setDict, annotDict):
    '''
    Aggregates annotations for sequences in annotDict up to the level of the setDict's
    keys (orthogroups).
    
    Parameters:
        setDict -- a dictionary with orthogroup names as keys and sets of gene IDs as values
        annotDict -- a dictionary with gene IDs as keys and sets of annotation terms as values
    Returns:
        aggregatedAnnots -- a dictionary with structure like:
                            {
                                "OG000001": set(["term1", "term2", ...]),
                                "OG000002": set(["term3", "term4", ...]),
                                ...
                            }
    '''
    aggregatedAnnots = {}
    for key, geneIDs in setDict.items():
        aggregatedAnnots[key] = set()
        for geneID in geneIDs:
            try:
                aggregatedAnnots[key].update(annotDict[geneID])
            except:
                raise ValueError(f"Gene ID {geneID} not found in annotation TSV file")
    return aggregatedAnnots

def main():
    # User input
    usage = """%(prog)s will receive an annotation TSV linking sequence identifiers to annotation terms
    alongside an OrthoFinder Orthogroups.tsv file. Annotations will be aggregated to the orthogroup level
    and output to a new TSV file.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="orthogroupsTSV",
                   required=True,
                   help="Specify the orthogroup TSV file")
    p.add_argument("-a", dest="annotationTSV",
                   required=True,
                   help="Specify the annotation TSV file")
    p.add_argument("-s", dest="speciesID",
                   required=True,
                   help="Specify the species ID as found in the orthogroup header")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for combined annotations")
    # Opts
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
    
    # Parse Orthogroups.tsv
    orthoDict = parse_orthogroups_tsv(args.orthogroupsTSV, args.speciesID)
    
    # Parse annotations
    annotDict = parse_annotation_file(args.annotationTSV, args.annotationDelimiter,
                                      args.blankCharacter, args.hasHeader)
    
    # Combine annotations to orthogroups
    aggregatedAnnots = aggregate_annotations(orthoDict, annotDict)
    
    # Write output file
    with open(args.outputFileName, "w") as outFile:
        for orthogroup, annotations in aggregatedAnnots.items():
            if annotations == set():
                annotations = args.blankCharacter
            else:
                annotations = args.annotationDelimiter.join(annotations)
            outFile.write(f"{orthogroup}\t{annotations}\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
