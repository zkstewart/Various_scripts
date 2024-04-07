#! python3
# rename_gff3_from_best_orthologs.py
# Script to rename the gene models in a GFF3 based on the results of
# orthofinder_best_orthologs.py, and specifically its manual_curation.tsv
# file which contains the information necessary to relate a gene to its
# best SOI ortholog.

# Load packages
import os, argparse, sys
from collections import Counter

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_GFF3IO

def validate_args(args):
    # Validate input locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.orthologsTSV):
        print('I am unable to locate the manual_curation.tsv[-like] file (' + args.orthologsTSV + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def extract_tsv_header(tsvFileName):
    with open(tsvFileName, 'r') as fileIn:
        header = fileIn.readline().strip().split('\t')
    return header[1:]

def parse_orthologs_tsv(fileName, headerOfInterest):
    '''
    Parameters:
        fileName -- a string indicating the location of the file to parse.
        headerOfInterest -- a string indicating the header value to pair with
                            sequence IDs.
    Returns:
        orthologsDict -- a dictionary with structure like:
                         {
                             "headerSeqID1": "sequenceID1",
                             "headerSeqID2": "sequenceID2",
                             ...
                         }
    '''
    orthologsDict = {}
    isFirstLine = True
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if isFirstLine:
                ofInterestIndex = sl.index(headerOfInterest)
                isFirstLine = False
            else:
                assert sl[ofInterestIndex] not in orthologsDict, \
                    f"Species of interest sequence '{sl[ofInterestIndex]}' occurs more than once!"
                orthologsDict[sl[ofInterestIndex]] = sl[0]
    return orthologsDict

def main():
    ### USER INPUT
    usage = """%(prog)s will parse the curation TSV resulting from orthofinder_best_orthologs.py
    and use it to rename the gene models in a GFF3 file. The TSV file should have a format wherein
    the first column contains the sequence ID to have its name modeled after, and another column
    (the header of which is specified with -s) containing your species whose GFF3 is to be renamed.
    It will write a new GFF3 as output. Note that any genes in your GFF3 which do not feature in
    the TSV file will be left unchanged. MOST IMPORTANTLY, note that it's expected that the header
    to be renamed contains parent feature IDs e.g., "gene" IDs and not "mRNA" IDs.
    """
    
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Specify the GFF3 file to be renamed.")
    p.add_argument("-t", dest="orthologsTSV",
                   required=True,
                   help="Specify the ortholog TSV file to use for renaming.")
    p.add_argument("-s", dest="soi",
                   required=True,
                   help="Column name of the species you are interested in renaming.")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output output, renamed GFF3 file name.")
    # Opts
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Extract header from TSV file
    header = extract_tsv_header(args.orthologsTSV)
    
    # Make sure SOI is in the header
    if args.soi not in header:
        print(f"SOI '{args.soi}' not found in OrthoGroups TSV file header.")
        print('Make sure your SOI is within the file header and try again; header shown below.')
        print(header)
        quit()
    
    # Parse the orthologs TSV file
    orthologsDict = parse_orthologs_tsv(args.orthologsTSV, args.soi)
    
    # Parse the GFF3 file
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File,
                             strict_parse = not args.relaxedParsing)
    
    # Iteratively update GFF3 file
    "Note for future: a lot of this code should be embedded logic in the GFF3 class"
    geneCounts = {}
    for parentType in gff3Obj.parentTypes:
        for parentFeature in gff3Obj.types[parentType]:
            # Skip if this gene is not being renamed
            if parentFeature.ID not in orthologsDict:
                continue
            
            # Generate a new ID for this gene
            idForRenaming = orthologsDict[parentFeature.ID]
            geneCounts.setdefault(idForRenaming, 0)
            geneCounts[idForRenaming] += 1
            
            newParentID = f"{idForRenaming}.copy{geneCounts[idForRenaming]}"
            parentFeature.update_id(parentFeature.ID, newParentID, gff3Obj)
            
            # Delete any transcript children features
            if hasattr(parentFeature, "transcript"):
                transcriptIDs = [ transcriptFeature.ID for transcriptFeature in parentFeature.transcript ]
                for transcriptID in transcriptIDs:
                    parentFeature.del_child(transcriptID, "transcript")
            
            # Update the children features more directly
            for index, childFeature in enumerate(parentFeature.children):
                newChildID = f"{newParentID}.{childFeature.type}{index+1}"
                childFeature.update_id(childFeature.ID, newChildID, gff3Obj)
                
                if hasattr(childFeature, "CDS"):
                    for cdsIndex, cdsFeature in enumerate(childFeature.CDS):
                        newCDSID = f"{newChildID}.CDS{cdsIndex+1}"
                        cdsFeature.update_id(cdsFeature.ID, newCDSID, gff3Obj)
                
                if hasattr(childFeature, "exon"):
                    for exonIndex, exonFeature in enumerate(childFeature.exon):
                        newExonID = f"{newChildID}.exon{exonIndex+1}"
                        exonFeature.update_id(exonFeature.ID, newExonID, gff3Obj)
    
    # Write the updated GFF3 file
    gff3Obj.write(args.outputFileName)
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
