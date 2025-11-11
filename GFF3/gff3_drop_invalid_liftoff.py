#! python3
# gff3_drop_invalid_liftoff.py
# Script to receive a GFF3 and a released nucleotide sequence
# file and will attempt to correct any GFF3 features where
# the coordinates don't match the real feature.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_GFF3IO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def main():
    # User input
    usage = """%(prog)s ...
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="gff3File",
                   required=True,
                   help="Specify the first GFF3 file name.")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output merged GFF3 file name.")
    # Opts (minor behavioural mods)
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3
    GFF3_obj = ZS_GFF3IO.GFF3(args.gff3File,
                              strict_parse = not args.relaxedParsing,
                              fix_duplicated_ids = False) # leave this for gff3_reformat.py
    
    # Write to file
    liftoffExclusions = 0
    childrenScolded = 0
    nonExclusions = 0
    noAttribute = 0
    with open(args.outputFileName, "w") as fileOut:
        for parentType in GFF3_obj.parentTypes:
            for parentFeature in GFF3_obj.types[parentType]:
                # Skip over invalid features
                if hasattr(parentFeature, "valid_ORFs"):
                    if parentFeature.valid_ORFs == "0":
                        liftoffExclusions += 1
                        continue
                else:
                    noAttribute += 1
                
                # Delete any child mRNAs which lack a valid ORF
                badChildren = []
                for childFeature in parentFeature.children:
                    if childFeature.valid_ORF == "False":
                        childrenScolded += 1
                        badChildren.append(childFeature.ID)
                for badChildID in badChildren:
                    parentFeature.del_child(badChildID, "mRNA")
                
                # Write current feature to file
                ZS_GFF3IO.GFF3._recursively_write_feature_details(parentFeature, fileOut)
                nonExclusions += 1
    
    # Emit filtering statistic
    print("# gff3_drop_invalid_liftoff.py output statistics:")
    print(f"# > {liftoffExclusions} genes were dropped due to having no valid ORFs.")
    print(f"# > {nonExclusions} genes were kept.")
    print(f"# > {childrenScolded} mRNAs were dropped due to having no valid ORF")
    print(f"# > {noAttribute} genes were skipped due to not having a valid_ORFs attribute.")
    
    # All done!
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
