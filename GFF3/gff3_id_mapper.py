#! python3
# gff3_id_mapping.py
# Uses the ZS_GFF3IO.GFF3 Class to parse a GFF3 file and
# allow flexible logic for mapping features to one or more
# attribute values.

import os, argparse, re, sys
from collections import OrderedDict

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
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

def gff3_id_mapper(forEach, map, to, gff3Obj):
    '''
    Provides a way to form a structured query of the GFF3 object to
    retrieve details in a list.
    
    Parameters:
        forEach -- a string indicating the type of GFF3 feature to check.
                   Whatever value is given here, must be a key within
                   the gff3Obj.types dictionary.
        map -- a string indicating the "primary key" that will be
               the left-most column of the mapping list.
        to -- a list containing strings of any valid keys that can be
              obtained from the GFF3 feature that is being iterated over
              with the forEach parameter.
        gff3Obj -- a ZS_GFF3IO.GFF3 instance.
    Returns:
        mapping -- a list containing sublists of mappings.
    '''
    # Conform and validate inputs
    forEach = gff3Obj._make_feature_case_appropriate(forEach)
    assert forEach in gff3Obj.types, \
        "'{0}' forEach value does not exist within the GFF3; \
        valid options include '{1}'".format(forEach, list(gff3Obj.types.keys()))
    
    # Perform the query operation
    mapping = []
    for feature in gff3Obj.types[forEach]:
        assert map in feature.__dict__, \
            "'{0}' map key does not exist within one or more features; \
            valid keys include '{1}'; check that letters have correct \
            upper/lower casing for {2}".format(map, list(feature.__dict__.keys()), feature)
        mappingDetail = feature.__dict__[map]
        
        toDetails = []
        for toKey in to:
            assert toKey in feature.__dict__, \
                "'{0}' to key does not exist within one or more features; \
                valid keys include '{1}'; check that letters have correct \
                upper/lower casing for {2}".format(toKey, list(feature.__dict__.keys()), feature)
            toDetails.append(feature.__dict__[toKey])
        
        mapping.append([mappingDetail] + toDetails)
    return mapping

if __name__ == "__main__":
    # User input
    usage = """%(prog)s accepts a GFF3 file and provides a simple way to format
    a statement that will extract relevant information mappings in a TSV output file.
    
    The statement consists of 3 parameters, namely being -forEach __ -map __ -to __.
    Logically, this should form a kind of sentence where we are saying e.g., "for each mRNA,
    map ID to Parent". In this example, the output file would give columns linking the mRNA 
    ID to its parent gene ID. Multiple -to __ keys are accepted, and each additional key will
    give an extra column in the output file. Separate -to __ keys with a space.
    
    Valid keys are any attribute key=value pairing in column 9 of the GFF3. Alternatively,
    the first 8 columns have their own key names which are: contig, source, ___ (specify in
    -forEach), start, end, score, strand, frame. Keys are case-sensitive!
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    ## Logical structure parameters
    p.add_argument("-forEach", dest="forEach", required=True,
        help="""
            The value here should relate to any feature type in column 3 of the GFF3;
            examples include gene or mRNA.
        """)
    p.add_argument("-map", dest="map", required=True,
        help="""
            The value here should relate to any valid key that will be the left-most column
            of the output and serve as the anchor for other details.
        """)
    p.add_argument("-to", dest="to", required=True, nargs="+",
        help="One or more valid keys separated with a space.")
    ## File I/O parameters
    p.add_argument("-g", dest="gff3File", required=False,
                help="Specify the location of the GFF3 file")
    p.add_argument("-o", dest="outputFileName", required=False,
                help="Specify the name and location for the output file")
    
    args = p.parse_args()
    validate_args(args)
    
    # Load GFF3
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, False) # non-strict parsing
    
    # Obtain data linkings
    mapping = gff3_id_mapper(args.forEach, args.map, args.to, gff3Obj)
    
    # Write output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("#{0}\t{1}\n".format(args.map, "\t".join(args.to))) # header
        for row in mapping:
            fileOut.write("{0}\n".format("\t".join(row)))
    
    print("Program completed successfully!")
