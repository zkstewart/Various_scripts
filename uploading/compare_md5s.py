#! python3
# compare_md5s.py
# Script to compare two MD5 checksum tables for differences

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.file1):
        raise FileNotFoundError(f"-i1 file '{args.file1}' does not exist or is not a file.")
    args.file1 = os.path.abspath(args.file1)
    
    if not os.path.isfile(args.file2):
        raise FileNotFoundError(f"-i2 file '{args.file2}' does not exist or is not a file.")
    args.file2 = os.path.abspath(args.file2)

def parse_md5_table(tsvFileName, hasHeader, hashIsRight=True):
    '''
    Parse a two-column TSV containing MD5 checksum values,
    returning a dictionary where keys are file basenames
    and value is the checksum as a string.
    '''
    md5Dict = {}
    with open(tsvFileName, "r") as fileIn:
        for line in fileIn:
            try:
                left, right = line.rstrip().split()
            except ValueError:
                raise ValueError(f"MD5 table file '{tsvFileName}' is not a two-column file as expected; " + 
                                 f"offending line is '{line.rstrip()}'")
            
            if hasHeader:
                hasHeader = False
            else:
                if hashIsRight:
                    md5Dict[os.path.basename(left)] = right
                else:
                    md5Dict[os.path.basename(right)] = left
    return md5Dict

def main():
    usage = """%(prog)s accepts two files containing MD5 checkums and compares
    these for any differences. Results are printed to stdout.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i1", dest="file1",
                   required=True,
                   help="First input file")
    p.add_argument("-i2", dest="file2",
                   required=True,
                   help="Second input file")
    # Opts
    p.add_argument("--i1HasHeader", dest="i1HasHeader",
                   required=False,
                   action="store_true",
                   help="First input file has a header; default is no header expected",
                   default=False)
    p.add_argument("--i1HashLeft", dest="i1HashLeft",
                   required=False,
                   action="store_true",
                   help="First input file has the hash in the left column; default is hash in right column",
                   default=False)
    p.add_argument("--i2HasHeader", dest="i2HasHeader",
                   required=False,
                   action="store_true",
                   help="Second input file has a header; default is no header expected",
                   default=False)
    p.add_argument("--i2HashLeft", dest="i2HashLeft",
                   required=False,
                   action="store_true",
                   help="Second input file has the hash in the left column; default is hash in right column",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse MD5 hash tables
    d1 = parse_md5_table(args.file1, args.i1HasHeader, hashIsRight=not args.i1HashLeft)
    d2 = parse_md5_table(args.file2, args.i2HasHeader, hashIsRight=not args.i2HashLeft)
    
    # Perform set comparisons
    d1Keys = set(d1.keys())
    d2Keys = set(d2.keys())
    
    sameKeys = d1Keys.intersection(d2Keys)
    d1Only = d1Keys.difference(d2Keys)
    d2Only = d2Keys.difference(d1Keys)
    
    # Alert user to incomparable files
    if len(sameKeys) == 0:
        print("No files found in common between -i1 and -i2")
    
    # Compare hashes where possible
    else:
        # Note any file incompatibilities
        if len(d1Only) != 0:
            formattedKeys = ", ".join([ f"'{x}'" for x in d1Only ])
            print(f"-i1 indicates files not found in second table, including: {formattedKeys}")
        
        if len(d2Only) != 0:
            formattedKeys = ", ".join([ f"'{x}'" for x in d2Only ])
            print(f"-i2 indicates files not found in first table, including: {formattedKeys}")
        
        # Note any hash differences
        diffHashes = []
        for key in sameKeys:
            hash1 = d1[key]
            hash2 = d2[key]
            if hash1 != hash2:
                diffHashes.append(key)
        
        if len(diffHashes) == 0:
            print("For files common to -i1 and -i2, there were no hash differences found")
        else:
            formattedKeys = ", ".join([ f"'{x}'" for x in diffHashes ])
            print(f"Hash differences found for files listed in -i1 and -i2, including: {formattedKeys}")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
