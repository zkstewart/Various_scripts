#! python3
# tabulate_md5s.py
# Simple script to read in .md5 files
# and tabulate the results for easy understanding

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.md5sDir):
        raise FileNotFoundError(f"Input directory '{args.md5sDir}' does not exist or is not a directory.")
    args.md5sDir = os.path.abspath(args.md5sDir)
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        raise FileExistsError(f"Output file '{args.outputFileName}' already exists. Please choose a different " + 
                              "name or move or delete the existing file.")

def parse_md5_file(file, md5Suffix=".md5"):
    '''
    Parameters:
        file -- a string containing the path to a .md5 file
                containing read md5 checksums.
        md5Suffix -- a string containing the suffix of the md5 file
                     to be used if we need to derive the original
                     file's name (default=='.md5')
    Returns:
        checksum -- a string containing the checksum
        name -- a string containing the name of the file
    '''
    checksum = None
    with open(file, "r") as fileIn:
        for line in fileIn:
            sl = line.strip("\r\n\t ").split()
            
            # Handle checksum only line
            if len(sl) == 1:
                checksum = sl[0]
                name = os.path.basename(file.rsplit(md5Suffix, maxsplit=1)[0]) # derive original name from md5 file
            # Handle checksum and name line
            elif len(sl) == 2:
                checksum = sl[0]
                name = os.path.basename(sl[1])
                if not os.path.isfile(sl[1]):
                    print(f"WARNING: File '{sl[1]}' specified in md5sum file '{file}' does not exist?")
            # Handle unexpected line format
            else:
                raise ValueError(f"Unexpected line format in file '{file}': {line.strip()}")
            break # we only need the first valid line
    
    if checksum != None:
        return checksum, name
    else:
        raise ValueError(f"The file '{file}' appears to be empty or is not in the expected MD5 format")

def get_md5_from_files(md5Files, md5Suffix=".md5"):
    '''
    Parameters:
        md5Files -- a list containing strings that point to .md5
                    files containing read md5 checksums.
        md5Suffix -- a string containing the suffix of the md5 file
                     to be passed along to the parse_md5_file function
                     (default=='.md5')
    Returns:
        md5sTable -- a list of strings, where each string
                     is a tab-separated line containing
                     the sample name and its checksums.
    '''
    md5sTable=["file\tchecksum"]
    for file in md5Files:
        checksum, name = parse_md5_file(file, md5Suffix=md5Suffix)
        md5sTable.append(f"{name}\t{checksum}")
    return md5sTable

## Main
def main():
    #### USER INPUT SECTION
    usage = """%(prog)s accepts a directory containing .md5 files and
    produces a tabulated output of checksums for each file.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="md5sDir",
                   required=True,
                   help="Input directory containing .md5 files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the tabulated results")
    # Opts
    p.add_argument("--suffix", dest="filesSuffix",
                   required=False,
                   help="""Optionally specify the suffix for the files
                   to be processed (default==.md5)""",
                   default=".md5")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get .md5 files list
    md5Files = []
    for file in os.listdir(args.md5sDir):
        if file.endswith(args.filesSuffix):
            md5Files.append(os.path.join(args.md5sDir, file))
    
    # Combine .md5 files
    md5sTable = get_md5_from_files(md5Files, md5Suffix=args.filesSuffix)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(md5sTable) + "\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
