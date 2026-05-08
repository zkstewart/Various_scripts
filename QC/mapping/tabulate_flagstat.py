#! python3
# tabulate_flagstat.py
# Simple script to read in .flagstat files
# and tabulate the results for easy understanding

import os, argparse, re

class DirectoryNotFoundError(Exception):
    pass

# Define functions
def validate_args(args):
    # Validate input file locations
    args.flagstatsDir = os.path.abspath(args.flagstatsDir)
    if not os.path.isdir(args.flagstatsDir):
        raise DirectoryNotFoundError(f"-i location '{args.flagstatsDir}' does not exist or is not a directory!")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.isfile(args.outputFileName):
        raise FileExistsError(f"-o file '{args.outputFileName}' already exists and will not be overwritten!")
    
    parentDir = os.path.dirname(args.outputFileName)
    if not os.path.isdir(parentDir):
        raise DirectoryNotFoundError(f"-o file '{args.outputFileName}' cannot be written as its parent directory " + 
                                     f"({parentDir}) does not exist!")

def parse_percentage_line(pctLine):
    '''
    Parameters:
        pctLine -- a string from a flagstat file listing a mapping percentage
    Returns:
        lineInfo -- a string for the type of mapping this line is detailing e.g.,
                    'primary mapped' or 'properly paired'
        pctValue -- a string indicating the percentage that was listed in the file
    '''
    _, relevantBit = pctLine.rstrip("\r\n ").split(" + ") # remove the first number
    relevantBit = relevantBit.split(" ", maxsplit=1)[1] # remove the second number
    lineInfo, lastBit = relevantBit.split(" (")
    pctValue = lastBit.split(" :")[0]
    return lineInfo.replace(" ", "_"), pctValue

def get_flagstats_from_files(flagstatFiles):
    '''
    Parameters:
        flagstatFiles -- a dictionary where keys are sample identifiers and values
                         are strings indicating the location of a flagstat file
    Returns:
        statsTable -- a list of lists with structure like:
                      [
                          [header_row_values],
                          [sample1, statistics],
                          [sample2, statistics],
                          ...
                      ]
    '''
    percentageRegex = re.compile(r"\d{1,3}\.\d{1,2}%")
    
    # Parse the first file to derive the table header
    egSample, egFile = list(flagstatFiles.items())[0] # all files are assumed to have the same format
    
    statsTable=["sample"]
    with open(egFile, "r") as fileIn:
        for line in fileIn:
            if "%" in line:
                lineInfo, pctValue = parse_percentage_line(line)
                statsTable[0] += f"\t{lineInfo}"
    
    # Parse all flagstat files
    for sampleName, file in flagstatFiles.items():
        with open(file, "r") as fileIn:
            statsTable.append([sampleName])
            
            for line in fileIn:
                if "%" in line:
                    lineInfo, pctValue = parse_percentage_line(line)
                    statsTable[-1].append(pctValue)
            statsTable[-1] = "\t".join(statsTable[-1])
    
    return statsTable

## Main
def main():
    # User input
    usage = """%(prog)s will search a provided directory for flagstat files,
    identified by their --suffix, and collate the contained mapping percentages
    to an output TSV table.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="flagstatsDir",
                   required=True,
                   help="Input directory containing flagstat files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the flagstat table")
    p.add_argument("--suffix", dest="fileSuffix",
                   required=False,
                   help="""Optionally, specify the suffix of the flagstat
                   files (default='.flagstat')""",
                   default=".flagstat")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get flagstat files list
    flagFiles = {}
    for file in os.listdir(args.flagstatsDir):
        if file.endswith(args.fileSuffix):
            prefix = os.path.basename(file).rsplit(args.fileSuffix, maxsplit=1)[0].rstrip("._ ")
            if prefix in flagFiles:
                raise ValueError(f"Somehow the file '{file}' has a prefix ({prefix}) which is not unique?")
            
            flagFiles[prefix] = os.path.join(args.flagstatsDir, file)
    
    # Combine flagstat files
    statsTable = get_flagstats_from_files(flagFiles)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(statsTable))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
