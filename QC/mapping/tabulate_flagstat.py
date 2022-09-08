#! python3
# tabulate_flagstat.py
# Simple script to read in .flagstat files
# and tabulate the results for easy understanding

import os, argparse, re

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.flagstatsDir):
        print('I am unable to locate the parent directory where .flagstat files are (' + args.flagstatsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def get_flagstats_from_files(flagstatFiles):
    '''
    Parameters:
        flagstatFiles -- a list containing strings that point to .flagstat
                         files containing read mapping statistics
    Returns:
        statsTable -- a list of lists with structure like:
                      [
                          [header_row_values],
                          [sample1, statistics],
                          [sample2, statistics],
                          ...
                      ]
    '''
    percentageRegex = re.compile(r"\d{1,2}\.\d{1,2}%")
    
    # Setup the output table format
    statsTable=["sample\tmapped_percentage"]
    
    # Parse all logs files
    for file in flagstatFiles:
        # Get the sample name from the directory path
        "Since we have the absolute path, we assume the parent directory is the sample name"
        sampleName = os.path.basename(file).split(".")[0]
        
        # Parse the file
        with open(file, "r") as fileIn:
            statsTable.append([sampleName])
            
            for line in fileIn:
                l = line.strip("\r\n ")
                # Skip irrelevant lines
                isRelevant = "%" in l
                if not isRelevant:
                    continue
                
                # Extract information from line
                percentage = percentageRegex.findall(l)[0]
                
                # All other files
                statsTable[-1].append(percentage)
            statsTable[-1] = "\t".join(statsTable[-1])
    
    return statsTable

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing .flagstat files 
    and writes a condensed table of read mapping percentages to the output file.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="flagstatsDir", required=True,
        help="Input directory containing .flagstat files")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the reads count table")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get .flagstat files list
    flagFiles = []
    for file in os.listdir(args.flagstatsDir):
        if file.endswith(".flagstat"):
            flagFiles.append(os.path.join(args.flagstatsDir, file))
    
    # Combine .flagstat files
    statsTable = get_flagstats_from_files(flagFiles)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(statsTable))
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
