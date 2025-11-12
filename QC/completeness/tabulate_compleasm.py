#! python3
# tabulate_compleasm.py
# Simple script to read in compleasm summary files
# and tabulate the results for easy understanding

import os, argparse, re
import pandas as pd

class DirectoryNotFoundError(Exception):
    pass

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.compleasmDir):
        raise DirectoryNotFoundError(f"-i directory '{args.compleasmDir}' could not be found or is not a directory")
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        raise FileExistsError(f"-o file '{args.outputFileName}' already exists and will not be overwritten")

def compleasm_summaries_to_pandas_df(summaryFiles):
    '''
    Parameters:
        summaryFiles -- a list containing strings that point to 'summary.txt' compleasm files
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
    summaryDict = {}
    for file in summaryFiles:
        # Get the sample name from the directory path
        "Since we have the absolute path, we assume the parent directory is the sample name"
        sampleName = os.path.basename(os.path.dirname(file))
        
        # Split off irrelevant suffix (if it exists)
        if sampleName.endswith("compleasm"):
            sampleName = sampleName.rsplit("compleasm", maxsplit=1)[0].strip("_.")
        
        # Parse the file
        summaryDict[sampleName] = {}
        with open(file, "r") as fileIn:            
            for line in fileIn:
                l = line.strip("\r\n ")
                if l.startswith("#"):
                    continue
                
                # Extract information from line
                key, value = l.split(":")
                percentage = value.split("%")[0]
                
                # Store in dict
                summaryDict[sampleName][key] = percentage
    
    # Convert dict to df
    statsDf = pd.DataFrame.from_dict(summaryDict, orient="index")
    return statsDf

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing compleasm output directories
    (which each contain a 'summary.txt' file) and writes a condensed table of summary
    statistics to the output file
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="compleasmDir",
                   required=True,
                   help="Input directory containing compleasm output subdirectories")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the reads count table")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get summary files list
    summaryFiles = []
    for subdir in [ os.path.join(args.compleasmDir, d) for d in os.listdir(args.compleasmDir) if os.path.isdir(os.path.join(args.compleasmDir, d)) ]:
        summaryFile = os.path.join(subdir, "summary.txt")
        if os.path.isfile(summaryFile):
            summaryFiles.append(summaryFile)
    
    # Combine summary files
    statsDf = compleasm_summaries_to_pandas_df(summaryFiles)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        statsDf.to_csv(fileOut, sep="\t")
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
