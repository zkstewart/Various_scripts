#! python3
# concat_star_readsPerGene.py
# Combines STAR output ReadsPerGene.out.tab files into a table that's amenable
# to DGE analysis

import os, argparse

class DirectoryNotFoundError(Exception):
    pass

# Define functions
def validate_args(args):
    # Validate input file locations
    args.starOutputsDir = os.path.abspath(args.starOutputsDir)
    if not os.path.isdir(args.starOutputsDir):
        raise DirectoryNotFoundError(f"-i location '{args.starOutputsDir}' does not exist or is not a directory!")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.isfile(args.outputFileName):
        raise FileExistsError(f"-o file '{args.outputFileName}' already exists and will not be overwritten!")
    
    parentDir = os.path.dirname(args.outputFileName)
    if not os.path.isdir(parentDir):
        raise DirectoryNotFoundError(f"-o file '{args.outputFileName}' cannot be written as its parent directory " + 
                                     f"({parentDir}) does not exist!")

def get_counts_from_files(countsFiles):
    '''
    Parameters:
        countsFiles -- a list containing strings that point to STAR output files
                       containing reads counts (ReadsPerGene.out.tab)
    Returns:
        statsTable -- a list of lists with structure like:
                      [
                          [header_row_values (\t separated sample names)],
                          [gene_id_1\tsample_1_count\tsample_2_count...],
                          ...
                      ]
    '''
    IRRELEVANT_LINES = [
        "N_unmapped", "N_multimapping",
        "N_noFeature", "N_ambiguous"
    ]
    fileCount = 0
    statsTable=["gene_id"]
    
    for file in countsFiles:
        # Get the sample name from the directory path
        "Since we have the absolute path, we assume the parent directory is the sample name"
        sampleName = os.path.basename(os.path.dirname(file))
        
        # Parse the file
        lineCount = 0
        with open(file, "r") as fileIn:
            statsTable[0] += "\t{0}".format(sampleName) # on first line, add to header row
            
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                if sl == []:
                    continue
                geneName, count, _, _ = sl # drop last 2 stranded counts, just work with unstranded count
                
                # Skip irrelevant lines
                if geneName in IRRELEVANT_LINES:
                    continue
                
                # First file being checked
                if fileCount == 0:
                    statsTable.append("\t".join([geneName, count]))
                    
                # All other files
                else:
                    statsTable[lineCount + 1] += "\t{0}".format(count)
                
                lineCount += 1
        fileCount += 1
    
    return statsTable

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing subdirectories that house
    STAR output files, notably the ReadsPerGene.out.tab file, and writes a condensed
    table of read counts to the specified output location.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="starOutputsDir",
                   required=True,
                   help="Input directory containing subdirectories where STAR was run")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the reads count table")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get ReadsPerGene.out.tab files list
    countsFiles = []
    for fileOrDir in os.listdir(args.starOutputsDir):
        if os.path.isdir(os.path.join(args.starOutputsDir, fileOrDir)):
            for subFileOrDir in os.listdir(os.path.join(args.starOutputsDir, fileOrDir)):
                if subFileOrDir == "ReadsPerGene.out.tab":
                    countsFiles.append(os.path.join(
                        args.starOutputsDir, fileOrDir, subFileOrDir
                    ))
    
    # Combine ReadsPerGene.out.tab files
    statsTable = get_counts_from_files(countsFiles)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(statsTable))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
