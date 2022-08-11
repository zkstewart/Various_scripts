#! python3
# tabulate_star_logs.py
# Simple script to read in the Log.final.out files
# produced by STAR after mapping multiple samples and
# tabulates these into a format that's easy to read.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.starOutputsDir):
        print('I am unable to locate the parent directory where STAR subdirectories are (' + args.starOutputsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def get_mapping_stats_from_files(logsFiles):
    '''
    Parameters:
        logsFiles -- a list containing strings that point to STAR output files
                     containing read mapping statistics (Log.final.out)
    Returns:
        statsTable -- a list of lists with structure like:
                      [
                          [header_row_values (\t separated sample names)],
                          [statistic_type_1\tsample_1_result\tsample_2_result...],
                          ...
                      ]
    '''
    RELEVANT_LINES = [
        r"Uniquely mapped reads %", r"% of reads mapped to multiple loci",
        r"% of reads mapped to too many loci", r"% of reads unmapped: too many mismatches",
        r"% of reads unmapped: too short", r"% of chimeric reads"
    ]
    
    # Setup the output table format
    statsTable=[""]
    for rl in RELEVANT_LINES:
        statsTable.append(rl)
    
    # Parse all logs files
    for file in logsFiles:
        # Get the sample name from the directory path
        "Since we have the absolute path, we assume the parent directory is the sample name"
        sampleName = os.path.basename(os.path.dirname(file))
        
        # Parse the file
        with open(file, "r") as fileIn:
            statsTable[0] += "\t{0}".format(sampleName) # on first line, add to header row
            
            for line in fileIn:
                l = line.strip("\r\n ")
                # Skip irrelevant lines
                isRelevant = any([True if l.startswith(rl) else False for rl in RELEVANT_LINES])
                if not isRelevant:
                    continue
                
                # Extract information from line
                statType, stat = l.split("|")
                statType = statType.strip("\t ")
                stat = stat.strip("\t ")
                
                # Get the index for this stat line
                tableIndex = RELEVANT_LINES.index(statType) + 1 # add one since statsTable starts with header row
                
                # All other files
                statsTable[tableIndex] += "\t{0}".format(stat)
    
    return statsTable

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing subdirectories that house
    STAR output files, notably the Log.final.out file, and writes a condensed
    table of read mapping percentages to the output file.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="starOutputsDir", required=True,
        help="Input directory containing subdirectories where STAR was run")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the reads count table")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get Log.final.out files list
    logsFiles = []
    for fileOrDir in os.listdir(args.starOutputsDir):
        if os.path.isdir(os.path.join(args.starOutputsDir, fileOrDir)):
            for subFileOrDir in os.listdir(os.path.join(args.starOutputsDir, fileOrDir)):
                if subFileOrDir == "Log.final.out":
                    logsFiles.append(os.path.join(
                        args.starOutputsDir, fileOrDir, subFileOrDir
                    ))
    
    # Combine ReadsPerGene.out.tab files
    statsTable = get_mapping_stats_from_files(logsFiles)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(statsTable))
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
