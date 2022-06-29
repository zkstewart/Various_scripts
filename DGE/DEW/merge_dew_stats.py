#! python3
# merge_dew_stats.py
# A simple python program to combine DEW output .stats files into
# a table that's amenable to my way of doing DGE analysis

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.statsDir):
        print('I am unable to locate the directory containing .stats files (' + args.statsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def get_stats_from_files(statsFiles):
    fileCount = 0
    statsTable=["gene_id"]
    
    for file in statsFiles:
        lineCount = 0
        
        with open(file, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                
                # First line of a new file
                if lineCount == 0:
                    if "_vs_" in file:
                        species = file.split("_vs_")[1].split(".")[0] # relevant file name exists between "_vs_" and "."
                    else:
                        species = file.split(".")[0] # assumes relevant file name precedes all "." characters
                    statsTable[0] += "\t{0}".format(species)
                
                # First file being checked
                if fileCount == 0:
                    statsTable.append("\t".join([sl[0], sl[2]])) # get rid of the second (gene length) and last (?? unsure) columns
                # All other files
                else:
                    statsTable[lineCount + 1] += "\t{0}".format(sl[2])
                
                lineCount += 1
        fileCount += 1
    
    return statsTable

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing DEW .stats files, and writes a condensed
    table of read counts to the specified output file.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="statsDir", required=True,
        help="Input directory containing .stats files")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the reads count table")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get .stats files list
    statsFiles = [os.path.join(args.statsDir, file) for file in os.listdir(args.statsDir) if file.endswith(".stats")]
    
    # Combine .stats files
    statsTable = get_stats_from_files(statsFiles)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(statsTable[:-1])) # get rid of the last asterisk row
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
