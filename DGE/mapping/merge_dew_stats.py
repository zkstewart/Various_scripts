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
    # Parse files into dictionary of counts
    statsDict = {}
    uniqueGeneIDs = set()
    uniqueSpecies = []
    
    for file in statsFiles:
        with open(file, "r") as fileIn:
            firstLine = True
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                modelID, _, count, _ = sl
                if modelID == "*":
                    continue
                uniqueGeneIDs.add(modelID)
                
                # First line of a new file
                if firstLine is True:
                    if "_vs_" in file:
                        species = file.split("_vs_")[1].split(".")[0] # relevant file name exists between "_vs_" and "."
                    else:
                        species = file.split(".")[0] # assumes relevant file name precedes all "." characters
                    uniqueSpecies.append(species)
                    statsDict[species] = {}
                    firstLine = False
                
                # All other lines
                statsDict[species][modelID] = count
    uniqueGeneIDs = list(uniqueGeneIDs)
    
    # Untangle dictionary into a list with proper formatting
    statsTable=["gene_id\t{0}".format("\t".join(uniqueSpecies))]
    for geneID in uniqueGeneIDs:
        statsTable.append(geneID)
    for species in uniqueSpecies:
        for i in range(1, len(statsTable)):
            geneID = uniqueGeneIDs[i-1] # our two lists are out of sync by 1
            try:
                statsTable[i] += "\t{0}".format(statsDict[species][geneID])
            except:
                print("DEW encountered a known bug")
                print("If you check the .stats files, they don't all have the same length")
                print("Realistically, this means the DEW counts cannot be trusted. Sorry.")
                quit()
    
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
