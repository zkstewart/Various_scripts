#! python3
# merge_reports.py
# Script to receive BSA reports generated from multiple datasets and
# combine them into a single report file for investigation.

import os, argparse

REPORT_COLUMNS = ["#contig", "gene_ID", "coords"]
INSERT_AFTER = "delta_SNPindex"

# Define functions
def validate_args(args):
    # Validate input file locations
    if len(args.reportFiles) < 2:
        print(f"Only received {len(args.reportFiles)} input files; must receive two or more!")
        quit()
    for index, reportFile in enumerate(args.reportFiles):
        if not os.path.isfile(reportFile):
            print(f'I am unable to locate report file #{index+1} ({reportFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate tags parameter
    if not len(args.reportFiles) == len(args.tags):
        print(f"Received {len(args.reportFiles)} input files but {len(args.tags)} tag identifiers")
        print("These values must be equivalent for this program to function correctly")
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def validate_report_files(reportFiles):
    '''
    Validating function that returns True if the files are compatible, and false if not.
    Compatibility is assessed simply by determining whether the file have identical
    header rows or not.
    '''
    isCompatible = True
    lastHeader = None
    for file in reportFiles:
        with open(file, "r") as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n ")
                if lastHeader != None and l != lastHeader:
                    isCompatible = False
                    print(f"Report file incompatibility detected! This header == {l}")
                    print(f"Last header == {lastHeader}")
                lastHeader = l
                break
    return isCompatible

def parse_multiple_report_files(reportFiles, tags):
    '''
    Function which receives multiple report files (which are assumed to be
    pre-validated) and stores data into a dictionary structure so as to allow
    sorted output to be generated later.
    '''
    reportDict = {}
    orderDict = {}
    for index, file in enumerate(reportFiles):
        thisTag = tags[index]
        with open(file, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                # Handle header rows
                if line.startswith("#"):
                    header = sl # all headers should be identical since we assume pre-validation
                    contigIndex = sl.index(REPORT_COLUMNS[0])
                    geneIndex = sl.index(REPORT_COLUMNS[1])
                    coordsIndex = sl.index(REPORT_COLUMNS[2])
                    continue
                # Handle content rows
                else:
                    # Extract relevant details
                    contig = sl[contigIndex]
                    geneID = sl[geneIndex]
                    start, end = map(int, sl[coordsIndex].split("-"))
                    
                    # Index data in our reporting dict
                    reportDict.setdefault(contig, {})
                    reportDict[contig].setdefault(geneID, [])
                    reportDict[contig][geneID].append(sl + [thisTag]) # embed the report's tag into each row
                    
                    # Index data in our ordering/sorting dict
                    orderDict.setdefault(contig, {})
                    orderDict[contig][geneID] = [start, end] # assumed to be the same across files
    
    # Generate sorted output rows
    outputRows = ["\t".join(header + ["source"])] # add ["source"] as header for the tags column
    for contigID in sorted(orderDict.keys()):
        # Figure out the gene order for this contig
        sortedGeneIDs = [
            k
                for k, v in sorted(orderDict[contigID].items(), # don't use item[0] since that's the gene ID
                                   key = lambda item: (item[1][0], item[1][1])) # sort by start,end
        ]
        
        # Store output in order
        for geneID in sortedGeneIDs:
            rows = reportDict[contigID][geneID]
            outputRows += map(lambda x: "\t".join(x), rows)
    
    return outputRows

def main():
    # User input
    usage = """%(prog)s receives two or more BSA report files and merges them into a single cohesive
    file. Tags are set for each input file so it's clear from which report they originate.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-r", dest="reportFiles", nargs="+",
                   required=True,
                   help="Specify multiple BSA report file names name")
    p.add_argument("-t", dest="tags", nargs="+",
                   required=True,
                   help="Specify an identifying tag for each BSA report file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Check that the input report files are compatible with each other
    isCompatible = validate_report_files(args.reportFiles)
    if not isCompatible:
        print("The input files do not all have identical header rows.")
        print("This means some of the files have been modified... hence, the columns will no longer match up.")
        print("Program will exit now; make sure you specify files with identical header rows!")
        quit()
    
    # Parse files into rows structure
    reportRows = parse_multiple_report_files(args.reportFiles, args.tags)
    
    # Write merged report file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(reportRows))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
