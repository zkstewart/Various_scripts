#! python3
# tair_addin.py
# Helper script to add TAIR BLAST results into the bulked segregant analysis
# report table.

import os, argparse

REPORT_COLUMNS = ["gene_ID", "mRNA_ID"]
INSERT_AFTER = "gene_name"

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.tairFile):
        print(f'I am unable to locate the TAIR annotation file ({args.tairFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.reportFile):
        print(f'I am unable to locate the report table file ({args.reportFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def main():
    # User input
    usage = """%(prog)s is a simple helper script for adding TAIR annotation
    details into the BSA report file.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="tairFile",
                   required=True,
                   help="Specify the TAIR annotation file name")
    p.add_argument("-r", dest="reportFile",
                   required=True,
                   help="Specify the report file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse SNP report to find gene IDs
    tairDict = {}
    with open(args.reportFile, 'r') as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            
            # Handle header line
            if line.startswith("#"):
                geneIndex = sl.index(REPORT_COLUMNS[0])
                mrnaIndex = sl.index(REPORT_COLUMNS[1])
                continue
            
            # Handle content lines
            geneID = sl[geneIndex]
            # mrnaID = sl[mrnaIndex] # shouldn't be needed
            tairDict[geneID] = None
    
    # Parse TAIR annotation table for relevant gene annotations
    with open(args.tairFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            
            # Skip header lines
            if line.startswith("#"):
                continue
            
            # Extract relevant information from this line
            mrnaID, urID, urName, urEvalue, tairID, tairName, tairEvalue, urUrl, tairUrl  = sl
            
            geneID = mrnaID.split(".")[0]
            if geneID not in tairDict:
                continue
            
            # Index data
            tairDict[geneID] = {
                "TAIR ID": tairID,
                "TAIR gene name": tairName,
                "TAIR URL": tairUrl,
                "UniRef ID": urID,
                "UniRef URL": urUrl
            }
    
    # Update annotation table
    with open(args.reportFile, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Update header line
            if line.startswith("#"):
                # Format insertion for this line
                insertIndex = sl.index(INSERT_AFTER)
                newHeaderBits = [
                    "TAIR_gene_name", "TAIR_hit_URL"
                ]
                
                # Write out updated line
                newLine = sl[0:insertIndex+1] + newHeaderBits + sl[insertIndex+1:]
                fileOut.write("\t".join(newLine) + "\n")
                continue
            
            # Update content lines
            else:
                # Get relevant data
                geneID = sl[geneIndex] # this variable still exists in scope
                tairAnnotation = tairDict[geneID]
                
                # Format insertion for this line
                newBits = [
                    tairAnnotation["TAIR gene name"],
                    tairAnnotation["TAIR URL"]
                ]
                
                # Write out updated line
                newLine = sl[0:insertIndex+1] + newBits + sl[insertIndex+1:]
                fileOut.write("\t".join(newLine) + "\n")

if __name__ == "__main__":
    main()
