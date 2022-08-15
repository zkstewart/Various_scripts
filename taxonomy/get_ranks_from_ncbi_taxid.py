#! python 3
# get_ranks_from_ncbi_taxid.py
# Script to receive a list of NCBI taxonomy IDs
# and convert that into a tabular format that lists
# the various taxonomic ranks above that.

import os, argparse
from ete3 import NCBITaxa

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.tsvFile):
        print('I am unable to locate the input TSV file (' + args.tsvFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if not 1 <= args.columnIndex:
        print("columnIndex only makes sense as an integer >=1 (columns are counted 1-based)")
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print('The specified output file already exists. This program will not allow overwriting.')
        print('Specify a different output file name or move/rename the existing file and try again.')

def is_a_tsv(fileName):
    numRows = None
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if sl == []:
                continue
            
            if numRows != None and len(sl) != numRows:
                return False
            
            numRows = len(sl)
    return True

def parse_taxids_from_tsv_by_colIndex(tsvFileName, columnIndex):
    taxIDs = []
    with open(tsvFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if sl == []:
                continue
            
            colContents = sl[columnIndex-1]
            
            # Skip header rows
            if colContents != ".":
                try:
                    int(colContents)
                except:
                    continue
            
            # Format gaps as None
            colContents = colContents if colContents != "." else None
            taxIDs.append(colContents)
    return taxIDs

def get_all_desired_ranks(NCBITaxa_obj, taxIDsList, desired_ranks):
    resultsDict = {} # lets us skip already-computed results
    resultsList = [] # stores results to return from function
    
    for taxID in taxIDsList:
        if taxID == None:
            resultsList.append([None]*len(desired_ranks))
            continue
        
        if taxID not in resultsDict:
            result = get_desired_ranks(NCBITaxa_obj, taxID, desired_ranks, asTaxID=False)
            resultsDict[taxID] = result
        resultsList.append(resultsDict[taxID])
    return resultsList

def get_desired_ranks(NCBITaxa_obj, taxID, desired_ranks, asTaxID=True):
    '''
    Credit to https://stackoverflow.com/questions/36503042/how-to-get-taxonomic-specific-ids-for-kingdom-phylum-class-order-family-gen
    '''
    lineage = NCBITaxa_obj.get_lineage(taxID)
    lineage2ranks = NCBITaxa_obj.get_rank(lineage)
    ranks2lineage = dict((rank, taxID) for (taxID, rank) in lineage2ranks.items())
    taxIdRanks = [ranks2lineage.get(rank, '<not present>') for rank in desired_ranks]
    if asTaxID is True:
        return taxIdRanks
    else:
        nameTranslationDict = NCBITaxa_obj.get_taxid_translator([t for t in taxIdRanks if t != "<not present>"])
        return [nameTranslationDict[taxID] if taxID in nameTranslationDict else None for taxID in taxIdRanks]

def main():
    usage = """%(prog)s receives a .tsv formatted file and, based upon the specified
    column index which contains NCBI taxonomy IDs or "." for blanks, appends new columns
    to the file containing the ranks from kingdom down to species.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="tsvFile", required=True,
                help="Specify the location of the input TSV file")
    p.add_argument("-o", dest="outputFileName", required=True,
                help="Output file name to be created")
    p.add_argument("--index", dest="columnIndex", required=True, type=int,
                help="Specify the 1-based column index where tax IDs are found")
    args = p.parse_args()
    validate_args(args)
    
    # Load in NCBI taxonomy structure
    RANKS_TO_GET = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ncbi = NCBITaxa()
    
    # Validate that the input file is a TSV
    isTsv = is_a_tsv(args.tsvFile)
    
    # Parse input TSV to retrieve NCBI tax IDs
    taxIDs = parse_taxids_from_tsv_by_colIndex(args.tsvFile, args.columnIndex)
    
    # Get all the taxonomic ranks for each line in the TSV file
    taxRanksList = get_all_desired_ranks(ncbi, taxIDs, RANKS_TO_GET)
    
    # Write to output now
    with open(args.tsvFile, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        firstLine = True
        ongoingCount = 0
        for line in fileIn:
            l = line.rstrip("\r\n")
            
            # Append to header
            if firstLine is True:
                l += "\t{0}".format("\t".join(RANKS_TO_GET))
                firstLine = False
            # Append results to rows
            else:
                rowTaxRanks = [t if t != None else "." for t in taxRanksList[ongoingCount]]
                l += "\t{0}".format("\t".join(rowTaxRanks))
                ongoingCount += 1
            
            fileOut.write(l + "\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
