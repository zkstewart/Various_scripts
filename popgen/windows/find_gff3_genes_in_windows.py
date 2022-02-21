#! python3
# find_gff3_genes_in_windows.py
# Script to take in a GFF3 file and a windows file
# (tab or comma delimited, at least three columns i.e., 
# chromosome : start : end). Numbers should be 1-based
# just like in a GFF3. Other columns will be ignored, and
# a header column can optionally be indicated.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.windowsFile):
        print('I am unable to locate the windows text file (' + args.windowsFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3):
        print('I am unable to locate the GFF3 file (' + args.gff3 + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def get_gene_boundaries(gff3File):
    geneRegions = {} # this will store ALL potential gene regions as simple geneID : range() pairs
    outputRegions = {} # this will store ONLY gene regions that we've identified as containing an mRNA
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            l = line.rstrip("\r\n").split("\t")
            chrom = l[0]
            annotType = l[2]
            start = int(l[3])
            end = int(l[4])
            details = l[8].strip("\"").split(';')
            detail_dict = {}
            for i in range(len(details)):
                if details[i] == '':
                    continue
                split_details = details[i].split('=', maxsplit=1)
                detail_dict[split_details[0]] = split_details[1]
            
            # Set up storage structure at the start of each chromosome
            if chrom not in geneRegions:
                geneRegions[chrom] = {}
                outputRegions[chrom] = {}
            
            # Work out the chrom -> gene -> mRNA pairings 
            if annotType == "gene":
                geneRegions[chrom][detail_dict["ID"]] = [start, end] # store potential gene boundaries
            elif annotType == "mRNA":
                geneID = detail_dict["Parent"]
                outputRegions[chrom][geneID] = geneRegions[chrom][geneID] # store a gene boundary the function should output
            
    return outputRegions

def get_window_boundaries(windowsFile, hasHeader=False, isCsv=False):
    windowBoundaries = {}
    skippedHeader = False
    with open(windowsFile, "r") as fileIn:
        for line in fileIn:
            # Handle header lines
            if hasHeader and not skippedHeader:
                skippedHeader = True
                continue
            
            # Extract relevant details
            if isCsv:
                l = line.rstrip("\r\n").split(",")
            else:
                l = line.rstrip("\r\n").split("\t")
            chrom = l[0]
            start = int(l[1])
            end = int(l[2])
            
            # Store in dictionary
            if chrom not in windowBoundaries:
                windowBoundaries[chrom] = []
            windowBoundaries[chrom].append([start, end])
            
    return windowBoundaries

def get_gene_window_overlaps(geneBoundaries, windowBoundaries):
    genes = []
    for chrom, wBoundaries in windowBoundaries.items():
        for wBoundary in wBoundaries:
            for geneID, gBoundary in geneBoundaries[chrom].items():
                gStart, gEnd = gBoundary
                wStart, wEnd = wBoundary
                if gStart <= wEnd and wStart <= gEnd: # head-to-tail comparison for overlap of two ranges
                    genes.append(geneID)
    return genes

## Main
def main():
    # User input
    usage = """%(prog)s reads in a GFF3 and a tab/comma-delimited file containing
    at least 3 columns (chromosome : start : end) indicating genomic windows
    that you wish to locate genes within. Start and end numbers should be 1-based
    like a GFF3, and any columns beyond the first three will be ignored. If the
    windows file has a header row, specify --header so it can be ignored. If it's
    comma delimited, specify --csv so the program knows it's not tab-delimited.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-w", dest="windowsFile", required=True,
        help="Input tab-delimited files indicating the genomic windows")
    p.add_argument("-g", dest="gff3", required=True,
        help="Input GFF3 file")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output text file name")
    ## Optional
    p.add_argument("--header", dest="hasHeader", required=False, action="store_true",
        help="Optionally indicate that the windows file has a header row to be ignored",
        default=False)
    p.add_argument("--csv", dest="isCsv", required=False, action="store_true",
        help="Optionally indicate that the windows file is comma delimited; assumed to be tab otherwise",
        default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Simple parse of GFF3 for CDS regions and their corresponding gene ID
    geneBoundaries = get_gene_boundaries(args.gff3)
    
    # Parse windows file to get windows of interest
    windowBoundaries = get_window_boundaries(args.windowsFile, args.hasHeader, args.isCsv)
    
    # Find genes which overlap windows
    genes = get_gene_window_overlaps(geneBoundaries, windowBoundaries)
    
    # Produce output file0
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(genes))
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
