#! python3
# add_qtl_info_to_proximity_reports.py
# Script to parse outputs from R/qtl and add relevant information
# into the SNP and gene proximity reports.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    for suffix in [".gene_proximity.tsv", ".snp_proximity.tsv"]:
        if not os.path.isfile(args.proximityReportPrefix + suffix):
            print(f'I am unable to locate one of the proximity report files ({args.proximityReportPrefix + suffix})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.scanoneFile == None and args.scanTwoFile == None and args.windowsFile == None:
        print("At least one of the scanone, scantwo, or windows file(s) must be provided!")
        print("Correct this issue and try again.")
        quit()
    if args.scanoneFile != None:
        if not os.path.isfile(args.scanoneFile):
            print(f'I am unable to locate the scanone file ({args.scanoneFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.scantwoFile != None:
        if not os.path.isfile(args.scantwoFile):
            print(f'I am unable to locate the scantwo file ({args.scantwoFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.windowsFile != None:
        if not os.path.isfile(args.windowsFile):
            print(f'I am unable to locate the windows file ({args.windowsFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate numeric inputs
    if args.kbRadius < 1:
        print("kbRadius must be a positive integer (greater than zero)")
    if args.pvalue < 0:
        print("kbRadius must be an integer greater than or equal to zero")
    if args.lod < 0:
        print("los must be an integer greater than or equal to zero")
    # Handle file overwrites
    for suffix in [".gene_proximity.tsv", ".snp_proximity.tsv"]:
        if os.path.isfile(args.outputFilePrefix + suffix):
            print(f'"{args.outputFilePrefix + suffix}" already exists. Delete/move/rename this file and run the program again.')
            quit()

def parse_scanone_file(scanoneFile, pvalueFilter=0.05, lodFilter=10.0):
    '''
    Expects a TSV format produced by me (ZKS) from R/qtl. The file has a single
    header row with split by tab looking like:
    
    ["", "chr", "pos", "lod", "pval"]
    
    The first column itself is formatted like:
    
        S{chr}_{bp}
    
    Where {chr} == the value in the next column over, and {bp} is the numeric position
    in the chromosome measured in basepairs.
    
    If the file hasn't had a P-value filter enforced, this function will do that here.
    
    Parameters:
        scanoneFile -- a string indicating the file location of the scanone output
                       from R/qtl.
        pvalueFilter -- a float indicating the P-value significance to enforce.
    Returns:
        scanoneDict -- a dictionary with structure like:
            {
                "chr1": [bp1, bp2] # bp == nucleotide bp as int value
            }
    '''
    scanoneDict = {}
    
    firstLine = True
    with open(scanoneFile, "r") as fileIn:
        for line in fileIn:
            if firstLine == True:
                firstLine = False
                continue
            else:
                # Get relevant information from this line
                chr_bp, chr, pos, lod, pval = line.rstrip("\r\n ").split("\t")
                bp = int(chr_bp.split("_")[-1])
                lod = float(lod)
                pval = float(pval)
                
                # Skip if pval or lod don't pass cutoffs
                if pval > pvalueFilter or lod < lodFilter:
                    continue
                
                # Store results otherwise
                scanoneDict.setdefault(chr, [])
                scanoneDict[chr].append(bp)
    
    return scanoneDict

def parse_scantwo_file(scantwoFile, pvalueFilter=0.05, lodFilter=10.0):
    '''
    Expects a TSV format produced by me (ZKS) from R/qtl. The file has a single
    header row with split by tab looking like:
    
    ["", "chr1", "chr2", "pos1f", "pos2f", "lod.full", "pval", ... "chr1_pos", "chr2_pos"]
    (intervening columns don't matter to us)
    
    If the file hasn't had a P-value filter enforced, this function will do that here.
    
    Parameters:
        scantwoFile -- a string indicating the file location of the scantwo output
                       from R/qtl.
        pvalueFilter -- a float indicating the P-value significance to enforce.
    Returns:
        scantwoDict -- a dictionary with structure like:
            {
                n: [ # == integer indexing this pair of QTL values
                    "chr1", "pos1", "chr2", "pos2"
                ],
                ...
            }
    '''
    scantwoDict = {}
    
    firstLine = True
    ongoingCount = 0
    with open(scantwoFile, "r") as fileIn:
        for line in fileIn:
            if firstLine == True:
                firstLine = False
                continue
            else:
                # Get relevant information from this line
                sl = line.rstrip("\r\n ").split("\t")
                chr1, chr2, lod, pval, pos1, pos2 = sl[1], sl[2], sl[5], sl[6], sl[-2], sl[-1]
                pos1, pos2 = int(pos1), int(pos2)
                lod = float(lod)
                pval = float(pval)
                
                # Skip if pval or lod don't pass cutoffs
                if pval > pvalueFilter or lod < lodFilter:
                    continue
                
                # Store results otherwise
                scantwoDict[ongoingCount] = [
                    chr1, pos1, chr2, pos2
                ]
                ongoingCount += 1
    
    return scantwoDict

def parse_windows_file(windowsFile, cols=["contig", "start", "end"]):
    '''
    Expects a TSV or CSV file containing, at minimum, three columns with a single
    header row with split by tab containing the columns indicated by the cols value.
    
    Parameters:
        windowsFile -- a string indicating the file location of a TSV or CSV file
                       containing at least the contig ID, start, and end positions
                       as numbers
        cols -- a list containing three values indicating the header column values
                for the contig, start, and end positions (in that order!)
    Returns:
        windowsDict -- a dictionary with structure like:
            {
                "chr1": [[bp1, bp2], [ ... ], ... ], # bp == nucleotide bp as int value
                "chr2": [ [ ... ], ... ],
                ...
            }
    '''
    windowsDict = {}
    
    firstLine = True
    with open(windowsFile, "r") as fileIn:
        for line in fileIn:
            # Parse TSV or CSV
            if "\t" in line:
                sl = line.rstrip("\r\n ").split("\t")
            else:
                sl = line.rstrip("\r\n ").split(",")
            
            # Handle header line
            if firstLine == True:
                contigIndex = sl.index(cols[0])
                startIndex = sl.index(cols[1])
                endIndex = sl.index(cols[2])
                
                firstLine = False
                continue
            # Handle content lines
            else:
                # Get relevant information from this line
                contig = sl[contigIndex]
                start = int(sl[startIndex])
                end = int(sl[endIndex])
                
                # Store results
                windowsDict.setdefault(contig, [])
                windowsDict[contig].append([start, end])
    
    return windowsDict

def main():
    # User input
    usage = """%(prog)s receives the proximity report files generated by
    snp_proximity_report.py and adds R/qtl results into the file.
    Specifically, it will receive one or both of the scanone and scantwo
    outputs and note any SNPs or genes containing SNPs which are within
    a N-kb window centred on the predicted QTL region. The size of the window
    can be specified here.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-p", dest="proximityReportPrefix",
                   required=True,
                   help="""Specify the prefix to the two proximity reports e.g.,
                   '/path/to/qtlseq_comparison' which should then match to two files
                   with '.gene_proximity.tsv' and '.snp_proximity.tsv' suffixes""")
    p.add_argument("-o", dest="outputFilePrefix",
                   required=True,
                   help="Specify the output file prefix (will create two new prox report files)")
    # Optional
    p.add_argument("-s1", dest="scanoneFile",
                   required=False,
                   help="Input scanone file produced by Rqtl")
    p.add_argument("-s2", dest="scantwoFile",
                   required=False,
                   help="Input scantwo file produced by Rqtl")
    p.add_argument("-w", dest="windowsFile",
                   required=False,
                   help="Input windows file containing at least 3 columns for contig, start, end")
    p.add_argument("--windowColumns", dest="windowColumns", nargs="+",
                   required=False,
                   help="""Optionally, specify the names of the three columns to parse from the windows file
                   (default == ["contig", "start", "end"]""",
                   default=["contig", "start", "end"])
    p.add_argument("--kbRadius", dest="kbRadius", type=int,
                   required=False,
                   help="Optionally, specify the radius of the window in kb (default==100)",
                   default=100)
    p.add_argument("--pvalue", dest="pvalue", type=float,
                   required=False,
                   help="Optionally, specify the P-value to enforce for QTL discovery (default==0.05)",
                   default=0.05)
    p.add_argument("--lod", dest="lod", type=float,
                   required=False,
                   help="Optionally, specify the minimum LOD to accept from a result (default==10.0)",
                   default=10)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get our proximity report file names
    snpProxFile = os.path.join(args.proximityReportPrefix + ".snp_proximity.tsv")
    geneProxFile = os.path.join(args.proximityReportPrefix + ".gene_proximity.tsv")
    
    # Parse the scanone file (if relevant)
    if args.scanoneFile != None:
        scanoneDict = parse_scanone_file(args.scanoneFile, args.pvalue, args.lod)
    else:
        scanoneDict = {}
    
    # Parse the scantwo file (if relevant)
    if args.scantwoFile != None:
        scantwoDict = parse_scantwo_file(args.scantwoFile, args.pvalue, args.lod)
    else:
        scantwoDict = {}
    
    # Parse the windows file (if relevant)
    if args.scantwoFile != None:
        windowsDict = parse_windows_file(args.windowsFile, args.windowColumns)
    else:
        windowsDict = {}
    
    # Write new SNP proximity file
    with open(snpProxFile, "r") as fileIn, open(args.outputFilePrefix + ".snp_proximity.tsv", "w") as fileOut:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            
            # Update header
            if l.startswith("#"):
                l += "\tscanone_window\tscantwo_window\ttsnp_window"
                fileOut.write(f"{l}\n")
            # Update content lines
            else:
                contig, pos = l.split("\t")[0:2]
                pos = int(pos)
                
                # Check to see if this position corresponds to a radius around a scanone or scantwo site
                isScanOne = False
                if contig.upper() in scanoneDict:
                    qtlPosList = scanoneDict[contig.upper()]
                    for qtlPos in qtlPosList:
                        lowerBoundary = qtlPos - (args.kbRadius * 1000)
                        upperBoundary = qtlPos + (args.kbRadius * 1000)
                        
                        if lowerBoundary <= pos <= upperBoundary:
                            isScanOne = True
                            break
                
                scanTwoNum = []
                for qtlNum, qtlValues in scantwoDict.items():
                    chr1, pos1, chr2, pos2 = qtlValues
                    
                    if chr1 == contig.upper():
                        lowerBoundary = pos1 - (args.kbRadius * 1000)
                        upperBoundary = pos1 + (args.kbRadius * 1000)
                        
                        if lowerBoundary <= pos <= upperBoundary:
                            scanTwoNum.append(qtlNum)
                    
                    elif chr2 == contig.upper():
                        lowerBoundary = pos2 - (args.kbRadius * 1000)
                        upperBoundary = pos2 + (args.kbRadius * 1000)
                        
                        if lowerBoundary <= pos <= upperBoundary:
                            scanTwoNum.append(qtlNum)
                
                # Check to see if it's contained within a window
                isInWindow = False
                if contig in windowsDict:
                    thisContigWindows = windowsDict[contig]
                    for windowStart, windowEnd in thisContigWindows:
                        if windowStart <= pos <= windowEnd:
                            isInWindow = True
                            break
                
                # Write line to file
                l += "\t{0}\t{1}\t{2}".format(
                    "y" if isScanOne is True else "n",
                    "; ".join([f"qtl{n}" for n in scanTwoNum]) if scanTwoNum != [] else ".",
                    "y" if isInWindow is True else "n"
                )
                fileOut.write(f"{l}\n")
    
    # Write new gene proximity file
    with open(geneProxFile, "r") as fileIn, open(args.outputFilePrefix + ".gene_proximity.tsv", "w") as fileOut:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            sl = l.split("\t")
            
            # Update header
            if l.startswith("#"):
                coordsIndex = sl.index("coords") # for retrieving coords from content lines
                contigIndex = sl.index("contig") # likewise
                l += "\tscanone_window\tscantwo_window\tsnp_window"
                fileOut.write(f"{l}\n")
            # Update content lines
            else:
                geneID = sl[0]
                start, end = list(map(int, sl[coordsIndex].split("-")))
                contig = sl[contigIndex]
                
                # Check to see if this gene overlaps any QTL regions
                isScanOne = False
                if contig.upper() in scanoneDict:
                    qtlPosList = scanoneDict[contig.upper()]
                    for qtlPos in qtlPosList:
                        lowerBoundary = qtlPos - (args.kbRadius * 1000)
                        upperBoundary = qtlPos + (args.kbRadius * 1000)
                        
                        if start <= upperBoundary and end >= lowerBoundary: # if it overlaps
                            isScanOne = True
                            break
                
                scanTwoNum = []
                for qtlNum, qtlValues in scantwoDict.items():
                    chr1, pos1, chr2, pos2 = qtlValues
                    
                    if chr1 == contig.upper():
                        lowerBoundary = pos1 - (args.kbRadius * 1000)
                        upperBoundary = pos1 + (args.kbRadius * 1000)
                        
                        if start <= upperBoundary and end >= lowerBoundary: # if it overlaps
                            scanTwoNum.append(qtlNum)
                    
                    elif chr2 == contig.upper():
                        lowerBoundary = pos2 - (args.kbRadius * 1000)
                        upperBoundary = pos2 + (args.kbRadius * 1000)
                        
                        if start <= upperBoundary and end >= lowerBoundary: # if it overlaps
                            scanTwoNum.append(qtlNum)
                
                # Check to see if it's contained within a window
                isInWindow = False
                if contig in windowsDict:
                    thisContigWindows = windowsDict[contig]
                    for windowStart, windowEnd in thisContigWindows:
                        
                        if (windowStart <= start <= windowEnd) or (windowStart <= end <= windowEnd):
                            isInWindow = True
                            break
                
                # Write line to file
                l += "\t{0}\t{1}\t{2}".format(
                    "y" if isScanOne is True else "n",
                    "; ".join([f"qtl{n}" for n in scanTwoNum]) if scanTwoNum != [] else ".",
                    "y" if isInWindow is True else "n"
                )
                fileOut.write(f"{l}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
