#! python3
# windows_vcf.py
# Script to subset a VCF to just the variants within provided
# window region(s).

import os, argparse, re, gzip
from contextlib import contextmanager

WINDOW_REGEX = re.compile(r"^(.+?):(\d+?)-(\d+)$")

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate windows format
    for window in args.windows:
        if WINDOW_REGEX.match(window) == None:
            print(f"'{window}' isn't formatted appropriately. Fix this and try again.")
            quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def filter_vcf_by_windows(vcf, windowsDict):
    '''
    Filter to remove variants from a ZS_VCFIO.VCF object which are
    NOT contained within a given window region.
    
    Parameters:
        vcf -- a ZS_VCFIO.VCF object
        windowsDict -- a dictionary with structure like:
                       {
                           'contigID1': [
                               [x1, y1],
                               [x2, y2],
                               ...
                           ],
                           'contigID2': [
                               ...
                           ],
                           ...
                       }
    '''
    
    
    posDrop = {}
    contigDrop = []
    for contigID, contigDict in vcf.items():
        # Fail if contig not in a window
        if not contigID in windowsDict:
            contigDrop.append(contigID)
            continue
        
        # Check if position should be failed
        for pos, posDict in contigDict.items():
            # Check if in window
            isInWindow = any(
                [ int(pos) >= start and int(pos) <= end
                  for start, end in windowsDict[contigID]
                ])
            
            # Fail if it's not
            if not isInWindow:
                posDrop.setdefault(contigID, set())
                posDrop[contigID].add(pos)
    
    # Eliminate contigs entirely if they had no windows
    for contigID in contigDrop:
        vcf.del_contig(contigID)
    
    # Eliminate variants if they failed the filter
    for contigID, posSet in posDrop.items():
        for pos in posSet:
            vcf.del_variant(contigID, pos)

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF alongside one or more window
    regions to OBTAIN from the VCF, outputting a subsetted VCF containing
    only variants found within the given window(s).
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input VCF file for filtering")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the filtered SNPs")
    p.add_argument("-w", dest="windows", nargs="+",
                   required=True,
                   help="""Indicate window(s) to get by providing one or more strings
                   with format 'CHROM:start-end'. Start and end should both be numeric
                   values, and CHROM should correspond to a contig ID. Multiple windows
                   can be provided for the same CHROM, but you should still provide
                   these separately.""")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse windows argument into dictionary structure
    windowsDict = {}
    for window in args.windows:
        chrom, start, end = WINDOW_REGEX.match(window).groups()
        windowsDict.setdefault(chrom, [])
        windowsDict[chrom].append([int(start), int(end)])
    
    # Filter VCF by window ranges to output file
    outputLines = 0
    with open_vcf_file(args.vcfFile) as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            # Handle header rows
            if line.startswith("#CHROM"): # insert this header comment just before the #CHROM line
                headerComment = f"##windows_vcf.py windows={args.windows}"
                fileOut.write(f"{headerComment}\n")
            if line.startswith("#"):
                fileOut.write(line)
            
            # Handle content rows:
            else:
                sl = line.rstrip("\r\n").split("\t")
                contigID, pos = sl[0:2]
                pos = int(pos)
                
                # Fail if contig not in a window
                if not contigID in windowsDict:
                    continue
                
                # Check if this position is in a window
                isInWindow = any([
                    pos >= start and pos <= end
                    for start, end in windowsDict[contigID]
                ])
                
                # Write to output if it is
                if isInWindow:
                    fileOut.write(line)
                    outputLines += 1
    
    # Write filtered output (if relevant)
    if outputLines == 0:
        print("In the process of filtering this VCF, we ended up with 0 SNPs remaining!")
        print("You should check to make sure your parameters make sense...")
        print("(Since there's no SNPs left, I'm going to clean up the empty header file now)")
        os.unlink(args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
