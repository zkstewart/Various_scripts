#! python3
# extract_vcf.py
# Simply extracts variants from a VCF and produces
# the modified output.

import os, argparse, gzip
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.tsvFile):
        print('I am unable to locate the TSV file (' + args.tsvFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()
    # Validate numeric inputs
    if len(args.columns) != 2:
        print("columns must be given exactly two values")
        quit()

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def parse_tsv_file(tsvFile, headerColumns):
    warnedAlready = False
    filterDict = {}
    firstLine = True
    with open(tsvFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            
            # Validate header on first line
            if firstLine == True:
                for i in range(len(headerColumns)):
                    header = headerColumns[i]
                    if header not in sl:
                        print(f"'{header}' wasn't found in the TSV file!")
                        print("Make sure your TSV file and --columns values match, then try again.")
                        quit()
                    elif i == 0:
                        contigIndex = sl.index(header)
                    else:
                        posIndex = sl.index(header)
                firstLine = False
            # Handle content lines
            else:
                contig, position = sl[contigIndex], sl[posIndex]
                filterDict.setdefault(contig, set())
                
                # Warn exactly once if we detect a duplicate
                if position in filterDict[contig]:
                    if warnedAlready == False:
                        print(f"WARNING: Position '{position}' already found in TSV for contig '{contig}'" +
                            "; if duplicates are expected, you can ignore this.")
                        warnedAlready = True
                
                # Store in dictionary for function output
                filterDict[contig].add(position)
    
    return filterDict

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF and a TSV with two columns indicating
    the contig ID and position for extraction from the VCF. Produces a new VCF
    with the same header to the original, but only with the lines specified.
    The TSV is expected to contain exactly one header row with the header values
    specified in this script.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input VCF file for filtering")
    p.add_argument("-t", dest="tsvFile",
                   required=True,
                   help="Input TSV file with contig and position columns")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the filtered VCF")
    ## Optional
    p.add_argument("--columns", dest="columns", nargs="+",
                   required=False,
                   help="""Optionally, provide the two column headers if they differ
                   to the default; order MUST BE contig column then position column
                   (default=="contig position")""",
                   default=["contig", "position"])
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse TSV file for contig and positions
    filterDict = parse_tsv_file(args.tsvFile, args.columns)
    
    # Extract relevant entries from VCF and write to file
    foundDict = {}
    with open_vcf_file(args.vcfFile) as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            # Handle header rows
            if line.startswith("#"):
                fileOut.write(line)
            
            # Handle content rows:
            else:
                sl = line.rstrip("\r\n").split("\t")
                chrom, pos = sl[0:2]
                
                # Handle hits
                if chrom in filterDict and pos in filterDict[chrom]:
                    fileOut.write(line)
                    foundDict.setdefault(chrom, set())
                    foundDict[chrom].add(pos)
                # Handle misses
                else: 
                    continue
    
    # Print any relevant warning messages
    for contig, posSet in filterDict.items():
        if contig not in foundDict:
            print(f"WARNING: TSV indicates variants in contig '{contig}'; none were found")
        else:
            foundSet = foundDict[contig]
            missingSet = posSet.difference(foundSet)
            if len(missingSet) != 0:
                print(f"WARNING: {len(missingSet)} variants on contig '{contig}' were not found in the VCF")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
