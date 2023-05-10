#! python3
# vcf_overlap_gff3_features.py
# Simply extracts variants from a VCF that overlap a specific
# kind of GFF3 feature and writes the filtered VCF

import os, argparse, gzip, sys
from contextlib import contextmanager

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 3 dirs up is where we find dependencies

from Function_packages import ZS_GFF3IO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the TSV file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
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
    usage = """%(prog)s reads in a VCF and a GFF3 file and, depending on
    the type of feature in the GFF3 you specify, extracts variants that overlap
    those features and outputs a subset VCF file.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input VCF file for filtering")
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the filtered VCF")
    ## Optional
    p.add_argument("--feature", dest="featureType",
                   required=False,
                   help="""Optionally, indicate which feature type overlaps should be
                   searched for within; default=="CDS", but "exon" is a good
                   option (note case sensitivity)""",
                   default="CDS")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3 and validate featureType parameter
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    if not args.featureType in gff3Obj.types:
        print(f"{args.featureType} is not recognised as a feature in your GFF3")
        print("Fix your argument and try again.")
    
    # NCLS indexing for finding overlaps
    gff3Obj.create_ncls_index(typeToIndex=args.featureType)
    
    # Extract relevant entries from VCF and write to file
    numFound = 0
    with open_vcf_file(args.vcfFile) as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            # Handle header rows
            if line.startswith("#"):
                fileOut.write(line)
            
            # Handle content rows:
            else:
                sl = line.rstrip("\r\n").split("\t")
                chrom, pos = sl[0:2]
                pos = int(pos)
                
                # Locate GFF3 overlaps
                matches = gff3Obj.ncls_finder(pos, pos, "contig", chrom)
                
                # Handle hits
                if len(matches) > 0:
                    fileOut.write(line)
                    numFound += 1
                # Handle misses
                else: 
                    continue
    
    # Print statistical information
    if numFound == 0:
        print("WARNING: No VCF lines were output!")
        print(f"That's because no variants overlapped your feature type '{args.featureType}'")
        print("Something may be wrong...")
    else:
        print(f"INFO: {numFound} variants were found which overlapped '{args.featureType}' features")
        print("Hopefully this is approximately what you expected!")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
