#! python3
# drop_filter_fails.py
# Script to parse a VCF that has been hard-filtered and
# drop any lines containing filter values i.e., which
# are not marked as "PASS"

import os, argparse

# Define functions
def validate_args(args):
    # Validate that all arguments have been provided
    for key, value in vars(args).items():
        if value == None:
            print(key + ' argument was not specified. Fix this and try again.')
            quit()
    # Validate input file locations
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def drop_failed_variants(vcfFile, outputFileName):
    with open(vcfFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        for line in fileIn:
            # Handle comment lines
            if line.startswith("#"): # Just write them all to file
                fileOut.write(line)
                continue
            # Handle variant lines
            sl = line.split("\t")
            if sl[6] == "PASS": # Only write PASS lines
                fileOut.write(line)

def main():
    # User input
    usage = """%(prog)s reads in a VCF that has been hard-filtered and drops any
    variant calls that are not marked as "PASS"
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcf",
        help="Input VCF file name")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the filtered VCF")
    args = p.parse_args()
    validate_args(args)

    # Perform conversion
    drop_failed_variants(args.vcf, args.outputFileName)

if __name__ == "__main__":
    main()
