#! python3
# get_outliers_vcf.py
# Script to handle BayeScan output results and 
# cut back a VCF to just the outlier SNPs

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.outliersFile):
        print('I am unable to locate the BayeScan outliers file (' + args.outliersFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCFV file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def parse_bayescan_outliers(outliersFile):
    outliers=[]
    with open(outliersFile, "r") as fileIn:
        for line in fileIn:
            l = line.strip(" \r\n").split()
            if l[0] == "prob": # this is the header line, so let's skip it
                continue
            else:
                outliers.append(int(l[0]))
    return outliers

def write_outlier_vcf(outliers, vcfFile, outputFileName):
    with open(vcfFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        ongoingCount = 1 # this will relate to the index of a BayeScan SNP
        for line in fileIn:
            if line.startswith("#"): # header line
                fileOut.write(line)
            else:
                if ongoingCount in outliers:
                    fileOut.write(line)
                ongoingCount += 1

def main():
    # User input
    usage = """%(prog)s will curate a VCF file to retain only SNPs that are
    detected as being an outlier by BayeScan. This script is handy since BayeScan
    only gives SNP indices rather than a meaningful identifier.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-b", dest="outliersFile",
                   required=True,
                   help="Input BayeScan outliers text file")
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input VCF file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output VCF file filtered to just outlier SNPs")
    args = p.parse_args()
    validate_args(args)
    
    # Get BayeScan outlier SNP indices
    outliers = parse_bayescan_outliers(args.outliersFile)
    
    # Parse and write the output in one
    write_outlier_vcf(outliers, args.vcfFile, args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
