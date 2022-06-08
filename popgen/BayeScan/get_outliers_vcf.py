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
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gesteFile):
        print('I am unable to locate the GESTE file (' + args.gesteFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.pgdVcfFile):
        print('I am unable to locate the PGDSpider VCF file (' + args.pgdVcfFile + ')')
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

def parse_geste_for_loci_num(gesteFile):
    numLoci = None
    with open(gesteFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("[loci]"):
                numLoci = int(line.rstrip("\r\n ").split("=")[1])
                break
    assert numLoci != None, "GESTE file doesn't have a [loci] line in it...?"
    return numLoci

def parse_pgd_vcf_for_loci(pgdVcfFile):
    '''
    Note that you can use the length of the returned list to tell you how
    many loci were found in this VCF.
    
    Returns:
        vcfLoci -- a list containing sublists with structure like:
                        [
                            [chrom, position, ref, alt],
                            ...
                        ]
    '''
    vcfLoci = []
    with open(pgdVcfFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"): # header line
                continue
            else:
                sl = line.rstrip("\r\n ").split("\t")
                vcfLoci.append([sl[0], sl[1], sl[3], sl[4]])
    return vcfLoci

def write_outlier_vcf(outliers, pgdLoci, vcfFile, outputFileName):
    wroteIndices=[]
    with open(vcfFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        ongoingCount = 1 # this will relate to the index of a BayeScan SNP
        for line in fileIn:
            # Skip header lines
            if line.startswith("#"):
                fileOut.write(line)
            # Handle info lines
            else:
                # Get relevant details
                sl = line.rstrip("\r\n ").split("\t")
                chrom, pos, ref, alt = sl[0], sl[1], sl[3], sl[4]
                
                # Skip lines that aren't used by BayeScan
                found = any(
                    [True
                        for _chrom, _pos, _ref, _alt in pgdLoci
                        if [_chrom, _pos] == [chrom, pos]
                        and sorted(_ref.split(",")) == sorted(ref.split(","))
                        and sorted(_alt.split(",")) == sorted(alt.split(",")) # sometimes PGDSpider sorts it differently...
                    ]
                )
                
                # Skip lines that aren't used by BayeScan
                if not found:
                    continue
                
                # If this line is valid and detected as an outlier, write to file
                if ongoingCount in outliers:
                    fileOut.write(line)
                    wroteIndices.append(ongoingCount)
                
                # Iterate our loci counter
                ongoingCount += 1
    assert wroteIndices == outliers, "Writing the VCF failed! The output file can't be used, sorry!"

def main():
    # User input
    usage = """%(prog)s will curate a VCF file to retain only SNPs that are
    detected as being an outlier by BayeScan. This script is handy since BayeScan
    only gives SNP indices rather than a meaningful identifier.

    Several inputs are needed. You will need 1) the BayeScan outliers text file, and
    2) the original VCF from which the BayeScan GESTE file was created with PGDSpider.
    
    To validate things and ensure everything is working correctly, you also need
    3) a VCF created with PGDSpider using the same settings as the GESTE file was made,
    and 4) the GESTE file itself.
    
    The reason for inputs 3) and 4) is to get the number of loci from the GESTE file
    that BayeScan actually used, and to ensure that we'll extract the correct loci
    from the original VCF.
    
    It's a whole hassle but doing this makes sure your results are valid!
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-b", dest="outliersFile",
                   required=True,
                   help="Input BayeScan outliers text file")
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input VCF file")
    p.add_argument("-pv", dest="pgdVcfFile",
                   required=True,
                   help="Input VCF file that has been run through PGDSpider")
    p.add_argument("-g", dest="gesteFile",
                   required=True,
                   help="Input geste file that BayeScan was run with")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output VCF file filtered to just outlier SNPs")
    args = p.parse_args()
    validate_args(args)
    
    # Get BayeScan outlier SNP indices
    outliers = parse_bayescan_outliers(args.outliersFile)
    
    # Parse PGDSpider VCF
    pgdLoci = parse_pgd_vcf_for_loci(args.pgdVcfFile)
    
    # Parse GESTE file
    gesteLociNum = parse_geste_for_loci_num(args.gesteFile)
    
    # Validate that our loci numbers match
    if len(pgdLoci) != gesteLociNum:
        print("Loci counts differ between the GESTE file ({0}) and PGDSpider VCF file ({1})".format(gesteLociNum, pgdLociNum))
        print("This is an irreconcilable issue. It means your input files aren't matched.")
        quit()
    
    # Write the output file
    write_outlier_vcf(outliers, pgdLoci, args.vcfFile, args.outputFileName)
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
