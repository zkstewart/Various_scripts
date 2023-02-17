#! python3
# get_outliers_vcf.py
# Script to handle BayeScan output results and 
# cut back a VCF to just the outlier SNPs

import os, argparse
from get_outliers_vcf import * # we'll overwrite validate_args() and main() locally

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.verifFile):
        print('I am unable to locate the BayeScan verif file (' + args.verifFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print('I am unable to locate the pops metadata file (' + args.metadataFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def validate_bayescan(verifDict, metadataDict):
    # Check that verifDict only has pops found in metadataDict
    metadataPops = set(metadataDict.keys())
    for index, popAlleles in verifDict.items():
        if not set(popAlleles.keys()) == metadataPops:
            raise ValueError(
                f"Pops in verif file ({set(popAlleles.keys())}) don't match metadata ({metadataPops})"
            )

def get_vcf_lines(vcfFile, verifDict, metadataDict):
    '''
    Parameters:
        vcfFile -- as string indicating the location of the VCF file used for
                   BayeScan outlier prediction
        verifDict -- a dictionary containing the parsed contents of the BayeScan
                     verif text file
        metadataDict -- a dictionary containing the parsed contents of the metadata
                        text file used for GESTE file creation
    Returns:
        linesToWrite -- a set containing 0-based indices of VCF content lines
                        that should be output since they are 1) are index-adjusted
                        to account for changes during GESTE file creation from
                        the original VCF.
    '''
    
    pops = set(metadataDict.keys())
    
    linesToWrite = set() # will hold onto the line numbers we should write as outliers
    
    verifCounter = 1
    vcfCounter = 0
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            # Handle header line
            if line.startswith("#CHROM"):
                sampleIDs = sl[9:]
            # Skip comment lines
            elif line.startswith("#"):
                continue
            # Handle body lines
            else:
                # Extract details from line
                chrom, pos, id, ref, alt, \
                    qual, filt, info, \
                    format = sl[0:9]
                samples = sl[9:]
                gtIndex = format.split(":").index("GT")
                
                # Tabulate allele frequency for each population
                alleleDict = {
                    pop: [0 for x in range(0, len(alt.split(",")) + 1)]
                        for pop in pops
                }
                
                for i in range(len(samples)):
                    # Extract details for this sample's VCF data
                    sampleID = sampleIDs[i]
                    samplePop = get_pop_from_key(sampleID, metadataDict)
                    
                    sampleGT = samples[i].split(":")[gtIndex]
                    sampleAlleles = sampleGT.split("/")
                    
                    # Skip uncalled genotypes
                    if sampleGT == "./." or sampleGT == ".":
                        continue
                    
                    # Store result in allele dict
                    for allele in sampleAlleles:
                        alleleDict[samplePop][int(allele)] += 1

                # Check to see if this lines up with the verif file
                verifAlleles = verifDict[verifCounter]
                
                if all([ set(verifAlleles[pop]) == set(alleleDict[pop])
                            for pop in pops
                ]):
                    linesToWrite.add(vcfCounter)
                    verifCounter += 1
                
                vcfCounter += 1
    return linesToWrite

def main():
    # User input
    usage = """%(prog)s will curate a VCF file to retain only SNPs that are
    part of the GESTE file that BayeScan has received as input.
    
    Several inputs are needed:
    1) the original VCF file prior to conversion to GESTE file,
    2) the BayeScan verif.txt file, and
    3) the pops metadata file used for GESTE file creation.
    
    We need to do this because, when converting to GESTE file for BayeScan input,
    some variants can be silently dropped. Hence, your original VCF may be different
    to the input and your variant indices will differ. It's a hassle but this script
    will make sure your results are valid!
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-vcf", dest="vcfFile",
                   required=True,
                   help="Input VCF file")
    p.add_argument("-verif", dest="verifFile",
                   required=True,
                   help="Input BayeScan verif text file")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Input pops metadata file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output VCF file filtered to just outlier SNPs")
    args = p.parse_args()
    validate_args(args)
    
    # Parse BayeScan verif file
    verifDict = parse_bayescan_verif(args.verifFile)
    
    # Parse metadata file
    metadataDict = parse_pops_metadata(args.metadataFile)
    
    # Validate that outliers, verif, and metadata match
    validate_bayescan(verifDict, metadataDict)
    
    # Figure out which lines we are going to write as output
    linesToWrite = get_vcf_lines(args.vcfFile, verifDict, metadataDict)
    
    # Write the output file
    with open(args.vcfFile, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        ongoingCount = 0
        for line in fileIn:
            if line.startswith("#"):
                fileOut.write(line)
            else:
                if ongoingCount in linesToWrite:
                    fileOut.write(line)
                ongoingCount += 1
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
