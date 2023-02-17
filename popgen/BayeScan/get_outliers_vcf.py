#! python3
# get_outliers_vcf.py
# Script to handle BayeScan output results and 
# cut back a VCF to just the outlier SNPs

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.outliersFstFile):
        print('I am unable to locate the BayeScan outliers FST file (' + args.outliersFstFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
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

def parse_bayescan_outliers(outliersFile):
    outliers=set()
    with open(outliersFile, "r") as fileIn:
        for line in fileIn:
            l = line.strip(" \r\n").split()
            if l[0] == "prob": # this is the header line, so let's skip it
                continue
            else:
                outliers.add(int(l[0]))
    return outliers

def parse_bayescan_verif(verifFile):
    verifDict = {}
    with open(verifFile, "r") as fileIn:
        for line in fileIn:
            if not line.startswith("Pop. "):
                continue
            else:
                sl = line.rstrip("\r\n ").split()
                pop = sl[1]
                locus = int(sl[3])
                nums = list(map(int, sl[5:]))
                
                verifDict.setdefault(locus, {})
                verifDict[locus][pop] = nums
    return verifDict

def parse_pops_metadata(metadataFile):
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if sl == []:
                continue
            else:
                sample, pop = sl
                metadataDict.setdefault(pop, set())
                metadataDict[pop].add(sample)
    return metadataDict

def validate_bayescan_outliers(outliers, verifDict, metadataDict):
    # Check that our outliers don't exceed the range of verifDict
    for outlierIndex in outliers:
        if outlierIndex not in verifDict:
            raise ValueError(
                f"{outlierIndex} outlier index not found in verif file!"
            )
    
    # Check that verifDict only has pops found in metadataDict
    metadataPops = set(metadataDict.keys())
    for index, popAlleles in verifDict.items():
        if not set(popAlleles.keys()) == metadataPops:
            raise ValueError(
                f"Pops in verif file ({set(popAlleles.keys())}) don't match metadata ({metadataPops})"
            )

def get_pop_from_key(sampleID, metadataDict):
    for key, value in metadataDict.items():
        if sampleID in value:
            return key
    return None

def get_vcf_lines_of_outliers(vcfFile, verifDict, metadataDict, outliers):
    '''
    Parameters:
        vcfFile -- as string indicating the location of the VCF file used for
                   BayeScan outlier prediction
        verifDict -- a dictionary containing the parsed contents of the BayeScan
                     verif text file
        metadataDict -- a dictionary containing the parsed contents of the metadata
                        text file used for GESTE file creation
        outliers -- a set containing BayeScan indices of outlier SNPs
    Returns:
        linesToWrite -- a set containing 0-based indices of VCF content lines
                        that should be output since they are 1) outliers, and 
                        2) are index-adjusted to account for changes during GESTE
                        file creation from the original VCF.
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
                    if verifCounter in outliers:
                        linesToWrite.add(vcfCounter) 
                    verifCounter += 1
            
                vcfCounter += 1
    return linesToWrite

def main():
    # User input
    usage = """%(prog)s will curate a VCF file to retain only SNPs that are
    detected as being an outlier by BayeScan. This script is handy since BayeScan
    only gives SNP indices rather than a meaningful identifier.
    
    Several inputs are needed:
    1) the original VCF file prior to conversion to GESTE file,
    2) the BayeScan FST file subset to only significant outlier loci,
    3) the BayeScan verif.txt file, and
    4) the pops metadata file used for GESTE file creation
    
    We need to do this because, when converting to GESTE file for BayeScan input,
    some variants can be silently dropped. Hence, your original VCF may be different
    to the input and your variant indices will differ. It's a hassle but this script
    will make sure your results are valid!
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-vcf", dest="vcfFile",
                   required=True,
                   help="Input VCF file")
    p.add_argument("-f", dest="outliersFstFile",
                   required=True,
                   help="Input BayeScan outliers FST file")
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
    
    # Get BayeScan outlier SNP indices
    outliers = parse_bayescan_outliers(args.outliersFstFile)
    
    # Parse BayeScan verif file
    verifDict = parse_bayescan_verif(args.verifFile)
    
    # Parse metadata file
    metadataDict = parse_pops_metadata(args.metadataFile)
    
    # Validate that outliers, verif, and metadata match
    validate_bayescan_outliers(outliers, verifDict, metadataDict)
    
    # Figure out which lines we are going to write as output
    linesToWrite = get_vcf_lines_of_outliers(args.vcfFile, verifDict, metadataDict, outliers)
    
    # Write the output file
    with open(args.vcfFile, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        ongoingCount = 1
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
