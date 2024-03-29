#! python3
# sample_vcf_to_snp_index.py
# Script to receive a Freebayes prediction VCF (performed per-sample)
# and convert it into a file that mimics the SNP index files produced
# by QTL-seq. The main point is to calculate delta SNP index which can
# be used to filter SNPs down to just the important ones for a BSA
# analysis.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print('I am unable to locate the metadata file (' + args.metadataFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

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

def vcf_to_snp_index(vcfFile, metadataDict):
    '''
    This function will parse a VCF file and return a dictionary which contains
    data that corresponds to the QTL-seq SNP index format (minus the p95 and
    p99 calculations).
    
    Note1: SNP-index (at a position) = Count of alternate base / Count of reads aligned.
    
    Note2: According to the QTL-seq publication, they filter such that "positions with
    read depth < 7 in both the bulks and SNP-index < 0.3 in either of the bulks were
    filtered out, and SNPs with homozygous alleles in both the bulks were used for
    delta SNP-index calculation."
    
    Note3: According to the QTL-seq publication, "Only SNP positions with delta SNP-index = -1
    (i.e., the allele called in high trait value-bulk was the same as that of the resistant
    parent while contrastingly different in low trait value-bulk)"
    
    Parameters:
        vcfFile -- a string indicating the file location of the VCF file
        metadataDict -- a dictionary with structure like:
                        {
                            'sampleID1': 'bulk1',
                            'sampleID2': 'bulk2',
                            'bulk1': ['sampleID1', 'sampleID2'],
                            'bulk2': ...,
                            ...
                        }
    Returns:
        snpIndexDict -- a dictionary with structure like:
                        {
                            # TBD
                        }
    '''
    assert len(metadataDict["parent"]) == 1, \
        "vcf_to_snp_index will only accept one parent sample!"
    
    PARENT_ID = list(metadataDict["parent"])[0]
    
    snpIndexDict = {}
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle the header line
            if line.startswith("#CHROM"):
                sampleIDs = sl[9:]
                # Validate that our metadata matches
                assert all([ s in sampleIDs for v in metadataDict.values() for s in v ]), \
                    "The metadata file doesn't match this VCF!"
                continue
            
            # Skip other comment lines
            elif line.startswith("#"):
                continue
            
            # Handle content lines
            else:
                # Extract data out of the line
                chrom, pos, id, ref, alt, \
                    qual, filter, info, format = sl[0:9]
                samples = sl[9:]
                
                sampleFormat = format.split(":")
                
                # Parse relevant fields for each sample according to FORMAT
                samplesDict = {}
                for i in range(len(sampleIDs)):
                    sampleID = sampleIDs[i]
                    samplesDict.setdefault(sampleID, {})
                    sampleDetails = sl[9+i].split(":")
                    
                    ## > GT
                    gt = sampleDetails[sampleFormat.index("GT")]
                    samplesDict[sampleID]["GT"] = "/".join(["0"]*gt.count(".")) if set(gt.split("/")) == set(".") else gt
                    
                    ## > AD
                    ad = sampleDetails[sampleFormat.index("AD")]
                    samplesDict[sampleID]["AD"] = ",".join(["0"]*gt.count(".")) if ad == "." else ad
                    samplesDict[sampleID]["AD"] = samplesDict[sampleID]["AD"].replace(".", "0")
                
                # Get parent genotype and skip if it's a heterozygote / ambiguous
                """BSA largely relies on the premise that one of the parents is a
                homozygote. Since this script assumes a normal BSA analysis with
                QTL-seq for example, that means we only have one parent's genotype.
                Hence, it must be a homozygote for this analysis to have meaning!"""
                parentGT = samplesDict[PARENT_ID]["GT"]
                parentGTSet = set(parentGT.split("/"))
                if "." in parentGT or len(parentGTSet) != 1:
                    continue
                parentGT = int(list(parentGT)[0]) # get the value out of the set as an int
                
                # Compute population depth and SNP index values
                popDataDict = {}
                try:
                    for pop in ["bulk1", "bulk2"]: # skip the parent
                        # Calculate allele depth for each population
                        popAD = [0]*(len(alt)+1) # this is how many values will be in the AD field
                        for sampleID in metadataDict[pop]:
                            for index, adValue in enumerate(samplesDict[sampleID]["AD"].split(",")):
                                popAD[index] += int(adValue)
                        
                        # Calculate SNP index from the allele depth values
                        """QTL-seq makes it so that the delta SNP-index == -1 is where
                        causative SNPs lie. To make this so, you need the high-value bulk to
                        always have a SNP-index that evaluates as 0.
                        Code below based off of: 
https://github.com/YuSugihara/QTL-seq/blob/4b68ca1882fa0621d0996542250cac5633375437/qtlseq/snpfilt.py
                        """
                        if parentGT == 0:
                            altAD = sum(popAD[1:])
                        else:
                            altAD = popAD[0]
                        
                        popSnpIndex = altAD / sum(popAD) # sum of alternate / count of reads aligned
                        
                        # Store it
                        popDataDict[pop] = {
                            "depth": sum(popAD),
                            "snpIndex": popSnpIndex
                        }
                except ZeroDivisionError: # if this occurs, the SNP is garbage and must be skipped
                    continue
                
                # Filter and skip positions according to QTL-seq publication
                bothLowDepth = all([ v["depth"] < 7 for v in popDataDict.values() ])
                if bothLowDepth:
                    continue
                
                # Calculate delta SNP index
                deltaSnpIndex = popDataDict["bulk2"]["snpIndex"] - popDataDict["bulk1"]["snpIndex"]
                
                # Store result
                snpIndexDict.setdefault(chrom, [])
                snpIndexDict[chrom].append(list(map(str, [
                    pos, "snp" if len(ref) == len(alt) else "indel",
                    popDataDict["bulk1"]["depth"],
                    popDataDict["bulk2"]["depth"],
                    -99, -99,
                    popDataDict["bulk1"]["snpIndex"],
                    popDataDict["bulk2"]["snpIndex"],
                    deltaSnpIndex
                ])))
    return snpIndexDict

def main():
    # User input
    usage = """%(prog)s receives a VCF produced by Freebayes (or other software) which contains
    predictions per-sample and, using a metadata population map file, will bulk the predictions
    together. We are expecting there to be two bulks (bulk1, bulk2) and a parent sample (parent)
    listed in the population file.
    
    Note that the p95 and p99 values output by this script are nonsense. The calculation
    done by QTL-seq to produce these is harder for me to implement than I think is worth it.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input bulked VCF file")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Input pops map file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    args = p.parse_args()
    validate_args(args)
    
    # Parse metadata file
    metadataDict = parse_pops_metadata(args.metadataFile)
    assert set(metadataDict.keys()) == {'bulk1', 'bulk2', 'parent'}, \
        "Metadata file should only have 3 populations: bulk1, bulk2, and parent!"
    
    # Parse VCF file and generate SNP index dictionary structure
    snpIndexDict = vcf_to_snp_index(args.vcfFile, metadataDict)
    
    # Write to output file
    with open(args.outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "CHROM", "POSI", "variant", "bulk1_depth", "bulk2_depth",
            "p99", "p95", "bulk1_SNPindex", "bulk2_SNPindex", "delta_SNPindex"
        ])))
        # Write content lines
        for contig, positionLines in snpIndexDict.items():
            for pLine in positionLines:
                fileOut.write("{0}\t{1}\n".format(contig, "\t".join(
                    pLine
                )))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
