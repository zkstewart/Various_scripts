#! python3
# convert_freebayes_vcf_to_bulked_qtlseqr_format.py
# Script to receive a Freebayes prediction VCF (performed per-sample)
# and convert it into files for the bulks (1 & 2) and for one or both
# parents. Produces re-formatted output amenable to loading in
# to the QTLseqr package.

import os, argparse
import pandas as pd
import numpy as np
from statistics import median, mean
from collections import Counter

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.pops):
        print('I am unable to locate the pops file (' + args.pops + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def parse_pop_map(popsFile):
    ACCEPTED_POPS = ["bulk1", "bulk2", "parent"]
    popDict = {}
    with open(popsFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l == "" or l.startswith("#"):
                continue
            
            sampleID, pop = l.split("\t") if "\t" in l else l.split(",")
            assert pop in ACCEPTED_POPS, \
                f"{pop} is not recognised as a valid population from {ACCEPTED_POPS}"
            
            # Index in dictionary
            popDict[sampleID] = pop
            popDict.setdefault(pop, [])
            popDict[pop].append(sampleID)
    return popDict

def bulk_vcf_to_file(vcfFile, popsDict, outputFileName):
    '''
    This function will parse a VCF file and return a dictionary which can be used
    for various purposes
    
    Parameters:
        vcfFile -- a string indicating the file location of the VCF file
        popsDict -- a dictionary with structure like:
                    {
                        'sampleID1': 'bulk1',
                        'sampleID2': 'bulk2',
                        'bulk1': ['sampleID1', 'sampleID2'],
                        'bulk2': ...,
                        ...
                    }
        outputFileName -- a string indicating the location to write the bulked
                          VCF to.
    '''
    PL_INTERPRETATION = {
        3: ["0/0", "0/1", "1/1"],
        4: ["0/0", "0/1", "1/1", "0/2", "1/2", "2/2"],
        5: ["0/0", "0/1", "1/1", "0/2", "1/2", "2/2",
            "0/3", "1/3", "2/3", "3/3"],
        6: ["0/0", "0/1", "1/1", "0/2", "1/2", "2/2",
            "0/3", "1/3", "2/3", "3/3", "0/4", "1/4",
            "2/4", "3/4", "4/4"]
    }
    GT_INTERPRETATION = { # this links to PL_INTERPRETATION
        2: 3,
        3: 4,
        4: 5,
        5: 6
    }
    
    with open(vcfFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        firstLine = True
        for line in fileIn:
            # Skip comment lines
            if line.startswith("#") and not line.startswith("#CHROM"):
                continue
            
            l = line.rstrip("\r\n ")
            sl = l.split("\t")
            
            # Handle the header line
            if firstLine == True:
                assert line.startswith("#CHROM") or line.startswith("CHROM"), \
                    "First non-comment line isn't the expected VCF header line format"
                
                sampleIDs = sl[9:]
                firstLine = False
            # Handle content lines
            else:
                chrom, pos, id, ref, alt, \
                    qual, filter, info, format = sl[0:9]
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
                    
                    ## > GL
                    gl = sampleDetails[sampleFormat.index("GL")]
                    samplesDict[sampleID]["GL"] = gl
                
                # Make into a Pandas dataframe
                samplesDF = pd.DataFrame.from_dict(samplesDict, orient="index")
                
                # # Fix GL when it's unrealistically big
                # "It's causing issues when using PL_INTERPRETATION"
                # uniqueAlleles = list(sorted(set([g for gt in set(samplesDF["GT"]) for g in gt.split("/")])))
                
                # # >> First, fix the genotype values if they're weird
                # foundFixes = False
                # for i in range(len(uniqueAlleles)):
                #     if int(uniqueAlleles[i]) != i:
                #         foundFixes = True
                #         for index, row in samplesDF.iterrows():
                #             splitRowGT = row["GT"].split("/")
                #             for x in range(len(splitRowGT)):
                #                 thisSplitGT = int(splitRowGT)[x]
                #                 if thisSplitGT == uniqueAlleles[i]:
                #                     thisSplitGT[x] = str(i)
                #             row["GT"] = "/".join(thisSplitGT)
                
                # # >> Then, fix the GL values
                # for index, row in samplesDF.iterrows():
                #     splitRowGL = row["GL"].split(",")
                
                # Impute GL fields
                longestGL = max([len(row.split(",")) for row in samplesDF["GL"]])
                
                # >> Get average GL for 0/0 genotypes
                averageGL = [[] for i in range(0, longestGL)]
                ## > Scenario where 0/0 genotypes exist
                for index, row in samplesDF.iterrows():
                    if row["GT"] == "0/0":
                        sGL = row["GL"].split(",")
                        for i in range(len(sGL)):
                            if sGL[i] != ".":
                                averageGL[i].append(float(sGL[i]))
                if [t for thing in averageGL for t in thing] != []:
                    averageGL = [sum(gl) / len(gl) if len(gl) > 0 else 0 for gl in averageGL]
                ## > Scenario where they do not exist
                else:
                    for index, row in samplesDF.iterrows():
                        sGL = row["GL"].split(",")
                        for i in range(len(sGL)):
                            if sGL[i] != ".":
                                averageGL[i].append(float(sGL[i]))
                    worstForImputation = min([mean(gl) for gl in averageGL if gl != []])
                    averageGL = [sum(gl) / len(gl) if len(gl) > 0 else worstForImputation for gl in averageGL]
                
                # >> Impute average GL for any missing genotypes
                samplesDF["GL"].replace(".", np.NaN, inplace=True)
                samplesDF["GL"].fillna(",".join(map(str, averageGL)), inplace=True)
                for index, row in samplesDF.iterrows():
                    if "." in row["GL"]:
                        splitRowGL = row["GL"].split(",")
                        for i in range(len(splitRowGL)):
                            if splitRowGL[i] == ".":
                                splitRowGL[i] = str(averageGL[i])
                        row["GL"] = ",".join(splitRowGL)
                                
                # Get bulked sample data frame
                bulk1_DF = samplesDF.loc[[sampleID for sampleID, group in popsDict.items() if group == "bulk1"],:]
                bulk2_DF = samplesDF.loc[[sampleID for sampleID, group in popsDict.items() if group == "bulk2"],:]
                parent_DF = samplesDF.loc[[sampleID for sampleID, group in popsDict.items() if group == "parent"],:]
                
                # Compute bulked values
                ## > AD
                alleles = list(sorted(set([allele for row in samplesDF["GT"] for allele in row.split("/")])))
                
                bulk1_allelesCounts = { allele: 0 for allele in alleles }
                for index, row in bulk1_DF.iterrows():
                    rowAlleles = row["GT"].split("/")
                    rowAd = row["AD"].split(",")
                    for i in range(len(rowAlleles)):
                        bulk1_allelesCounts[rowAlleles[i]] += int(rowAd[i])
                bulk1_AD = [bulk1_allelesCounts[a] for a in alleles]
                
                bulk2_allelesCounts = { allele: 0 for allele in alleles }
                for index, row in bulk2_DF.iterrows():
                    rowAlleles = row["GT"].split("/")
                    rowAd = row["AD"].split(",")
                    for i in range(len(rowAlleles)):
                        bulk2_allelesCounts[rowAlleles[i]] += int(rowAd[i])
                bulk2_AD = [bulk2_allelesCounts[a] for a in alleles]
                
                parent_allelesCounts = { allele: 0 for allele in alleles }
                for index, row in parent_DF.iterrows():
                    rowAlleles = row["GT"].split("/")
                    rowAd = row["AD"].split(",")
                    for i in range(len(rowAlleles)):
                        parent_allelesCounts[rowAlleles[i]] += int(rowAd[i])
                parent_AD = [parent_allelesCounts[a] for a in alleles]
                
                ## > GL
                bulk1_GL = [
                    bulk1_DF["GL"].map(lambda x: x.split(",")[i]).astype(float).mean()
                        for i in range(longestGL)
                ]
                bulk2_GL = [
                    bulk2_DF["GL"].map(lambda x: x.split(",")[i]).astype(float).mean()
                        for i in range(longestGL)
                ]
                parent_GL = [
                    parent_DF["GL"].map(lambda x: x.split(",")[i]).astype(float).mean()
                        for i in range(longestGL)
                ]
                
                ## > Normalise GL
                bulk1_GL = [_gl - max(bulk1_GL) for _gl in bulk1_GL]
                bulk2_GL = [_gl - max(bulk2_GL) for _gl in bulk2_GL]
                parent_GL = [_gl - max(parent_GL) for _gl in parent_GL]
                
                ## > PL (needs to be calculated from GL)
                "math from https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/12688-samtools-option-to-get-genotype-likelihood-gl#post160524"
                bulk1_PL = [abs(-_gl*10) for _gl in bulk1_GL]
                bulk2_PL = [abs(-_gl*10) for _gl in bulk2_GL]
                parent_PL = [abs(-_gl*10) for _gl in parent_GL]
                
                ## > GQ
                "GQ derived from PL according to https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format"
                bulk1_GQ = median(bulk1_PL)
                bulk2_GQ = median(bulk2_PL)
                parent_GQ = median(parent_PL)
                
                ## > DP
                bulk1_DP = sum(bulk1_AD)
                bulk2_DP = sum(bulk2_AD)
                parent_DP = sum(parent_AD)
                
                ## > GT
                bulk1_GT = PL_INTERPRETATION[len(bulk1_PL)][bulk1_PL.index(0.0)]
                bulk2_GT = PL_INTERPRETATION[len(bulk2_PL)][bulk2_PL.index(0.0)]
                parent_GT = PL_INTERPRETATION[len(parent_PL)][parent_PL.index(0.0)]
                
                # Write line
                outLine = [
                    *sl[0:8],
                    "GT:DP:AD:GL:PL:GQ",
                    "{parent_GT}:{parent_DP}:{parent_AD}:{parent_GL}:{parent_PL}:{parent_GQ}".format(
                        parent_GT=parent_GT, parent_DP=parent_DP, parent_AD=",".join(map(str, parent_AD)),
                        parent_GL=",".join(map(str, parent_GL)), parent_PL=",".join(map(str, parent_PL)),
                        parent_GQ=parent_GQ
                    ),
                    "{bulk1_GT}:{bulk1_DP}:{bulk1_AD}:{bulk1_GL}:{bulk1_PL}:{bulk1_GQ}".format(
                        bulk1_GT=bulk1_GT, bulk1_DP=bulk1_DP, bulk1_AD=",".join(map(str, bulk1_AD)),
                        bulk1_GL=",".join(map(str, bulk1_GL)), bulk1_PL=",".join(map(str, bulk1_PL)),
                        bulk1_GQ=bulk1_GQ
                    ),
                    "{bulk2_GT}:{bulk2_DP}:{bulk2_AD}:{bulk2_GL}:{bulk2_PL}:{bulk2_GQ}".format(
                        bulk2_GT=bulk2_GT, bulk2_DP=bulk2_DP, bulk2_AD=",".join(map(str, bulk2_AD)),
                        bulk2_GL=",".join(map(str, bulk2_GL)), bulk2_PL=",".join(map(str, bulk2_PL)),
                        bulk2_GQ=bulk2_GQ
                    )
                ]
                fileOut.write("\t".join(outLine) + "\n")

def main():
    # User input
    usage = """%(prog)s receives a VCF produced by Freebayes which contains predictions
    per-sample and, using a metadata population map file, will bulk the predictions together.
    We are expecting there to be two bulks (bulk1, bulk2) and a parent sample (parent) listed
    in the population file.
    
    Note: Samples not found in the metadata file will be ignored. Use this as a feature if you
    want to exclude some samples.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input bulked VCF file")
    p.add_argument("-p", dest="pops",
                   required=True,
                   help="Input pops map file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    args = p.parse_args()
    validate_args(args)
    
    # Parse population map file
    popsDict = parse_pop_map(args.pops)
    
    # Parse VCF file and bulk each line directly to output file
    bulk_vcf_to_file(args.vcf, popsDict, args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
