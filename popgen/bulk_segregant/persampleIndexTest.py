#! python3
# persampleIndexTest.py
# Script to test SNP proportions across two bulks
# and computes the SNP-index-like value to allow for
# assessments of allele frequencies and its potential
# correlation with the phenotype

import os, argparse, pickle, sys
import pandas as pd

# Load functions from other scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 1 dir back is where we find windows
import VCF

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.vcfFile):
        print(f'I am unable to locate the input VCF file ({args.vcfFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print(f'I am unable to locate the metadata file ({args.metadataFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric values
    if args.topPct <= 0:
        print("topPct must be a value greater than zero")
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print('The specified output file already exists. This program will not allow overwriting.')
        print('If you want to re-do an analysis, stop now and move/delete the existing file to fix this issue.')
        quit()

def parse_bulk_metadata(metadataFile):
    '''
    Parses a metadata file with two columns. The first is the sample ID,
    the second contains a string indicating whether it belongs to bulk1 or bulk2;
    in other words, bulk1 represents one extreme phenotype and bulk2 the other.
    
    Parameters:
        metadataFile -- a string pointing to the location of the metadata TSV
    Returns:
        metadataDict -- a dictionary with structure like:
                        {
                            'bulk1': set([
                                'sampleID1',
                                'sampleID2',
                                ...
                            ]),
                            'bulk2': set([
                                'sampleID10',
                                'sampleID11',
                                ...
                            ]),
                        }
    '''
    ALLOWED_VALUES = ["bulk1", "bulk2"]
    
    # Parse file
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if sl == []:
                continue
            else:
                sample, bulkGroup = sl
                
                if bulkGroup not in ALLOWED_VALUES:
                    print(f"'{bulkGroup}' is not recognised as a bulk grouping")
                    print("Fix your metadata file and try again!")
                    quit()
                else:
                    metadataDict.setdefault(bulkGroup, set())
                    if sample in metadataDict[bulkGroup]:
                        print(f"'{sample}' appears to be a duplicated ID!")
                        print("Fix your metadata file and try again.")
                        quit()

                    metadataDict[bulkGroup].add(sample)
    
    return metadataDict

def get_vcf_snpindex_density(vcfFile, metadataDict):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        metadataDict -- a dictionary with structure like:
                        {
                            'bulk1': set([
                                'sample1',
                                'sample2',
                                ...
                            ]),
                            'bulk2': set([
                                'sample10',
                                'sample11',
                                ...
                            ])
                        }
        windowSize -- an integer value indicating what size to bin genes in
    '''
    indexDict = {}
    with VCF.open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
            if line.startswith("#"):
                continue
            
            # Handle content lines
            if len(sl) >= 10:
                # Parse out relevant details from this line
                chrom, pos = sl[0:2]
                
                # Calculate the SNP-index-like value for this SNP
                snpIndex = VCF.snp_index_from_vcf_line(sl, samples, metadataDict)
                
                # Store value
                indexDict.setdefault(chrom, [])
                indexDict[chrom].append([pos, snpIndex])
      
    return indexDict

def main():
    usage = """%(prog)s receives a VCF and creates SNP density plots per
    chromosome. It specifically calculates the average SNP-index-like value
    over each window and plots that value. Hence, it may help to visualise
    where in the genome regions associated with bulks exist. As such, we
    require a metadata file as input
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Specify the location of the input VCF file")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Specify the metadata file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file where results will be written")
    # Opts
    p.add_argument("--top_pct", dest="topPct",
                   required=False,
                   type=float,
                   help="""Optionally, specify the top X percent of results
                   to output to file; those not meeting this threshold will be
                   excluded (default==100.0; all SNPs to be output)""",
                   default=100.0)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse metadata file
    metadataDict = parse_bulk_metadata(args.metadataFile)
    
    # Calculate SNP-index-like for each variant in the VCF
    pickleFile = args.vcfFile + "_index.pkl"
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            indexDict = pickle.load(pickleIn)
    else:
        indexDict = get_vcf_snpindex_density(args.vcfFile, metadataDict)
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(indexDict, pickleOut)
    
    # Tabulate SNP-index-like results and extract the top X percent if relevant
    dfDict = {"contig": [], "position": [], "index": []}
    for contigID, posList in indexDict.items():
        for position, corrValue in posList:
            dfDict["contig"].append(contigID)
            dfDict["position"].append(int(position))
            dfDict["index"].append(corrValue)
    df = pd.DataFrame(dfDict)
    
    if args.topPct != 100.0:
        quantile = 1 - (args.topPct / 100)
        df = df[df["index"].ge(df["index"].quantile(quantile))]
    
    # Create output TSV file
    df.to_csv(args.outputFileName, sep="\t", header=True, index=False)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
