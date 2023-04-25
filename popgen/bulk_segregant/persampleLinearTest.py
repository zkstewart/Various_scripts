#! python3
# persampleLinearTest.py
# Script to test SNP proportions across two bulks
# to assess whether a SNP has a putatively causal role
# in determining a quantitative phenotype

import os, argparse, pickle, gzip
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from contextlib import contextmanager

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

def parse_pheno_metadata(metadataFile):
    '''
    Parses a metadata file with two columns. The first is the sample ID,
    the second contains an integer or float value representing a phenotype.
    
    Parameters:
        metadataFile -- a string pointing to the location of the metadata TSV
    Returns:
        metadataDict -- a dictionary with structure like:
                        {
                            'sampleID1': phenotypeFloatValue,
                            'sampleID2': phenotypeFloatValue,
                            ...
                        }
    '''
    # Parse file
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if sl == []:
                continue
            else:
                sample, quantValue = sl
                
                try:
                    quantValue = float(quantValue)
                except:
                    print(f"{quantValue} in metadata file isn't an int or float value.")
                    print("Phenotype metadata is expected to conform to this assumption.")
                    print("Program must end now.")
                    quit()
                
                metadataDict[sample] = quantValue
    
    return metadataDict

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def get_phenotype_correlations(vcfFile, metadataDict):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        metadataDict -- a dictionary with structure like:
                        {
                            'sampleID1': phenotypeFloatValue,
                            'sampleID2': phenotypeFloatValue,
                            ...
                        }
    Returns:
        corrDict -- a dictionary with structure like:
                    {
                        'chrom1': [
                            [pos1, corrValue1],
                            [pos2, corrValue2],
                            ...
                        ],
                        'chrom2': [
                            ...
                        ],
                        ...
                    }
    '''
    corrDict = {}
    with open_vcf_file(vcfFile) as fileIn:
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
                chrom, pos, alt = *sl[0:2], sl[4]
                
                # Skip if line isn't biallelic
                if "," in alt:
                    continue
                
                # Calculate the SNP-index-like value for this SNP
                correlation = linear_model_from_vcf_line(sl, samples, metadataDict)
                
                # Store value
                corrDict.setdefault(chrom, [])
                corrDict[chrom].append([pos, correlation])
    
    return corrDict

def linear_model_from_vcf_line(sl, samples, metadataDict):
    '''
    Calculates the SNP-index-like value from a VCF line.
    
    Parameters:
        sl -- a split line from a VCF file that is guaranteed not to
              be a header/comment line
        samples -- a list containing strings corresponding to sample IDs
                   parsed out of the VCF #CHROM line
        metadataDict -- a dictionary with structure like:
                        {
                            'sampleID1': phenotypeFloatValue,
                            'sampleID2': phenotypeFloatValue,
                            ...
                        }
    '''
    ZYGOTE_ORDER = ["0/0", "0/1", "1/1"]
    
    # Parse out relevant details from this line
    format = sl[8].split(":")
    
    # Extract sample genotype values
    gtIndex = format.index("GT")
    gtDict = { samples[i]: sl[9+i].split(":")[gtIndex] for i in range(len(samples)) } # sl[9+i].split(":") gives the sample data in FORMAT layout
    
    # Get genotypes excluding blanks
    gts = set([ v for v in gtDict.values() if v != r"./." ])
    if len(gts) == 1: # if there's only 1 GT, we can't perform any correlation analysis
        return 0.0
    
    # Order genotypes homozygote 0/0 -> heterozygote 0/1 -> homozygote 1/1
    gts = sorted(list(gts), key = lambda x: ZYGOTE_ORDER.index(x))
    gtPhenoDict = { gt: [] for gt in gts } # set up storage area
    
    # Store phenotype values under their respective genotype
    for sampleID, phenoValue in metadataDict.items():
        sampleGT = gtDict[sampleID]
        if sampleGT in gtPhenoDict:
            gtPhenoDict[sampleGT].append(phenoValue)
    
    # Formulate values as X and Y axis numpy arrays
    x = np.array(
        [
            i
            for i in range(len(gtPhenoDict.keys()))
            for value in gtPhenoDict[list(gtPhenoDict.keys())[i]]
        ]
    ).reshape((-1, 1))
    
    y = np.array(
        [
            value
            for i in range(len(gtPhenoDict.keys()))
            for value in gtPhenoDict[list(gtPhenoDict.keys())[i]]
        ]
    )
    
    # Build linear model to test the R-squared fit
    model = LinearRegression().fit(x, y)
    r_sq = model.score(x, y)
    
    return r_sq

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
    metadataDict = parse_pheno_metadata(args.metadataFile)
    
    # Tally SNPs over windows per contig
    pickleFile = args.vcfFile + "_linear.pkl"
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            corrDict = pickle.load(pickleIn)
    else:
        corrDict = get_phenotype_correlations(args.vcfFile, metadataDict)
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(corrDict, pickleOut)
    
    # Tabulate correlation results and extract the top X percent if relevant
    dfDict = {"contig": [], "position": [], "r_squared": []}
    for contigID, posList in corrDict.items():
        for position, corrValue in posList:
            dfDict["contig"].append(contigID)
            dfDict["position"].append(int(position))
            dfDict["r_squared"].append(corrValue)
    df = pd.DataFrame(dfDict)
    
    if args.topPct != 100.0:
        quantile = 1 - (args.topPct / 100)
        df = df[df["r_squared"].ge(df["r_squared"].quantile(quantile))]
    
    # Create output TSV file
    df.to_csv(args.outputFileName, sep="\t", header=True, index=False)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
