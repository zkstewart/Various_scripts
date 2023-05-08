#! python3
# persampleLinearWindows.py
# Script to test SNP proportions across two bulks
# to assess whether a SNP has a putatively causal role
# in determining a quantitative phenotype. Takes the extra
# step of trying to plot out windows where these SNPs
# congregate in the genome.

import os, argparse, pickle, gzip, math
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from scipy.ndimage.filters import gaussian_filter1d
from sklearn.linear_model import LinearRegression
from contextlib import contextmanager

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.vcfFile):
        print(f'I am unable to locate the input VCF file ({args.vcfFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        print(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print(f'I am unable to locate the metadata file ({args.metadataFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.windowSize < 1:
        print("windowSize must be a positive integer")
        quit()
    if args.minimumContigSize < (args.windowSize * 2):
        print("minimumContigSize must be at least 2x the window size")
        quit()
    if args.smoothingSigma < 0:
        print("smoothingSigma must be a float greater than or equal to zero")
        quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
        print('This program will not allowing overwriting. If you want to re-do an analysis, stop now and fix this issue.')

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
    Returns:
        corrDict -- a dictionary with structure like:
                    {
                        'contigID1': [
                            ['pos1', r_sq_correlation],
                            ['pos2', r_sq_correlation],
                            ...
                        ],
                        'contigID2': [
                            ...
                        ],
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

def get_density_from_pheno_corr(corrDict, lengthsDict, windowSize=100000):
    '''
    Parameters:
        corrDict -- a dictionary with structure like:
                    {
                        'contigID1': [
                            ['pos1', r_sq_correlation],
                            ['pos2', r_sq_correlation],
                            ...
                        ],
                        'contigID2': [
                            ...
                        ],
                        ...
                    }
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        windowSize -- an integer value indicating what size to bin genes in
    '''
    summingDict = {}
    for chrom, corrPairs in corrDict.items():
        for pos, corr in corrPairs:
            windowChunkIndex = math.floor(int(pos) / windowSize)

            # Add correlation value
            summingDict.setdefault(chrom, {
                "corrs": [ 0
                    for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                ],
                "counts": [ 0
                    for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                ],
                }
            )
            summingDict[chrom]["corrs"][windowChunkIndex] += corr
            summingDict[chrom]["counts"][windowChunkIndex] += 1
    
    # Average the correlation value per window
    densityDict = {}
    for chrom in summingDict.keys():
        densityDict.setdefault(chrom, [])
        for windowChunkIndex in range(len(summingDict[chrom]["corrs"])):
            try:
                average = summingDict[chrom]["corrs"][windowChunkIndex] / summingDict[chrom]["counts"][windowChunkIndex]
            except:
                average = 0.0
            densityDict[chrom].append(average)
    
    return densityDict

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
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Specify the metadata file name")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where plot files will be written")
    # Opts
    p.add_argument("--window_size", dest="windowSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the size of the window to sum
                   SNPs within (default=50000)""",
                   default=50000)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the minimum size of contigs which
                   should have plots created for (default=200000)"; this value
                   must exceed --window_size by at least 2x""",
                   default=200000)
    p.add_argument("--smoothing", dest="smoothingSigma",
                   type=float,
                   required=False,
                   help="""TURNED OFF: Optionally, specify how much smoothing should be applied
                   to the plot; should be float value >= 0.0 (default=1.0)""",
                   default=1.0)
    
    args = p.parse_args()
    validate_args(args)
    os.makedirs(args.outputDirectory, exist_ok=True)
    
    # Parse metadata file
    metadataDict = parse_pheno_metadata(args.metadataFile)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Tally SNPs over windows per contig
    pickleFile = args.vcfFile + "_linear.pkl"
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            corrDict = pickle.load(pickleIn)
    else:
        corrDict = get_phenotype_correlations(args.vcfFile, metadataDict)
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(corrDict, pickleOut)
    
    # Get the average of correlation values in windows
    densityDict = get_density_from_pheno_corr(corrDict, lengthsDict, args.windowSize)
    
    # Create plot per contig
    numContigsProcessed = 0
    numContigsPlotted = 0
    for contigID, length in lengthsDict.items():
        if length >= args.minimumContigSize:
            numContigsProcessed += 1
            
            # Derive our output file name and skip if already existing
            fileOut = os.path.join(args.outputDirectory, f"{contigID}.png")
            if os.path.isfile(fileOut):
                print(f"WARNING: Plot for '{contigID}' already found in output directory; skipping...")
                continue
            
            # Skip if we found no SNPs on this contig
            if not contigID in densityDict:
                print(f"WARNING: '{contigID}' is in the VCF header but has no SNPs associated " +
                      "with it; skipping...")
                continue
            
            # Get density values
            densityList = densityDict[contigID]
            
            # Smooth the curve for better visualisation
            "Smoothing is a mistake for this kind of analysis"
            # smoothedDensityList = gaussian_filter1d(densityList, sigma=args.smoothingSigma)
            smoothedDensityList = densityList
            
            # Configure plot
            kbpWindowSize = round(args.windowSize / 1000, 2)
            fig = plt.figure(figsize=(10,6))
            ax = plt.axes()
            
            ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows)", fontweight="bold")
            ax.set_ylabel("Mean R^2 value per window", fontweight="bold")
            ax.set_title(f"{contigID} SNP-index plot", fontweight="bold")
            
            ax.plot(smoothedDensityList)
            
            # Save output file
            plt.savefig(fileOut)
            plt.close()
            numContigsPlotted += 1
    
    # Raise relevant warnings
    if numContigsProcessed == 0:
        print(f"WARNING: We didn't find any contigs which exceeded {args.minimumContigSize}bp in size")
        print("Hence, no output files have been generated! Maybe you should fix your --minimum_contig value?")
    elif numContigsPlotted == 0:
        print("WARNING: We ended up skipping every contig! This means the program has already run to completion previously.")
        print("Hence, no new output files have been generated! Maybe you should delete the existing folder to restart?")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
