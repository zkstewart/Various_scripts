#! python3
# snpIndexDensityPlot.py
# Script to create visualisations of SNP-index-like
# statistic density for assessing hypotheses of SNP
# selection across chromosomes.

import os, argparse, math, gzip, sys, pickle
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.ndimage.filters import gaussian_filter1d
from contextlib import contextmanager

# Load functions from other scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find windows
import bulk_segregant

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

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def get_vcf_snpindex_density(vcfFile, lengthsDict, metadataDict, windowSize=100000):
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
                chrom, pos = sl[0:2]
                
                # Calculate the SNP-index-like value for this SNP
                snpIndex = snp_index_from_vcf_line(sl, samples, metadataDict)
                windowChunkIndex = math.floor(int(pos) / windowSize)
                
                # Add SNP-index
                indexDict.setdefault(chrom, {
                    "indices": [ 0
                        for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                    ],
                    "counts": [ 0
                        for windowChunk in range(math.ceil(lengthsDict[chrom] / windowSize))
                    ],
                    }
                )
                indexDict[chrom]["indices"][windowChunkIndex] += snpIndex
                indexDict[chrom]["counts"][windowChunkIndex] += 1
    
    # Average the SNP-index per window
    densityDict = {}
    for chrom in indexDict.keys():
        densityDict.setdefault(chrom, [])
        for windowChunkIndex in range(len(indexDict[chrom]["indices"])):
            try:
                average = indexDict[chrom]["indices"][windowChunkIndex] / indexDict[chrom]["counts"][windowChunkIndex]
            except:
                average = 0.0
            densityDict[chrom].append(average)
      
    return densityDict

def snp_index_from_vcf_line(sl, samples, metadataDict):
    '''
    Calculates the SNP-index-like value from a VCF line.
    
    Parameters:
        sl -- a split line from a VCF file that is guaranteed not to
              be a header/comment line
        samples -- a list containing strings corresponding to sample IDs
                   parsed out of the VCF #CHROM line
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
    '''
    # Parse out relevant details from this line
    format = sl[8].split(":")
    
    # Extract sample genotype values
    gtIndex = format.index("GT")
    gtDict = { samples[i]: sl[9+i].split(":")[gtIndex] for i in range(len(samples)) } # sl[9+i].split(":") gives the sample data in FORMAT layout
    
    # Get genotypes excluding blanks
    gts = set([ v for v in gtDict.values() if v != r"./." ])
    if len(gts) == 1: # if there's only 1 GT, we always have a difference ratio of 0
        return 0.0
    
    # Set up data storage for allele counts
    alleles = set([ allele for gt in gts for allele in gt.split("/") ])
    
    alleleDict = {}
    for pop in metadataDict.keys():
        alleleDict[pop] = { allele:0 for allele in alleles }
    
    # Tally alleles for each population
    anyFound = False
    for key, value in gtDict.items():
        # Skip if GT is blank
        if value == r"./.":
            continue
        
        # Figure out what population this sample belongs to
        pop = [ k for k,v in metadataDict.items() if key in v ]
        if len(pop) != 1:
            continue
        else:
            anyFound = True
        pop = pop[0]
        
        # Figure out and store this genotype
        for allele in value.split("/"):
            alleleDict[pop][allele] += 1
    
    if anyFound == False:
        print("ERROR: Your metadata file doesn't match your VCF")
        print("In other words, we failed to find ANY samples from your metadata file in the VCF")
        print("This is an irreconcilable error and we must exit out now...")
        quit()
    
    # Calculate the proportions of each allele for this locus
    bulk1Sum = sum(alleleDict["bulk1"].values())
    bulk2Sum = sum(alleleDict["bulk2"].values())
    
    if bulk1Sum == 0 or bulk2Sum == 0: # if this happens, we can't meaningfully calculate anything
        return 0.0
    
    alleleProportions = {
        "bulk1": { allele: (alleleCount / bulk1Sum) for allele, alleleCount in alleleDict["bulk1"].items() },
        "bulk2": { allele: (alleleCount / bulk2Sum) for allele, alleleCount in alleleDict["bulk2"].items() }
    }
    
    # Derive our difference ratio value
    proportionCommon = sum([
        min(alleleProportions["bulk1"][str(allele)], alleleProportions["bulk2"][str(allele)])
        for allele in range(len(alleleProportions["bulk1"].keys()))
    ])
    differenceRatio = 1 - proportionCommon
    
    return differenceRatio

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
    metadataDict = bulk_segregant.parse_pops_metadata(args.metadataFile)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Tally SNPs over windows per contig
    pickleFile = args.vcfFile + "_indexDensity.pkl"
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as pickleIn:
            densityDict = pickle.load(pickleIn)
    else:
        densityDict = get_vcf_snpindex_density(args.vcfFile, lengthsDict, metadataDict, args.windowSize)
        with open(pickleFile, "wb") as pickleOut:
            pickle.dump(densityDict, pickleOut)
    
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
            ax.set_ylabel("Mean SNP-index-like value per window", fontweight="bold")
            ax.set_title(f"{contigID} SNP-index plot", fontweight="bold")
            
            ax.plot(smoothedDensityList)
            
            # Save output file
            plt.savefig(fileOut)
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
