#! python3
# snpDensityPlot.py
# Script to create visualisations of SNP density
# for assessing hypotheses of SNP distribution
# across chromosomes.

import os, argparse, math, sys, re
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.ndimage.filters import gaussian_filter1d
from contextlib import contextmanager

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_GFF3IO, ZS_VCFIO

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
    
    # Handle numeric parameters
    if args.width < 0:
        raise ValueError("width must be a positive integer")
    if args.height < 0:
        raise ValueError("height must be a positive integer")
    
    # Handle regions
    for region in args.regions:
        if not re.match(r"^.+:\d+:\d+$", region):
            raise ValueError(f"Region '{region}' is improperly formatted" +
                             "Please provide regions in the format 'contig:start:end' and try again.")
    
    # Check for logical consistency of optional argument combination
    if args.onePlot and args.regions != []:
        raise ValueError("You can't provide both --onePlot and --regions; please choose one and try again.")
    if args.onePlot and args.mode == "gene":
        raise ValueError("--onePlot is not supported for gene mode.")
    if args.filter != []:
        if not args.weightByOccurrence and not args.skipMonoallelic:
            raise ValueError("You have provided a list of samples to filter on, " + 
                             "but you have not enabled either the --weightByOccurrence " + 
                             "or --skipMonoallelic flags. Please enable at least one of these " + 
                             "flags to proceed, or omit --filter values.")
    
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
    elif os.path.exists(args.outputDirectory):
        raise ValueError(f"The specified output directory ({args.outputDirectory}) is not a directory!")
    else:
        os.makedirs(args.outputDirectory, exist_ok=True)
        print(f"Output directory created at '{args.outputDirectory}' as part of argument validation.")
    
    # Handle mode-specific arguments
    if args.mode == "line":
        if args.wmaSize < 1:
            raise ValueError("wmaSize must be a positive integer")
        if args.lineWidth < 1:
            raise ValueError("lineWidth must be a positive integer")
        if args.windowSize < 1:
            raise ValueError("windowSize must be a positive integer")
        if args.minimumContigSize < (args.windowSize * 2):
            raise ValueError(f"minimumContig must be an integer >= 2* the window size ({(args.windowSize * 2)})")
    elif args.mode == "histogram":
        if args.binSize < 1:
            raise ValueError("binSize must be an integer >= 1")
    elif args.mode == "gene":
        if not os.path.isfile(args.gff3File):
            raise ValueError(f"I am unable to locate the input GFF3 file ({args.gff3File})")
    else:
        raise ValueError("Invalid mode specified")

def tally_variants_within_features(vcfFile, gff3, onlyCDS=False, weightByOccurrence=False,
                                   skipMonoallelic=False, filterSamples=[]):
    '''
    To future self or others: the interaction of skipMonoallelic and filterSamples can mean that you
    don't get an exact deduction of the number of variants that occurs in a sample. This is because
    if you e.g., filter all but one sample, you'd only count their heterozygous variants because
    of skipMonoallelic. This can be tricky to account for when testing.
    
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        gff3 -- a ZS_GFF3IO.GFF3 object
        onlyCDS -- OPTIONAL; a boolean flag to indicate whether to only count
                   variants within CDS regions (default == False)
        weightByOccurrence -- OPTIONAL; a boolean flag to indicate whether to
                              count a variant each time it occurs in a sample
                              (default == False)
        skipMonoallelic -- OPTIONAL; a boolean flag to indicate whether to skip
                           counting a variant when ALL samples have a non-variant
                           allele (default == False)
        filterSamples -- OPTIONAL; a list of strings indicating the sample IDs
                         to filter on. Only relevant if you are using the
                         --weightByOccurrence flag or the --skipMonoallelic flag
                         (default == [])
    Returns:
        tallyDict -- a dictionary with structure like:
                     {
                         "geneID1": [snpCount, geneLength],
                         "geneID2": [snpCount, geneLength],
                         ...
                     }
    '''
    # Establish tally dictionary
    tallyDict = {}
    for geneFeature in gff3.types["gene"]:
        geneID = geneFeature.ID
        if onlyCDS:
            longestMrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(geneFeature)
            geneLength = sum([ cdsFeature.end - cdsFeature.start + 1 for cdsFeature in longestMrnaFeature.CDS ])
        else:
            geneLength = geneFeature.end - geneFeature.start + 1
        tallyDict[geneID] = [0, geneLength]
    
    # Populate tally dictionary
    with ZS_VCFIO.open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # replace quotations may help with Excel files
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
                for sample in filterSamples:
                    if not sample in samples:
                        raise ValueError(f"Sample indicated for filtration ({sample}) not found in VCF header!")
                filterIndices = set([ samples.index(sample) for sample in filterSamples ])
                continue
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            chrom = sl[0]
            pos = int(sl[1])
            sampleDetails = sl[9:]
            
            fieldsDescription = sl[8]
            if ":" not in fieldsDescription:
                gtIndex = 0
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Iterate through any gene features this variant overlaps
            geneFeatures = gff3.ncls_finder(pos, pos, "contig", chrom)
            for geneFeature in geneFeatures:
                geneID = geneFeature.ID
                
                # Narrow down to variants within the CDS region if onlyCDS
                if onlyCDS:
                    mrnaFeatures = [
                        childFeature
                        for childFeature in geneFeature.children
                        if childFeature.type == "mRNA"
                    ]
                    
                    cdsFeatures = [
                        childFeature
                        for mrnaFeature in mrnaFeatures
                        for childFeature in mrnaFeature.children
                        if childFeature.type == "CDS"
                    ]
                    
                    isWithinCDS = any([ cdsFeature.start <= pos <= cdsFeature.end for cdsFeature in cdsFeatures ])
                    if not isWithinCDS:
                        continue
                
                # Get sample genotypes if relevant
                if skipMonoallelic or weightByOccurrence:
                    genotypes = [
                        sampleDetails[i].split(":")[gtIndex].replace("|", "/").split("/")
                        for i in range(len(sampleDetails))
                        if i not in filterIndices
                    ]
                
                # Skip if we're skipping monoallelic and all samples are the same
                if skipMonoallelic:
                    genotypeSet = set([ g for genotype in genotypes for g in genotype ]).difference(".") # remove missing data
                    if len(genotypeSet) == 1:
                        continue
                
                # See how many times we should count this variant if relevant
                if weightByOccurrence:
                    snpCount = sum([ 1 for genotype in genotypes if set(genotype).difference(".") != {"0"} ])
                else:
                    snpCount = 1
                
                # Add to the tally
                tallyDict[geneID][0] += snpCount

def get_vcf_density(vcfFile, lengthsDict, windowSize=100000, weightByOccurrence=False,
                    skipMonoallelic=False, filterSamples=[]):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        windowSize -- OPTIONAL; an integer value indicating what size to bin variants
                      within (default == 100000)
        weightByOccurrence -- OPTIONAL; a boolean flag to indicate whether to
                              count a variant each time it occurs in a sample
                              (default == False) or just once per variant
        skipMonoallelic -- OPTIONAL; a boolean flag to indicate whether to skip
                           counting a variant when ALL samples have a non-variant
                           allele (default == False); applies AFTER filterSamples
        filterSamples -- OPTIONAL; a list of strings indicating the sample IDs
                         to filter on. Only relevant if you are using the
                         --weightByOccurrence flag or the --skipMonoallelic flag
                         (default == [])
    '''
    densityDict = {}
    with ZS_VCFIO.open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # replace quotations may help with Excel files
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
                for sample in filterSamples:
                    if not sample in samples:
                        raise ValueError(f"Sample indicated for filtration ({sample}) not found in VCF header!")
                filterIndices = set([ samples.index(sample) for sample in filterSamples ])
                continue
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            chrom = sl[0]
            pos = int(sl[1])
            sampleDetails = sl[9:]
            
            fieldsDescription = sl[8]
            if ":" not in fieldsDescription:
                gtIndex = 0
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Get sample genotypes if relevant
            if skipMonoallelic or weightByOccurrence:
                genotypes = [
                    sampleDetails[i].split(":")[gtIndex].replace("|", "/").split("/")
                    for i in range(len(sampleDetails))
                    if i not in filterIndices
                ]
            
            # Skip if we're skipping monoallelic and all samples are the same
            if skipMonoallelic:
                genotypeSet = set([ g for genotype in genotypes for g in genotype ]).difference(".") # remove missing data
                if len(genotypeSet) == 1:
                    continue
            
            # See how many times we should count this variant if relevant
            if weightByOccurrence:
                snpCount = sum([ 1 for genotype in genotypes if set(genotype).difference(".") != {"0"} ])
            else:
                snpCount = 1
            
            # Store the variant in the relevant window tally
            if len(sl) >= 10:
                contigID, position = sl[0:2]
                densityDict.setdefault(contigID,
                    [ 0
                    for windowChunk in range(math.ceil(lengthsDict[contigID] / windowSize))
                    ]
                )
                windowChunkIndex = math.floor(int(position) / windowSize)
                densityDict[contigID][windowChunkIndex] += snpCount
    return densityDict

def normalise_density(densityList):
    '''
    Using min-max normalisation
    
    Parameters:
        densityList -- a list containing integers of any length
    '''
    minValue, maxValue = min(densityList), max(densityList)
    normalisedDensityList = [ (x - minValue) / (maxValue - minValue) for x in densityList ]
    return normalisedDensityList

def ideo_per_contig(densityDict, wmaSize, width, height, outputDirectory, plotPDF, showDots, linewidth):
    '''
    Parameters:
        densityDict -- a dictionary with structure like:
                       {
                           'contig1': [intDensity1, intDensity2, ...],
                           'contig2': [intDensity1, intDensity2, ...],
                            ...
                        }
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
        showDots -- a boolean flag indicating whether to show dots of each data point
                    in addition to the line plot
        linewidth -- an integer value indicating the line width
    '''
     

def main():
    usage = """%(prog)s receives a VCF and creates SNP density plots per
    chromosome.
    
    Note: if you want to filter the contigs to plot, you can filter the genome FASTA
    to only contain the contigs you want to plot.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Specify the location of the input VCF file")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where plot files will be written")
    ## Optional (styling)
    p.add_argument("--width", dest="width",
                    type=int,
                    required=False,
                    help="""Optionally, specify the output plot width (default=10)""",
                    default=10)
    p.add_argument("--height", dest="height",
                    type=int,
                    required=False,
                    help="""Optionally, specify the output plot height (default=6)""",
                    default=6)
    p.add_argument("--onePlot", dest="onePlot",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want a single plot to be
                    produced with all chromosomes""",
                    default=False)
    p.add_argument("--pdf", dest="plotPDF",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want outputs to be
                    in PDF format instead of PNG format""",
                    default=False)
    ## Optional (behaviour)
    p.add_argument("--weightByOccurrence", dest="weightByOccurrence",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want to count a variant
                    each time it occurs in a sample""",
                    default=False)
    p.add_argument("--skipMonoallelic", dest="skipMonoallelic",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want to skip counting a
                    variant when ALL samples have a non-variant allele""",
                    default=False)
    p.add_argument("--filter", dest="filter",
                    required=False,
                    nargs="+",
                    help="""Optionally, indicate one or more samples to filter on. Only
                    relevant if you are using the --weightByOccurrence flag or the
                    --skipMonoallelic flag""",
                    default=[])
    p.add_argument("--regions", dest="regions",
                    required=False,
                    nargs="+",
                    help="""Optionally, indicate one or more regions to plot in greater detail
                    by providing the contig ID and start and end positions in bp (e.g.
                    contig1:10000:20000)""",
                    default=[])
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser()
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    lineparser = subparsers.add_parser("line",
                                       parents=[p],
                                       add_help=False,
                                       help="Create line plots of SNP density")
    lineparser.set_defaults(func=linemain)
    
    histoparser = subparsers.add_parser("histogram",
                                        aliases=["histo"],
                                        parents=[p],
                                        add_help=False,
                                        help="Create histograms of SNP density")
    histoparser.set_defaults(func=histomain)
    
    ideoparser = subparsers.add_parser("ideogram",
                                       aliases=["ideo"],
                                       parents=[p],
                                       add_help=False,
                                       help="Create ideogram of SNP density per window")
    ideoparser.set_defaults(func=ideomain)
    
    geneparser = subparsers.add_parser("gene",
                                        parents=[p],
                                        add_help=False,
                                        help="Create broken bar plots of SNP density per gene")
    geneparser.set_defaults(func=genemain)
    
    # Line-subparser arguments
    lineparser.add_argument("--wmaSize", dest="wmaSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the number of previous values to consider
                            during weighted moving average calculation (default=5)""",
                            default=5)
    lineparser.add_argument("--windowSize", dest="windowSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the size of the window to sum
                            SNPs within (default=100000)""",
                            default=100000)
    lineparser.add_argument("--lineWidth", dest="lineWidth",
                            type=int,
                            required=False,
                            help="""Optionally, specify the line width (default=1)""",
                            default=1)
    lineparser.add_argument("--showDots", dest="showDots",
                            required=False,
                            action="store_true",
                            help="""Optionally, provide this flag if you want plots to show dots
                            of each data point in addition to the line plot""",
                            default=False)
    
    # Histogram-subparser arguments
    histoparser.add_argument("--binSize", dest="binSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the bin size to count variants
                            within (default=10000)""",
                            default=10000)
    
    # Gene-subparser arguments
    geneparser.add_argument("-g", dest="gff3File",
                            required=True,
                            help="Specify the location of the input GFF3 file")
    geneparser.add_argument("--onlyCDS", dest="onlyCDS",
                            required=False,
                            action="store_true",
                            help="""Optionally, provide this flag if you want to only count
                            variants within CDS regions""",
                            default=False)
    
    args = subParentParser.parse_args()
    validate_args(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Split into mode-specific functions
    if args.mode == "line":
        linemain(args, lengthsDict)
    elif args.mode in ["histogram", "histo"]:
        histomain(args, lengthsDict)
    elif args.mode == ["ideogram", "ideo"]:
        ideomain(args, lengthsDict)
    elif args.mode == "gene":
        genemain(args, lengthsDict)
    
    print("Program completed successfully!")

def linemain(args, lengthsDict):
    # Calculate SNP density
    densityDict = get_vcf_density(args.vcfFile, lengthsDict, args.windowSize,
                                  args.weightByOccurrence, args.skipMonoallelic,
                                  args.filter)
    
    # Create plots
    if args.onePlot:
        line_horizontal(densityDict, args.wmaSize,
                            args.width, args.height,
                            args.outputDirectory, args.plotPDF,
                            args.showDots, args.linewidth)
    elif args.regions != []:
        line_regions(densityDict, args.regions, args.wmaSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF,
                         args.showDots, args.linewidth)
    else:
        line_per_contig(densityDict, args.wmaSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.showDots, args.linewidth)

def histomain(args, lengthsDict):
    # Calculate SNP density
    densityDict = get_vcf_density(args.vcfFile, lengthsDict, args.windowSize,
                                  args.weightByOccurrence, args.skipMonoallelic,
                                  args.filter)
    
    # Create plots
    if args.onePlot:
        histo_horizontal(densityDict, args.binSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF,
                         args.showDots, args.linewidth)
    elif args.regions != []:
        histo_regions(densityDict, args.regions, args.binSize,
                      args.width, args.height,
                      args.outputDirectory, args.plotPDF,
                      args.showDots, args.linewidth)
    else:
        histo_per_contig(densityDict, args.binSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF,
                         args.showDots, args.linewidth)

def ideomain(args, lengthsDict):
    # Calculate SNP density
    densityDict = get_vcf_density(args.vcfFile, lengthsDict, args.windowSize,
                                  args.weightByOccurrence, args.skipMonoallelic,
                                  args.filter)
    
    # Create plots
    if args.onePlot:
        ideo_horizontal(densityDict, args.wmaSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.showDots, args.linewidth)
    elif args.regions != []:
        ideo_regions(densityDict, args.regions, args.wmaSize,
                     args.width, args.height,
                     args.outputDirectory, args.plotPDF,
                     args.showDots, args.linewidth)
    else:
        ideo_per_contig(densityDict, args.wmaSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.showDots, args.linewidth)

def genemain(args, lengthsDict):
    # Parse GFF3 with NCLS indexing
    gff3 = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    gff3.create_ncls_index(typeToIndex="gene")
    
    # Calculate SNP density
    tallyDict = tally_variants_within_features(args.vcfFile, gff3,
                                               args.onlyCDS, args.weightByOccurrence,
                                               args.skipMonoallelic, args.filter)
    
    # Create plots
    if args.onePlot:
        raise NotImplementedError("One plot not supported for gene mode")
    elif args.regions != []:
        gene_regions(tallyDict, args.regions, args.wmaSize,
                     args.width, args.height,
                     args.outputDirectory, args.plotPDF,
                     args.showDots, args.linewidth)
    else:
        gene_per_contig(tallyDict, args.wmaSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.showDots, args.linewidth)

def junkyard():
    # Tally SNPs over windows per contig
    densityDict = get_vcf_density(args.vcfFile, lengthsDict, args.windowSize)
    
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
            normalisedDensityList = normalise_density(densityList)
            
            # Smooth the curve for better visualisation
            smoothedDensityList = gaussian_filter1d(normalisedDensityList, sigma=args.smoothingSigma)
            
            # Configure plot
            kbpWindowSize = round(args.windowSize / 1000, 2)
            fig = plt.figure(figsize=(10,6))
            ax = plt.axes()
            
            ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows)", fontweight="bold")
            ax.set_ylabel("Min-max normalised SNP number per window", fontweight="bold")
            ax.set_title(f"{contigID} SNP density plot", fontweight="bold")
            
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
