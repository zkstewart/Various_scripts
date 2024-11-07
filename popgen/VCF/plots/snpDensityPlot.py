#! python3
# snpDensityPlot.py
# Script to create visualisations of SNP density
# for assessing hypotheses of SNP distribution
# across chromosomes.

import os, argparse, math, sys, re, pickle
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

from Bio import SeqIO
from hashlib import md5

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
    if not args.mode in ["genegram", "gene"]:
        if args.onePlot and args.regions != []:
            raise ValueError("You can't provide both --onePlot and --regions; please choose one and try again.")
    
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
    elif args.mode in ["histogram", "histo"]:
        if args.binSize < 1:
            raise ValueError("binSize must be an integer >= 1")
    elif args.mode in ["ideogram", "ideo"]:
        if args.binSize < 1:
            raise ValueError("binSize must be an integer >= 1")
    elif args.mode in ["genegram", "gene"]:
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
    return tallyDict

def bin_vcf_variants(vcfFile, lengthsDict, windowSize=100000, weightByOccurrence=False,
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
    Returns:
        binDict -- a dictionary with structure like:
                   {
                       'contig1': [binInt1, binInt2, ...],
                       'contig1': [...],
                       ...
                   }
    '''
    binDict = {}
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
                binDict.setdefault(contigID,
                    [ 0
                    for windowChunk in range(math.ceil(lengthsDict[contigID] / windowSize))
                    ]
                )
                windowChunkIndex = math.floor(int(position) / windowSize)
                binDict[contigID][windowChunkIndex] += snpCount
    return binDict

def get_vcf_variants(vcfFile, weightByOccurrence=False,
                    skipMonoallelic=False, filterSamples=[]):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
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
    Returns:
        dotsX -- a dict pairing chromosome IDs (keys) to a list of integers indicating
                 the position for each SNP (values)
        dotsY -- a dict pairing chromosome IDs (keys) to a list of floats indicating the
                 count for each SNP (values)
    '''
    dotsX, dotsY = [], []
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
            contig = sl[0]
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
            
            # Store data
            if not contig in dotsX:
                dotsX[contig] = []
                dotsY[contig] = []
            dotsX[contig].append(pos)
            dotsY[contig].append(snpCount)
    return dotsX, dotsY

def normalise_density(densityList):
    '''
    Using min-max normalisation
    
    Parameters:
        densityList -- a list containing integers of any length
    '''
    minValue, maxValue = min(densityList), max(densityList)
    normalisedDensityList = [ (x - minValue) / (maxValue - minValue) for x in densityList ]
    return normalisedDensityList

def histo_horizontal(binDict, binSize, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        binDict -- a dictionary with structure like:
                   {
                       'contig1': [intDensity1, intDensity2, ...],
                       'contig2': [intDensity1, intDensity2, ...],
                       ...
                   }
        binSize -- an integer value indicating the bin size to count variants
                   within
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''
    # Derive our output file name and skip if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    fileOut = os.path.join(outputDirectory, f"onePlot.histogram.{fileSuffix}")
    if os.path.isfile(fileOut):
        raise FileExistsError(f"'onePlot.histogram.{fileSuffix}' already found in output directory")
    
    # Get each contigs' plot data
    plotData = []
    for contigID in binDict.keys():
        y = binDict[contigID]
        x = np.arange(0, len(y))
        
        plotData.append([contigID, x, y])
    
    # Configure plot
    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    gs = fig.add_gridspec(1, len(plotData), hspace=0)
    axes = gs.subplots(sharey='row')
    
    ## Set the figure title
    stepSize = binSize / 1000 # convert to Kbp
    fig.suptitle(f"SNP occurrence histogram", fontweight="bold")
    fig.supxlabel(f"Bin number ({stepSize} Kbp step and width)", fontweight="bold")
    fig.supylabel(f"Number of variants", fontweight="bold")
    
    # Plot the data into each axis
    for ax, (contigID, x, y) in zip(axes, plotData):
        ax.set_title(contigID)
        ax.bar(x, y, zorder=0)
    for ax in fig.get_axes():
        ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    plt.close()

def histo_regions(binDict, regions, binSize, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        binDict -- a dictionary with structure like:
                   {
                       'contig1': [binInt1, binInt2, ...],
                       'contig1': [...],
                       ...
                   }
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        binSize -- an integer value indicating the bin size to count variants
                   within
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''    
    # Parse out regions for plotting
    regions = [ region.split(":") for region in regions ]
    
    # Convert start and end values to int
    regions = [ (contigID, int(start), int(end)) for contigID, start, end in regions ]
    
    # Plot each region
    for contig, start, end in regions:
        if not contig in binDict:
            continue
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.{start}_to_{end}.histogram.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Histogram plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Figure out which bins fall within this region
        binStart = start // binSize
        binEnd = end // binSize
        regionBins = binDict[contig][binStart:binEnd+1]
        
        # Raise error if incorrect number of bins found
        if len(regionBins) != binEnd - binStart + 1:
            raise ValueError(f"Region '{contig, start, end}' only has {len(regionBins)} bins " +
                             f"but should have {binEnd - binStart + 1} bins. This probably means " +
                             "that your region coordinates are incorrect.")
        
        # Get plotting values
        x = np.arange(binStart, binEnd+1)
        y = regionBins
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        stepSize = binSize / 1000 # convert to Kbp
        ax.set_xlabel(f"Bin number ({stepSize} Kbp step and width)", fontweight="bold")
        ax.set_ylabel(f"Number of variants", fontweight="bold")
        ax.set_title(f"{contig} variant occurrence histogram", fontweight="bold")
        
        # Plot histogram
        ax.bar(x, y, zorder=0)
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def histo_per_contig(binDict, binSize, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        binDict -- a dictionary with structure like:
                   {
                       'contig1': [binInt1, binInt2, ...],
                       'contig1': [...],
                       ...
                   }
        binSize -- an integer value indicating the bin size to count variants
                   within
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''
    for contig in binDict.keys():
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.histogram.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Histogram plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Get plotting values
        y = binDict[contig]
        x = np.arange(0, len(y))
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        stepSize = binSize / 1000 # convert to Kbp
        ax.set_xlabel(f"Window number ({stepSize} Kbp step and width)", fontweight="bold")
        ax.set_ylabel(f"Number of variants", fontweight="bold")
        ax.set_title(f"{contig} variant occurrence histogram", fontweight="bold")
        
        # Plot ideogram
        ax.bar(x, y, zorder=0)
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def ideo_horizontal(binDict, lengthsDict, binSize, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        binDict -- a dictionary with structure like:
                   {
                       'contig1': [binInt1, binInt2, ...],
                       'contig1': [...],
                       ...
                   }
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        binSize -- an integer value indicating the bin size to count variants
                   within
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''
    SPACING = 0.1
    
    # Get contig ordering by length
    contigOrder = sorted(lengthsDict.keys(), key=lambda x: lengthsDict[x]) # shortest to longest; gets reversed during plotting
    
    # Derive our output file name and skip if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    fileOut = os.path.join(outputDirectory, f"onePlot.ideogram.{fileSuffix}")
    if os.path.isfile(fileOut):
        raise FileExistsError(f"'onePlot.ideogram.{fileSuffix}' already found in output directory")
    
    # Get minimum and maximum density values for whole data set
    minDensity = min([ min(values) for values in binDict.values() ])
    maxDensity = max([ max(values) for values in binDict.values() ])
    longestDensity = max([ len(values) for values in binDict.values() ])
    
    # Colour map density values
    cmap = plt.cm.viridis
    norm = matplotlib.colors.Normalize(vmin=minDensity, vmax=maxDensity)
    
    # Configure plot
    fig = plt.figure(figsize=(width, height), tight_layout=True)
    ax = plt.axes()
    
    stepSize = binSize / 1000 # convert to Kbp
    ax.set_xlabel(f"Window number ({stepSize} Kbp step and width)", fontweight="bold")
    ax.set_ylabel("") # no y-axis label
    ax.set_xlim(0, longestDensity) # longestDensity implicitly adds +1 which includes the last bin
    
    # Plot ideograms per contig
    contigLabels = []
    ongoingCount = 0.5 # this centers the contig label
    for contig in contigOrder:
        if not contig in binDict:
            continue
        
        # Get plotting values
        y = binDict[contig]
        xranges = [ (x, 1) for x in np.arange(0, len(y)) ]
        
        # Plot contig ideogram
        ax.broken_barh(xranges, (ongoingCount+SPACING, 1-(SPACING*2)), facecolors=cmap(norm(y)))
        
        # Iterate values
        contigLabels.append(contig)
        ongoingCount += 1
    ax.set_ylim(0.5+SPACING, ongoingCount-SPACING)
    
    # Indicate contig labels
    ax.set_yticks(range(1, len(contigLabels)+1))
    ax.set_yticklabels(contigLabels)
    
    # Show the colour scale legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation="vertical", label="SNP number")
    
    # Save output file
    plt.savefig(fileOut)
    plt.close()

def ideo_regions(binDict, regions, binSize, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        binDict -- a dictionary with structure like:
                   {
                       'contig1': [binInt1, binInt2, ...],
                       'contig1': [...],
                       ...
                   }
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        binSize -- an integer value indicating the bin size to count variants
                   within
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''    
    # Parse out regions for plotting
    regions = [ region.split(":") for region in regions ]
    
    # Convert start and end values to int
    regions = [ (contigID, int(start), int(end)) for contigID, start, end in regions ]
    
    # Plot each region
    for contig, start, end in regions:
        if not contig in binDict:
            continue
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.{start}_to_{end}.ideogram.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Ideogram plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Figure out which bins fall within this region
        binStart = start // binSize
        binEnd = end // binSize
        regionBins = binDict[contig][binStart:binEnd+1]
        
        # Raise error if incorrect number of bins found
        if len(regionBins) != binEnd - binStart + 1:
            raise ValueError(f"Region '{contig, start, end}' only has {len(regionBins)} bins " +
                             f"but should have {binEnd - binStart + 1} bins. This probably means " +
                             "that your region coordinates are incorrect.")
        
        # Get plotting values
        stepSize = binSize / 1000 # convert to Kbp
        xranges = [ (x, 1) for x in np.arange(binStart, binEnd+1) ]
        y = regionBins
        
        # Colour map density values
        cmap = plt.cm.viridis
        norm = matplotlib.colors.Normalize(vmin=min(y), vmax=max(y))
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Window number ({stepSize} Kbp step and width)", fontweight="bold")
        ax.get_yaxis().set_visible(False)
        ax.set_xlim(binStart, binEnd+1) # +1 to include the last bin
        ax.set_ylim(1, 2)
        
        # Plot ideogram
        ax.broken_barh(xranges, (1, 2), facecolors=cmap(norm(y)))
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="SNP number")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def ideo_per_contig(binDict, binSize, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        binDict -- a dictionary with structure like:
                   {
                       'contig1': [binInt1, binInt2, ...],
                       'contig1': [...],
                       ...
                   }
        binSize -- an integer value indicating the bin size to count variants
                   within
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''    
    for contig in binDict.keys():
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.ideogram.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Ideogram plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Get plotting values
        y = binDict[contig]
        xranges = [ (x, 1) for x in np.arange(0, len(y)) ]
        stepSize = binSize / 1000 # convert to Kbp
        
        # Colour map density values
        cmap = plt.cm.viridis
        norm = matplotlib.colors.Normalize(vmin=min(y), vmax=max(y))
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Window number ({stepSize} Kbp step and width)", fontweight="bold")
        ax.set_ylabel(contig, fontweight="bold")
        plt.tick_params(left = False, labelleft = False) 
        ax.set_xlim(0, len(y)) # len(y) implicitly adds +1 which includes the last bin
        ax.set_ylim(1, 2)
        
        # Plot ideogram
        ax.broken_barh(xranges, (1, 2), facecolors=cmap(norm(y)))
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="SNP number")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def gene_regions(tallyDict, gff3, regions, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        tallyDict -- a dictionary with structure like:
                     {
                         "geneID1": [snpCount, geneLength],
                         "geneID2": [snpCount, geneLength],
                         ...
                     }
        gff3 -- a ZS_GFF3IO.GFF3 object
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''
    SPACING = 0.5
    
    # Parse out regions for plotting
    regions = [ region.split(":") for region in regions ]
    
    # Convert start and end values to int
    regions = [ (contigID, int(start), int(end)) for contigID, start, end in regions ]
    
    # Plot each region
    for contig, start, end in regions:
        if not contig in tallyDict:
            continue
        
        # Get gene IDs within this region
        geneFeatures = gff3.ncls_finder(start, end, "contig", contig)
        
        geneIDs = [
            [geneFeature.ID, geneFeature.start, geneFeature.end]
            for geneFeature in geneFeatures
        ]
        geneIDs.sort(key=lambda x: (x[1], x[2]))
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.{start}_to_{end}.genegram.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Genegram plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Get plotting values
        xranges, y = [] , []
        ongoingCount = 1
        for geneID, geneStart, geneEnd in geneIDs:
            snpCount, geneLength = tallyDict[geneID]
            density = snpCount / geneLength
            y.append(density)
            
            xrange = [ongoingCount-SPACING, 1]
            xranges.append(xrange)
            ongoingCount += 1
        
        # Colour map density values
        cmap = plt.cm.viridis
        norm = matplotlib.colors.Normalize(vmin=min(y), vmax=max(y))
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Gene number (in sequence within region)", fontweight="bold")
        ax.set_ylabel(f"{contig} ({start} to {end}bp)", fontweight="bold")
        plt.tick_params(left = False, labelleft = False) 
        ax.set_xlim(1-SPACING, ongoingCount+SPACING)
        #ax.set_yticks(range(0, len(xranges)+1))
        
        # Plot genegram
        ax.broken_barh(xranges, (1, 2), facecolors=cmap(norm(y)))
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="Variant density")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def gene_per_contig(tallyDict, gff3, width, height, outputDirectory, plotPDF):
    '''
    Parameters:
        tallyDict -- a dictionary with structure like:
                     {
                         "geneID1": [snpCount, geneLength],
                         "geneID2": [snpCount, geneLength],
                         ...
                     }
        gff3 -- a ZS_GFF3IO.GFF3 object
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
    '''
    SPACING = 0.5
    
    for contig in gff3.contigs:
        # Get gene IDs along this contig
        geneIDs = [
            [geneFeature.ID, geneFeature.start, geneFeature.end]
            for geneFeature in gff3.types["gene"]
            if geneFeature.contig == contig
        ]
        geneIDs.sort(key=lambda x: (x[1], x[2]))
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.genegram.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Genegram plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Skip if no genes found
        if len(geneIDs) == 0:
            print(f"WARNING: No genes found for '{contig}'; skipping...")
            continue
        
        # Get plotting values
        xranges, y = [] , []
        ongoingCount = 1
        for geneID, geneStart, geneEnd in geneIDs:
            snpCount, geneLength = tallyDict[geneID]
            density = snpCount / geneLength
            y.append(density)
            
            xrange = [ongoingCount-SPACING, 1]
            xranges.append(xrange)
            ongoingCount += 1
        
        # Colour map density values
        cmap = plt.cm.viridis
        norm = matplotlib.colors.Normalize(vmin=min(y), vmax=max(y))
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Gene number (in sequence along chromosome)", fontweight="bold")
        ax.set_ylabel(contig, fontweight="bold")
        plt.tick_params(left = False, labelleft = False) 
        ax.set_xlim(1-SPACING, ongoingCount+SPACING)
        
        # Plot genegram
        ax.broken_barh(xranges, (1, 2), facecolors=cmap(norm(y)))
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="Variant density")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def main():
    usage = """%(prog)s receives a VCF and genome FASTA to creates SNP density plots per
    chromosome. It can produce line plots, histograms, and ideograms of SNP density.
    It also has an ideogram-like plot for SNP density per gene referred to as a 'genegram'.
    
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
    
    geneparser = subparsers.add_parser("genegram",
                                       aliases=["gene"],
                                       parents=[p],
                                       add_help=False,
                                       help="Create broken bar plots of SNP density per gene")
    geneparser.set_defaults(func=genemain)
    
    # Line-subparser arguments
    lineparser.add_argument("--onePlot", dest="onePlot",
                            required=False,
                            action="store_true",
                            help="""Optionally, provide this flag if you want a single plot to be
                            produced with all chromosomes""",
                            default=False)
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
    histoparser.add_argument("--onePlot", dest="onePlot",
                             required=False,
                             action="store_true",
                             help="""Optionally, provide this flag if you want a single plot to be
                             produced with all chromosomes""",
                             default=False)
    histoparser.add_argument("--binSize", dest="binSize",
                             type=int,
                             required=False,
                             help="""Optionally, specify the bin size to count variants
                             within (default=10000)""",
                             default=10000)
    
    # Ideogram-subparser arguments
    ideoparser.add_argument("--onePlot", dest="onePlot",
                            required=False,
                            action="store_true",
                            help="""Optionally, provide this flag if you want a single plot to be
                            produced with all chromosomes""",
                            default=False)
    ideoparser.add_argument("--binSize", dest="binSize",
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
    elif args.mode in ["ideogram", "ideo"]:
        ideomain(args, lengthsDict)
    elif args.mode in ["genegram", "gene"]:
        genemain(args)
    
    print("Program completed successfully!")

def linemain(args, lengthsDict):
    # Figure out what our pickle file should be called
    hashString = f"{args.vcfFile}" + \
                 f"{'.wbo' if args.weightByOccurrence else ''}" + \
                 f"{'.sma' if args.skipMonoallelic else ''}.pkl" + \
                 f"{args.regions}{args.filter}"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Get SNP density from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            dotsX, dotsY = pickle.load(fileIn)
    else:
        dotsX, dotsY = get_vcf_variants(args.vcfFile, lengthsDict, args.windowSize,
                                        args.weightByOccurrence, args.skipMonoallelic,
                                        args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump([dotsX, dotsY], fileOut)
    
    # Create plots
    if args.onePlot:
        line_horizontal(dotsX, dotsY, args.wmaSize,
                            args.width, args.height,
                            args.outputDirectory, args.plotPDF,
                            args.showDots, args.linewidth)
    elif args.regions != []:
        line_regions(dotsX, dotsY, args.regions, args.wmaSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF,
                         args.showDots, args.linewidth)
    else:
        line_per_contig(dotsX, dotsY, args.wmaSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.showDots, args.linewidth)

def histomain(args, lengthsDict):
    # Figure out what our pickle file should be called
    hashString = f"{args.vcfFile}" + \
                 f"{'.wbo' if args.weightByOccurrence else ''}" + \
                 f"{'.sma' if args.skipMonoallelic else ''}.pkl" + \
                 f"{args.regions}{args.filter}"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Get SNP density from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            binDict = pickle.load(fileIn)
    else:
        binDict = bin_vcf_variants(args.vcfFile, lengthsDict, args.binSize,
                                   args.weightByOccurrence, args.skipMonoallelic,
                                   args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump(binDict, fileOut)
    
    # Create plots
    if args.onePlot:
        histo_horizontal(binDict, args.binSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF)
    elif args.regions != []:
        histo_regions(binDict, args.regions, args.binSize,
                      args.width, args.height,
                      args.outputDirectory, args.plotPDF)
    else:
        histo_per_contig(binDict, args.binSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF)

def ideomain(args, lengthsDict):
    # Figure out what our pickle file should be called
    hashString = f"{args.vcfFile}" + \
                 f"{'.wbo' if args.weightByOccurrence else ''}" + \
                 f"{'.sma' if args.skipMonoallelic else ''}.pkl" + \
                 f"{args.regions}{args.filter}"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Get SNP density from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            binDict = pickle.load(fileIn)
    else:
        binDict = bin_vcf_variants(args.vcfFile, lengthsDict, args.binSize,
                                   args.weightByOccurrence, args.skipMonoallelic,
                                   args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump(binDict, fileOut)
    
    # Create plots
    if args.onePlot:
        ideo_horizontal(binDict, lengthsDict, args.binSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF)
    elif args.regions != []:
        ideo_regions(binDict, args.regions, args.binSize,
                     args.width, args.height,
                     args.outputDirectory, args.plotPDF)
    else:
        ideo_per_contig(binDict, args.binSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF)

def genemain(args):
    # Figure out what our pickle file should be called
    hashString = f"{args.vcfFile}" + \
                 f"{'.wbo' if args.weightByOccurrence else ''}" + \
                 f"{'.sma' if args.skipMonoallelic else ''}.pkl" + \
                 f"{args.regions}{args.filter}" + \
                 f"{args.gff3File}" + \
                 f"{'.cds' if args.onlyCDS else ''}"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Parse GFF3 with NCLS indexing
    gff3 = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    gff3.create_ncls_index(typeToIndex="gene")
    
    # Get SNP tally from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            tallyDict = pickle.load(fileIn)
    else:
        tallyDict = tally_variants_within_features(args.vcfFile, gff3,
                                                   args.onlyCDS, args.weightByOccurrence,
                                                   args.skipMonoallelic, args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump(tallyDict, fileOut)
    
    # Create plots
    if args.regions != []:
        gene_regions(tallyDict, gff3, args.regions,
                     args.width, args.height,
                     args.outputDirectory, args.plotPDF)
    else:
        gene_per_contig(tallyDict, gff3,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF)

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
