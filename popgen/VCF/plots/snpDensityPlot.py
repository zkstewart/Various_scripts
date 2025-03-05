#! python3
# snpDensityPlot.py
# Script to create visualisations of SNP density
# for assessing hypotheses of SNP distribution
# across chromosomes.

import os, argparse, math, sys, re, pickle
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Bio import SeqIO
from hashlib import md5
from contextlib import nullcontext

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_GFF3IO, ZS_VCFIO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError(f'Unable to locate the input VCF file ({args.vcfFile})')
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f'Unable to locate the input genome FASTA file ({args.genomeFasta})')
    if hasattr(args, "gff3File") and args.gff3File != None:
        if not os.path.isfile(args.gff3File):
            raise FileNotFoundError(f"Unable to locate the input GFF3 file '{args.gff3File}'")
        else:
            args.gff3Obj = ZS_GFF3IO.GFF3(args.gff3File) # parsing now to raise errors early
            args.gff3Obj.create_ncls_index("gene")
    else:
        args.gff3Obj = None
    
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
    
    if args.mode == "line":
        if not args.weightByOccurrence:
            raise ValueError("You must enable the --weightByOccurrence flag for line mode since " +
                             "a line plot would otherwise have no meaning.")
    
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
    elif os.path.exists(args.outputDirectory):
        raise ValueError(f"The specified output directory ({args.outputDirectory}) is not a directory!")
    else:
        os.makedirs(args.outputDirectory, exist_ok=True)
        print(f"Output directory created at '{args.outputDirectory}' as part of argument validation.")
    
    # Handle colour palette
    if args.colourMap == "viridis":
        args.cm = plt.cm.viridis
    elif args.colourMap == "Greys":
        args.cm = plt.cm.Greys
    elif args.colourMap == "GnBu":
        args.cm = plt.cm.GnBu
    elif args.colourMap == "RdBu":
        args.cm = plt.cm.RdBu
    else:
        raise ValueError("Invalid colour map specified")
    
    # Handle mode-specific arguments
    if args.mode == "line":
        if args.wmaSize < 1:
            raise ValueError("wmaSize must be a positive integer")
        if args.lineWidth < 1:
            raise ValueError("lineWidth must be a positive integer")
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

def WMA(s, period):
    """
    See https://stackoverflow.com/questions/74518386/improving-weighted-moving-average-performance
    
    Parameters:
        s -- a numpy array of values to smooth
        period -- an integer value indicating the number of previous values to consider
                  during weighted moving average calculation
    Returns:
        sw -- a pandas Series of the smoothed values
    """
    w = np.arange(period)+1
    w_s = w.sum()
    
    try:
        swv = np.lib.stride_tricks.sliding_window_view(s.flatten(), window_shape=period)
    except ValueError:
        "Less data points than period size causes this error"
        return None
    sw = (swv * w).sum(axis=1) / w_s
    
    # Need to now return it as a normal series
    sw = np.concatenate((np.full(period - 1, np.nan), sw))
    try:
        sw[0:period] = sw[period] # set first n=period values to be same as first smoothed value
    except:
        "len(sw)==1 causes this error"
        return None
    return pd.Series(sw)

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
            # Skip non-mRNA genes
            if not hasattr(geneFeature, "mRNA"):
                continue
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
                
                # Skip genes that were skipped in the initial tally
                "Can occur if onlyCDS is enabled and the gene has no mRNA children"
                if not geneID in tallyDict:
                    continue
                
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

def count_variants_along_features(vcfFile, gff3, weightByOccurrence=False,
                                   skipMonoallelic=False, filterSamples=[]):
    '''
    Counts variants by position along a feature, rather than simply tallying the number
    of variants per feature.
    
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        gff3 -- a ZS_GFF3IO.GFF3 object
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
        
        # Skip non-mRNA genes
        if not hasattr(geneFeature, "mRNA"):
            continue
        
        longestMrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(geneFeature)
        coords, zeros = [], []
        for cdsFeature in longestMrnaFeature.CDS:
            coords.append((cdsFeature.start, cdsFeature.end))
            zeros.append(np.zeros(cdsFeature.end - cdsFeature.start + 1, np.uint)) # value will only be a positive integer
        zeros = np.concatenate(zeros)
        
        tallyDict[geneID] = [coords, zeros]
    
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
                coords, zeros = tallyDict[geneID]
                
                # Narrow down to variants within the CDS region
                isWithinCDS = any([ start <= pos <= end for start, end in coords ])
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
                
                # Locate the variant position in the zeros array
                ongoingCount = 0
                for start, end in coords:
                    if start <= pos <= end:
                        ongoingCount += pos - start
                        zeros[ongoingCount] += snpCount
                        break
                    ongoingCount += end - start + 1
                
                # Add to the tally
                zeros[ongoingCount] += snpCount
    return tallyDict

def bin_vcf_variants(vcfFile, windowSize=100000, weightByOccurrence=False,
                    skipMonoallelic=False, filterSamples=[]):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
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
    lengthsDict = {}
    with ZS_VCFIO.open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # replace quotations may help with Excel files
            
            # Handle header lines
            if line.startswith("##contig"):
                contigID, contigLength = sl[0].split("ID=")[1].split(",length=")
                lengthsDict[contigID] = int(contigLength.rstrip(">"))
            elif line.startswith("#CHROM"):
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
    dotsX, dotsY = {}, {}
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

def line_horizontal(dotsX, dotsY, lengthsDict, wmaSize, width, height,
                    outputDirectory, plotPDF, showDots, lineWidth, createTSV):
    '''
    Parameters:
        dotsX -- a dict pairing chromosome IDs (keys) to a list of integers indicating
                 the position for each SNP (values)
        dotsY -- a dict pairing chromosome IDs (keys) to a list of integers indicating the
                 count for each SNP (values)
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        wmaSize -- an integer value indicating the window size for the weighted moving average
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
        showDots -- a boolean flag indicating whether to show the dots on the plot
        lineWidth -- an integer value indicating the width of the line to plot
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
    '''
    # Get contig ordering by length
    contigOrder = sorted(lengthsDict.keys(), key=lambda x: lengthsDict[x], reverse=True) # longest to shortest
    
    # Derive our output file name and skip if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    fileOut = os.path.join(outputDirectory, f"onePlot.line.{fileSuffix}")
    if os.path.isfile(fileOut):
        raise FileExistsError(f"'onePlot.line.{fileSuffix}' already found in output directory")
    
    # Get each contigs' plot data
    plotData = []
    for contig in contigOrder:
        # Skip if we found no SNPs on this contig
        if not contig in dotsY:
            print(f"WARNING: '{contig}' is in the genome FASTA but has no SNPs associated " +
                    "with it; skipping...")
            continue
        
        # Get plotting values
        x = np.array(dotsX[contig]) / 1000000 # convert to Mbp
        y = np.array(dotsY[contig])
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contig}' has too few data points to smooth; skipping...")
            continue
        
        # Store dot values
        plotData.append([contig, x, y, smoothedY])
    
    # Produce the figure axes
    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    gs = fig.add_gridspec(1, len(plotData), hspace=0)
    axes = gs.subplots(sharey='row')
    
    ## Set the figure title
    fig.suptitle(f"Variant occurrence plot", fontweight="bold")
    fig.supxlabel(f"Chromosomal position (Mbp)", fontweight="bold")
    fig.supylabel(f"Weighted moving average of variant occurrence (WMA size = {wmaSize})", fontweight="bold")
    
    # Plot the data into each axis
    with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
        # Write TSV header if applicable
        if createTSV:
            fileOutTSV.write("contigID\tposition\tvariant_occurrence\tsmoothed_variant_occurrence\n")
        
        # Plot each contig
        for ax, (contigID, x, y, smoothedY) in zip(axes, plotData):
            # Set plot title
            ax.set_title(contigID)
            
            # Plot dots (if applicable)
            if showDots:
                ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
            
            # Plot line
            ax.plot(x, smoothedY, zorder=1, linewidth=lineWidth)
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal, smoothedYVal in zip(x*1000000, y, smoothedY): # convert back to bp
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\t{smoothedYVal}\n")
        
        for ax in fig.get_axes():
            ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    plt.close()

def line_regions(dotsX, dotsY, regions, wmaSize, width, height, outputDirectory,
                 plotPDF, showDots, lineWidth, createTSV):
    '''
    Parameters:
        dotsX -- a dict pairing chromosome IDs (keys) to a list of integers indicating
                 the position for each SNP (values)
        dotsY -- a dict pairing chromosome IDs (keys) to a list of integers indicating the
                 count for each SNP (values)
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        wmaSize -- an integer value indicating the window size for the weighted moving average
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
        showDots -- a boolean flag indicating whether to show the dots on the plot
        lineWidth -- an integer value indicating the width of the line to plot
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
    '''
    # Parse out regions for plotting
    regions = [ region.split(":") for region in regions ]
    
    # Convert start and end values to int
    regions = [ (contigID, int(start), int(end)) for contigID, start, end in regions ]
    
    # Plot each region
    for contig, start, end in regions:
        if not contig in dotsX:
            continue
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.{start}_to_{end}.line.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Line plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Get values within this region
        regionValues = [ [x, y] for x, y in zip(dotsX[contig], dotsY[contig]) if x >= start and x <= end ]
        x = np.array([ x / 1000000 for x, y in regionValues ]) # convert to Mbp
        y = np.array([ y for x, y in regionValues ])
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contig}' has too few data points to smooth with WMA size of {wmaSize}; skipping...")
            continue
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position (Mbp)", fontweight="bold")
        ax.set_ylabel(f"Weighted moving average of variant occurrence (WMA size = {wmaSize})", fontweight="bold")
        ax.set_title(contig, fontweight="bold")
        
        # Plot dots (if applicable)
        if showDots:
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\tposition\tvariant_occurrence\tsmoothed_variant_occurrence\n")
                
            # Plot the contig
            ax.plot(x, smoothedY, zorder=1, linewidth=lineWidth)
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal, smoothedYVal in zip(x*1000000, y, smoothedY):
                    fileOutTSV.write(f"{contig}\t{xVal}\t{yVal}\t{smoothedYVal}\n")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def line_per_contig(dotsX, dotsY, wmaSize, width, height, outputDirectory,
                    plotPDF, showDots, lineWidth, createTSV):
    '''
    Parameters:
        dotsX -- a dict pairing chromosome IDs (keys) to a list of integers indicating
                 the position for each SNP (values)
        dotsY -- a dict pairing chromosome IDs (keys) to a list of integers indicating the
                 count for each SNP (values)
        wmaSize -- an integer value indicating the window size for the weighted moving average
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
        showDots -- a boolean flag indicating whether to show the dots on the plot
        lineWidth -- an integer value indicating the width of the line to plot
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
    '''
    for contig, x in dotsX.items():
        
        # Derive our output file name and skip if already existing
        fileSuffix = "pdf" if plotPDF else "png"
        fileOut = os.path.join(outputDirectory, f"{contig}.line.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Line plot for '{contig}' already found in output directory; skipping...")
            continue
        
        # Get plotting values
        x = np.array(x) / 1000000 # convert to Mbp
        y = np.array(dotsY[contig])
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contig}' has too few data points to smooth with WMA size of {wmaSize}; skipping...")
            continue
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Chromosomal position (Mbp)", fontweight="bold")
        ax.set_ylabel(f"Weighted moving average of variant occurrence (WMA size = {wmaSize})", fontweight="bold")
        ax.set_title(contig, fontweight="bold")
        
        # Plot dots (if applicable)
        if showDots:
            ax.scatter(x, y, color="red", s=3, alpha=0.5, zorder=0)
        
        # Plot line
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\tposition\tvariant_occurrence\tsmoothed_variant_occurrence\n")
            
            # Plot the contig
            ax.plot(x, smoothedY, zorder=1, linewidth=lineWidth)
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal, smoothedYVal in zip(x*1000000, y, smoothedY):
                    fileOutTSV.write(f"{contig}\t{xVal}\t{yVal}\t{smoothedYVal}\n")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def histo_horizontal(binDict, binSize, width, height, outputDirectory,
                     plotPDF, createTSV):
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
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
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
    with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
        # Write TSV header if applicable
        if createTSV:
            fileOutTSV.write("contigID\tbin_number\tvariant_occurrence\n")
        
        # Plot each contig
        for ax, (contigID, x, y) in zip(axes, plotData):
            ax.set_title(contigID)
            ax.bar(x, y, zorder=0)
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal in zip(x, y):
                    fileOutTSV.write(f"{contigID}\t{xVal}\t{yVal}\n")
        for ax in fig.get_axes():
            ax.label_outer()
    
    # Save output file
    plt.savefig(fileOut)
    plt.close()

def histo_regions(binDict, regions, binSize, width, height, outputDirectory,
                  plotPDF, createTSV):
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
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
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
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\tbin_number\tvariant_occurrence\n")
            
            # Plot the region
            ax.bar(x, y, zorder=0)
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal in zip(x, y):
                    fileOutTSV.write(f"{contig}\t{xVal}\t{yVal}\n")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def histo_per_contig(binDict, binSize, width, height, outputDirectory,
                     plotPDF, createTSV):
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
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
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
        
        # Plot histogram
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\tbin_number\tvariant_occurrence\n")
            
            # Plot the contig
            ax.bar(x, y, zorder=0)
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal in zip(x, y):
                    fileOutTSV.write(f"{contig}\t{xVal}\t{yVal}\n")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def ideo_horizontal(binDict, lengthsDict, binSize, width, height, outputDirectory,
                    plotPDF, createTSV, gff3Obj, cmap):
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
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
        gff3Obj -- an instance of the GFF3 class from ZS_GFF3IO.py to annotate gene
                   locations on the plot OR None
        cmap -- a matplotlib colour map to use for the plot
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
    
    with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
        # Write TSV header if applicable
        if createTSV:
            fileOutTSV.write("contigID\twindow_number\tvariant_occurrence\n")
        
        # Plot each contig
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
            
            # Plot gene locations if applicable
            if gff3Obj is not None:
                geneFeatures = gff3Obj.ncls_finder(0, lengthsDict[contig], "contig", contig)
                mrnaFeatures = [
                    ZS_GFF3IO.GFF3.longest_isoform(geneFeature)
                    for geneFeature in geneFeatures
                    if hasattr(geneFeature, "mRNA")
                ]
                mrnaCentres = [
                    int((mrnaFeature.start + mrnaFeature.end) / 2) / binSize # convert to fractional bin position
                    for mrnaFeature in mrnaFeatures
                ]
                ax.scatter(mrnaCentres, [ongoingCount+SPACING]*len(mrnaCentres), marker="v", color="#EA8527")
                
                ongoingCount += 1
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal in enumerate(y):
                    fileOutTSV.write(f"{contig}\t{xVal}\t{yVal}\n")
    ax.set_ylim(0.5+SPACING, ongoingCount-SPACING)
    
    # Indicate contig labels
    if gff3Obj is not None:
        ax.set_yticks(np.arange(1, ongoingCount + 0.5, 2))
    else:
        ax.set_yticks(range(1, len(contigLabels)+1))
    ax.set_yticklabels(contigLabels)
    
    # Show the colour scale legend
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation="vertical", label="SNP number")
    
    # Save output file
    plt.savefig(fileOut)
    plt.close()

def ideo_regions(binDict, lengthsDict, regions, binSize, width, height, outputDirectory,
                 plotPDF, createTSV, gff3Obj, cmap):
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
        regions -- a list of strings indicating regions to plot in greater detail with format
                   'contigID:startPos:endPos'
        binSize -- an integer value indicating the bin size to count variants
                   within
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
        cmap -- a matplotlib colour map to use for the plot
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
        norm = matplotlib.colors.Normalize(vmin=min(y), vmax=max(y))
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Window number ({stepSize} Kbp step and width)", fontweight="bold")
        ax.get_yaxis().set_visible(False)
        ax.set_xlim(binStart, binEnd+1) # +1 to include the last bin
        
        if gff3Obj is not None:
            ax.set_ylim(1, 3)
        else:
            ax.set_ylim(1, 2)
        
        # Plot ideogram
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\twindow_number\tvariant_occurrence\n")
            
            # Plot the region
            ax.broken_barh(xranges, (1, 2), facecolors=cmap(norm(y)))
            
            # Plot gene locations if applicable
            if gff3Obj is not None:
                geneFeatures = gff3Obj.ncls_finder(0, lengthsDict[contig], "contig", contig)
                mrnaFeatures = [
                    ZS_GFF3IO.GFF3.longest_isoform(geneFeature)
                    for geneFeature in geneFeatures
                    if hasattr(geneFeature, "mRNA")
                ]
                mrnaCentres = [
                    int((mrnaFeature.start + mrnaFeature.end) / 2) / binSize # convert to fractional bin position
                    for mrnaFeature in mrnaFeatures
                ]
                ax.scatter(mrnaCentres, [2]*len(mrnaCentres), marker="v", color="#EA8527")
            
            # Write TSV data if applicable
            if createTSV:
                for xRange, yVal in zip(xranges, y):
                    fileOutTSV.write(f"{contig}\t{xRange[0]}\t{yVal}\n")
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="SNP number")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def ideo_per_contig(binDict, lengthsDict, binSize, width, height, outputDirectory,
                    plotPDF, createTSV, gff3Obj, cmap):
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
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
        cmap -- a matplotlib colour map to use for the plot
    '''
    for contig in lengthsDict.keys():
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
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\twindow_number\tvariant_occurrence\n")
            
            # Plot the contig
            ax.broken_barh(xranges, (1, 1), facecolors=cmap(norm(y)))
            
            # Plot gene locations if applicable
            if gff3Obj is not None:
                geneFeatures = gff3Obj.ncls_finder(0, lengthsDict[contig], "contig", contig)
                mrnaFeatures = [
                    ZS_GFF3IO.GFF3.longest_isoform(geneFeature)
                    for geneFeature in geneFeatures
                    if hasattr(geneFeature, "mRNA")
                ]
                mrnaCentres = [
                    int((mrnaFeature.start + mrnaFeature.end) / 2) / binSize # convert to fractional bin position
                    for mrnaFeature in mrnaFeatures
                ]
                ax.scatter(mrnaCentres, [1.99]*len(mrnaCentres), marker="v", color="#EA8527", s=10**2)
            
            # Write TSV data if applicable
            if createTSV:
                for xVal, yVal in enumerate(y):
                    fileOutTSV.write(f"{contig}\t{xVal}\t{yVal}\n")
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="SNP number")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def gene_regions(tallyDict, gff3, regions, width, height, outputDirectory,
                 plotPDF, createTSV, cmap):
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
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
        cmap -- a matplotlib colour map to use for the plot
    '''
    SPACING = 0.5
    
    # Parse out regions for plotting
    regions = [ region.split(":") for region in regions ]
    
    # Convert start and end values to int
    regions = [ (contigID, int(start), int(end)) for contigID, start, end in regions ]
    
    # Plot each region
    for contig, start, end in regions:
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
            # Skip genes that were skipped in the initial tally
            "Can occur if onlyCDS is enabled and the gene has no mRNA children"
            if not geneID in tallyDict:
                continue
            
            snpCount, geneLength = tallyDict[geneID]
            density = snpCount / geneLength
            y.append(density)
            
            xrange = [ongoingCount-SPACING, 1]
            xranges.append(xrange)
            ongoingCount += 1
        
        # Colour map density values
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
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\tgene_number\tgene_id\tvariant_density\n")
            
            # Plot the region
            ax.broken_barh(xranges, (1, 2), facecolors=cmap(norm(y)))
            
            # Write TSV data if applicable
            if createTSV:
                for geneValue, xRange, yVal in zip(geneIDs, xranges, y):
                    geneID = geneValue[0]
                    xVal = int(xRange[0]) # implicitly subtracts SPACING
                    fileOutTSV.write(f"{contig}\t{xVal}\t{geneID}\t{yVal}\n")
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="Variant density")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def gene_per_contig(tallyDict, gff3, width, height, outputDirectory,
                    plotPDF, createTSV, cmap):
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
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
        cmap -- a matplotlib colour map to use for the plot
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
            # Skip genes that were skipped in the initial tally
            "Can occur if onlyCDS is enabled and the gene has no mRNA children"
            if not geneID in tallyDict:
                continue
            
            snpCount, geneLength = tallyDict[geneID]
            density = snpCount / geneLength
            y.append(density)
            
            xrange = [ongoingCount-SPACING, 1]
            xranges.append(xrange)
            ongoingCount += 1
        
        # Colour map density values
        norm = matplotlib.colors.Normalize(vmin=min(y), vmax=max(y))
        
        # Configure plot
        fig = plt.figure(figsize=(width, height), tight_layout=True)
        ax = plt.axes()
        
        ax.set_xlabel(f"Gene number (in sequence along chromosome)", fontweight="bold")
        ax.set_ylabel(contig, fontweight="bold")
        plt.tick_params(left = False, labelleft = False) 
        ax.set_xlim(1-SPACING, ongoingCount+SPACING)
        
        # Plot genegram
        with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
            # Write TSV header if applicable
            if createTSV:
                fileOutTSV.write("contigID\tgene_number\tgene_id\tvariant_density\n")
            
            # Plot the contig
            ax.broken_barh(xranges, (1, 2), facecolors=cmap(norm(y)))
            
            # Write TSV data if applicable
            if createTSV:
                for geneValue, xyVal in zip(geneIDs, enumerate(y)):
                    geneID = geneValue[0]
                    xVal, yVal = xyVal
                    fileOutTSV.write(f"{contig}\t{xVal}\t{geneID}\t{yVal}\n")
        
        # Show the colour scale legend
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="Variant density")
        
        # Save output file
        plt.savefig(fileOut)
        plt.close()

def gene_models(tallyDict, gff3, width, height, outputDirectory,
                plotPDF, createTSV, cmap):
    '''
    Parameters:
        tallyDict -- a dictionary with structure like:
                     {
                         "geneID1": [(coord1, coord2, ...), np.array([pos1Count, pos2Count, ...])],
                         "geneID2": [ ... ],
                         ...
                     }
        gff3 -- a ZS_GFF3IO.GFF3 object
        width -- an integer value indicating the width of the output plot
        height -- an integer value indicating the height of the output plot
        outputDirectory -- a string indicating the directory to write output plots to
        plotPDF -- a boolean flag indicating whether to output plots in PDF format
        createTSV -- a boolean flag indicating whether to output a TSV file of the data
        cmap -- a matplotlib colour map to use for the plot
    '''
    SPACING = 0.1
    
    # Derive our output file name and skip if already existing
    fileSuffix = "pdf" if plotPDF else "png"
    filePrefix = os.path.basename(gff3.fileLocation).rsplit('.gff3', maxsplit=1)[0]
    fileOut = os.path.join(outputDirectory, f"{filePrefix}.genemodels.{fileSuffix}")
    if os.path.isfile(fileOut):
        print(f"WARNING: Gene model plot for '{gff3.fileLocation}' already found in output directory; skipping...")
        return
    
    # Get representative gene models in this GFF3 object
    mrnaFeatures = [
        ZS_GFF3IO.GFF3.longest_isoform(geneFeature)
        for geneFeature in gff3.types["gene"]
        if hasattr(geneFeature, "mRNA")
    ]
    
    # Raise error if no genes found
    if len(mrnaFeatures) == 0:
        raise ValueError(f"No genes with mRNA subfeatures found in provided GFF3 file!")
    
    # Sort gene models by length
    mrnaFeatures.sort(key=lambda x: len(tallyDict[x.Parent][1]))
    
    # Get the xlim and colour map values for the plot
    xlimMax = 0
    cmapMin = np.inf
    cmapMax = 0
    for mrnaFeature in mrnaFeatures:
        coords, zeros = tallyDict[mrnaFeature.Parent]
        xlimMax = max(xlimMax, len(zeros))
        cmapMin = min(cmapMin, min(zeros))
        cmapMax = max(cmapMax, max(zeros))
    
    # Colour map density values
    norm = matplotlib.colors.Normalize(vmin=cmapMin, vmax=cmapMax)
    
    # Configure plot
    fig = plt.figure(figsize=(width, height), tight_layout=True)
    ax = plt.axes()
    
    # Set xlim and xlabel
    ax.set_xlim(0.5, xlimMax+0.5) # 1-based indexing; position centres on integers, stretch 0.5 either side
    ax.set_xlabel(f"Nucleotide position (bp)", fontweight="bold")
    
    # Plot gene models
    geneLabels = []
    with open(fileOut.replace(f".{fileSuffix}", ".tsv"), "w") if createTSV else nullcontext() as fileOutTSV:
        # Write TSV header if applicable
        if createTSV:
            fileOutTSV.write("gene_id\tposition\tvariant_density\n")
        
        # Iterate over each gene model
        for colNum, mrnaFeature in enumerate(mrnaFeatures):
            # Get plotting values
            coords, zeros = tallyDict[mrnaFeature.Parent]
            xranges = [ (x+0.5, 1) for x in np.arange(0, len(zeros)) ]
            
            # Plot gene gram
            ax.broken_barh(xranges, (colNum+SPACING, 1-(SPACING*2)), facecolors=cmap(norm(zeros)))
            
            # Plot intron splice sites
            # spliceSites = []
            # ongoingCount = 0
            # for coordIndex, (start, end) in enumerate(coords):
            #     if coordIndex != len(coords) - 1:
            #         exonLength = end - start + 1
            #         spliceSites.append(ongoingCount + exonLength)
            #         ongoingCount += exonLength
            # ax.scatter(spliceSites, [colNum+1]*len(spliceSites), marker="v", color="#EA8527", s=5**2)
            
            # Store gene name
            geneLabels.append(mrnaFeature.Parent)
            
            # Write TSV data if applicable
            if createTSV:
                for position, variant_density in enumerate(zeros):
                    fileOutTSV.write(f"{mrnaFeature.Parent}\t{position+1}\t{variant_density}\n")
    ax.set_ylim(0+SPACING, colNum+1-SPACING)
    
    # Indicate contig labels
    ax.set_yticks(np.arange(0.5, colNum + 0.5 + 1))
    ax.set_yticklabels(geneLabels)
    
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
    p.add_argument("--colour", dest="colourMap",
                    required=False,
                    choices=["viridis", "Greys", "GnBu", "RdBu"],
                    help="""Optionally, specify the colour scheme to use for the plot;
                    default is 'viridis'; refer to
                    https://matplotlib.org/stable/users/explain/colors/colormaps.html
                    for examples""",
                    default="viridis")
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
    p.add_argument("--tsv", dest="createTSV",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want to have the plotted
                    data output as a TSV file""",
                    default=False)
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
    ideoparser.add_argument("--gff3", "-g", dest="gff3File",
                            required=False,
                            help="""Optionally, specify the location of an input GFF3 file
                            if you want to annotate gene locations""")
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
    geneparser.add_argument("--models", dest="geneModels",
                            required=False,
                            action="store_true",
                            help="""Optionally, provide this flag if you want to show gene model
                            diagrams indicating the position of variants; note that you should
                            use a GFF3 subset to only include models you want to plot!""",
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
        histomain(args)
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
                 f"{args.filter}_lineplot"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Get SNP density from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            dotsX, dotsY = pickle.load(fileIn)
    else:
        dotsX, dotsY = get_vcf_variants(args.vcfFile,
                                        args.weightByOccurrence, args.skipMonoallelic,
                                        args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump([dotsX, dotsY], fileOut)
    
    # Create plots
    if args.onePlot:
        line_horizontal(dotsX, dotsY, lengthsDict, args.wmaSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.showDots, args.lineWidth,
                        args.createTSV)
    elif args.regions != []:
        line_regions(dotsX, dotsY, args.regions, args.wmaSize,
                     args.width, args.height,
                     args.outputDirectory, args.plotPDF,
                     args.showDots, args.lineWidth,
                     args.createTSV)
    else:
        line_per_contig(dotsX, dotsY, args.wmaSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.showDots, args.lineWidth,
                        args.createTSV)

def histomain(args):
    # Figure out what our pickle file should be called
    hashString = f"{args.vcfFile}" + \
                 f"{'.wbo' if args.weightByOccurrence else ''}" + \
                 f"{'.sma' if args.skipMonoallelic else ''}.pkl" + \
                 f"{args.filter}"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Get SNP density from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            binDict = pickle.load(fileIn)
    else:
        binDict = bin_vcf_variants(args.vcfFile, args.binSize,
                                   args.weightByOccurrence, args.skipMonoallelic,
                                   args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump(binDict, fileOut)
    
    # Create plots
    if args.onePlot:
        histo_horizontal(binDict, args.binSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF,
                         args.createTSV)
    elif args.regions != []:
        histo_regions(binDict, args.regions, args.binSize,
                      args.width, args.height,
                      args.outputDirectory, args.plotPDF,
                      args.createTSV)
    else:
        histo_per_contig(binDict, args.binSize,
                         args.width, args.height,
                         args.outputDirectory, args.plotPDF,
                         args.createTSV)

def ideomain(args, lengthsDict):
    # Figure out what our pickle file should be called
    hashString = f"{args.vcfFile}" + \
                 f"{'.wbo' if args.weightByOccurrence else ''}" + \
                 f"{'.sma' if args.skipMonoallelic else ''}.pkl" + \
                 f"{args.filter}"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Get SNP density from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            binDict = pickle.load(fileIn)
    else:
        binDict = bin_vcf_variants(args.vcfFile, args.binSize,
                                   args.weightByOccurrence, args.skipMonoallelic,
                                   args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump(binDict, fileOut)
    
    # Create plots
    if args.onePlot:
        ideo_horizontal(binDict, lengthsDict, args.binSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.createTSV, args.gff3Obj, args.cm)
    elif args.regions != []:
        ideo_regions(binDict, lengthsDict, args.regions, args.binSize,
                     args.width, args.height,
                     args.outputDirectory, args.plotPDF,
                     args.createTSV, args.gff3Obj, args.cm)
    else:
        ideo_per_contig(binDict, lengthsDict, args.binSize,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.createTSV, args.gff3Obj, args.cm)

def genemain(args):
    # Figure out what our pickle file should be called
    hashString = f"{args.vcfFile}" + \
                 f"{'.wbo' if args.weightByOccurrence else ''}" + \
                 f"{'.sma' if args.skipMonoallelic else ''}.pkl" + \
                 f"{args.filter}" + \
                 f"{args.gff3File}" + \
                 f"{'.cds' if args.onlyCDS else ''}" + \
                 f"{'.models' if args.geneModels else ''}"
    
    hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
    pickleFile = os.path.join(args.outputDirectory, f"snpDensityPlot_{hash}.pkl")
    
    # Get SNP tally from pickle or via calculation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            tallyDict = pickle.load(fileIn)
    else:
        if args.geneModels:
            tallyDict = count_variants_along_features(args.vcfFile, args.gff3Obj,
                                                      args.weightByOccurrence,
                                                      args.skipMonoallelic, args.filter)
        else:
            tallyDict = tally_variants_within_features(args.vcfFile, args.gff3Obj,
                                                       args.onlyCDS, args.weightByOccurrence,
                                                       args.skipMonoallelic, args.filter)
        with open(pickleFile, "wb") as fileOut:
            pickle.dump(tallyDict, fileOut)
    
    # Create plots
    if args.geneModels:
        gene_models(tallyDict, args.gff3Obj,
                    args.width, args.height,
                    args.outputDirectory, args.plotPDF,
                    args.createTSV, args.cm)
    elif args.regions != []:
        gene_regions(tallyDict, args.gff3Obj, args.regions,
                     args.width, args.height,
                     args.outputDirectory, args.plotPDF,
                     args.createTSV, args.cm)
    else:
        gene_per_contig(tallyDict, args.gff3Obj,
                        args.width, args.height,
                        args.outputDirectory, args.plotPDF,
                        args.createTSV, args.cm)

if __name__ == "__main__":
    main()
