#! python3
# genome_annotation_circlos.py
# Script to plot genome annotation data in a circular plot using pyCirclize

import os, argparse, math
import numpy as np
from Bio import SeqIO, SeqFeature
from matplotlib.patches import Patch

from pycirclize import Circos
from pycirclize.parser import Gff

colourPalette = ['#377eb8', '#ff7f00', '#4daf4a',
                 '#f781bf', '#a65628', '#984ea3'] #,
                 #'#999999', '#e41a1c', '#dede00']

def validate_args(args):
    # Validate input files
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f"-f '{args.genomeFasta}' is not a file!")
    if not os.path.isfile(args.genomeGFF3):
        raise FileNotFoundError(f"-g '{args.genomeGFF3}' is not a file!")
    if args.repeatMaskerOut != None:
        if not os.path.isfile(args.repeatMaskerOut):
            raise FileNotFoundError(f"--repeats '{args.repeatMaskerOut}' is not a file!")
        if args.repeatFamilies == []:
            raise ValueError(f"--repeats specified but no --families provided!")
    
    # Validate that number of repeat families is supported
    if (len(args.repeatFamilies)-1) > len(colourPalette):
        raise ValueError(f"Number of repeat families ({len(args.repeatFamilies)}) exceeds the number " + 
                         f"of supported colours ({len(colourPalette)})")
    
    # Validate numeric values
    if args.binSize < -1:
        raise ValueError(f"--binSize '{args.binSize}' must be greater than or equal to -1")
    elif args.binSize == 0:
        raise ValueError(f"If setting --binSize, you must provide a positive integer")
    
    # Validate output file
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o output file '{args.outputFileName}' already exists!")
    args.outputFileName = os.path.abspath(args.outputFileName)
    
    if not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"-o parent directory '{os.path.dirname(args.outputFileName)}' does not exist!")

def bin_values(features, contigSize, binSize):
    '''
    Parameters:
        features -- a list of SeqFeature objects
        contigSize -- an integer value indicating the length of the contig
        binSize -- an integer value indicating the size of the bins/windows
    Returns:
        x -- a numpy array of bin midpoints
        y -- a numpy array of the number of features in each bin
    '''
    if binSize > contigSize:
        raise ValueError(f"binSize '{binSize}' is greater than contigSize '{contigSize}'")
    
    # Get the cleanest division of the contig into the bin size
    remainder = (contigSize / binSize) % 1
    leftover = (remainder * binSize) / (contigSize // binSize)
    adjustedBinSize = binSize + math.ceil(leftover)
    
    y = np.array([
        0 for _ in range(math.ceil((contigSize - 1) / adjustedBinSize)) # -1 for 0-based indexing
    ])
    for feature in features:
        startBin = feature.location.start // adjustedBinSize
        endBin = feature.location.end // adjustedBinSize
        for binIndex in range(startBin, endBin+1):
            y[binIndex] += 1
    return np.arange(0, contigSize, adjustedBinSize), y

def plot_gene_features(seqid2features, sector, colour, currentPosition,
                       showOutline, TRACK_HEIGHT):
    track = sector.add_track((currentPosition-TRACK_HEIGHT, currentPosition))
    if showOutline:
        track.axis()
     
    for feature in seqid2features[sector.name]:
        if feature.type == "gene":
            track.genomic_features([feature], fc=colour)
    return track

def plot_gene_density(seqid2features, sector, colour, currentPosition, contigSize, binSize,
                      showOutline, TRACK_HEIGHT):
    x, y = bin_values([ feature for feature in seqid2features[sector.name] if feature.type == "gene" ],
                      contigSize, binSize)
    
    track = sector.add_track((currentPosition-TRACK_HEIGHT, currentPosition))
    if showOutline:
        track.axis()
    
    track.bar(x, y, vmax=max(y), width=binSize, fc=colour, align="edge")
    return track

def plot_repeat_features(seqid2features, sector, colour, currentPosition, repeatType,
                         showOutline,TRACK_HEIGHT):
    track = sector.add_track((currentPosition-TRACK_HEIGHT, currentPosition))
    if showOutline:
        track.axis()
     
    for feature in seqid2repeats[seqid]:
        if feature.qualifiers["class"] == repeatType or feature.qualifiers["family"] == repeatType:
            track.genomic_features([feature], fc=colour)
    return track

def plot_repeat_density(seqid2features, sector, colour, currentPosition, repeatType,
                        contigSize, binSize, showOutline, TRACK_HEIGHT):
    x, y = bin_values([ feature for feature in seqid2features[sector.name]
                       if feature.qualifiers["class"] == repeatType or feature.qualifiers["family"] == repeatType ],
                      contigSize, binSize)
    
    track = sector.add_track((currentPosition-TRACK_HEIGHT, currentPosition))
    if showOutline:
        track.axis()
    
    track.bar(x, y, vmax=max(y), width=binSize, fc=colour, align="edge")
    return track

def main():
    usage = """%(prog)s produces a circos plot representation of a genome and its gene annotations
    with optional repeat annotation through RepeatMasker output.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Required arguments
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-g", dest="genomeGFF3",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify the location to write the output file; this must
                   end with '.pdf', '.png', or '.svg'""")
    # Optional arguments
    p.add_argument("--legend", dest="legendLocation",
                   required=False,
                   choices=["top", "centre", "left"],
                   help="""Optionally, specify the location of the legend as being at
                   the 'top', 'centre', or 'left'; default is 'centre'""",
                   default="centre")
    p.add_argument("--binSize", dest="binSize",
                   required=False,
                   type=int,
                   help="""Optionally, specify the bin size for represnting features
                   by their density; default is -1 which means no binning occurs
                   and each feature is represented individually""",
                   default=-1)
    p.add_argument("--repeats", dest="repeatMaskerOut",
                   required=False,
                   help="Optionally, specify the location of the RepeatMasker .out file")
    p.add_argument("--families", dest="repeatFamilies",
                   required=False,
                   nargs="+",
                   help="""If you are plotting RepeatMasker output, specify
                   the repeat classes and/or families to plot.""",
                   default = [])
    p.add_argument("--outline", dest="sectorOutline",
                   required=False,
                   action="store_true",
                   help="Provide this flag to outline each sector in the plot",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    seqid2size = { record.id:len(record) for record in genomeRecords }
    
    # Raise error if no contigs are found in genome FASTA
    if seqid2size == {}:
        raise ValueError(f"No contigs found in genome FASTA '{args.genomeFasta}'; is it actually a FASTA file?")
    
    # Parse GFF3 file
    parser = Gff(args.genomeGFF3)
    seqid2features = parser.get_seqid2features(feature_type=None)
    
    # Parse RepeatMasker output
    seqid2repeats = {}
    if args.repeatMaskerOut != None:
        with open(args.repeatMaskerOut, "r") as fileIn:
            for line in fileIn:
                # Skip header lines
                sl = line.strip().split()
                if sl == [] or sl[0] == "SW" or sl[0] == "score":
                    continue
                
                # Parse relevant info
                seqid, start, end, repeatID, classFam = *sl[4:7], sl[9], sl[10]
                if "/" in sl[10]:
                    repeatClass, repeatFam = classFam.split("/")
                else:
                    repeatClass, repeatFam = sl[10], sl[10]
                if not seqid in seqid2size:
                    continue
                
                # Store repeat info if requested
                if repeatClass in args.repeatFamilies or repeatFam in args.repeatFamilies:
                    if seqid not in seqid2repeats:
                        seqid2repeats[seqid] = []
                    
                    # Format and store a SeqFeature
                    feature = SeqFeature.SeqFeature(
                        location=SeqFeature.FeatureLocation(int(start), int(end)),
                        type="repeat",
                        strand=None,
                        id=repeatID,
                        qualifiers={"class": repeatClass, "family": repeatFam}
                    )
                    
                    seqid2repeats[seqid].append(feature)
    
    # Initialise circos instance
    space = 0 if len(seqid2size) == 1 else 2
    if args.legendLocation == "top":
        circos = Circos(seqid2size, space=space,
                        start=10 if len(args.repeatFamilies) < 5 else 20,
                        end=350 if len(args.repeatFamilies) < 5 else 340,
                        endspace=False)
    elif args.legendLocation == "centre":
        circos = Circos(seqid2size, space=space)
    elif args.legendLocation == "left":
        circos = Circos(seqid2size, space=space,
                        start=-80,
                        end=260,
                        endspace=False)
    
    # Plot each chromosome / sector
    TRACK_GAP = 2
    OUTER_HEIGHT = 0.3
    TRACK_HEIGHT = 10 if len(args.repeatFamilies) < 5 else 8
    for sector in circos.sectors:
        handles = [] # it is okay if this gets reset each loop, we just need the last one
        currentPosition = 95
        trackIndex = 0
        
        # Plot outer chromosome track
        sector.text(sector.name, size=10)
        outer_track = sector.add_track((currentPosition-OUTER_HEIGHT, currentPosition))
        outer_track.axis(fc="black")
        major_interval = 10000000
        minor_interval = int(major_interval / 10)
        if sector.size > minor_interval:
            outer_track.xticks_by_interval(major_interval, label_formatter=lambda v: f"{v / 1000000:.0f} Mb")
            outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)
        currentPosition -= (OUTER_HEIGHT + TRACK_GAP)
        
        # Plot gene track data
        if args.binSize == -1:
            geneTrack = plot_gene_features(seqid2features, sector, colourPalette[trackIndex], currentPosition,
                                           args.sectorOutline, TRACK_HEIGHT)
        else:
            geneTrack = plot_gene_density(seqid2features, sector, colourPalette[trackIndex], currentPosition,
                                          seqid2size[sector.name], args.binSize, args.sectorOutline, TRACK_HEIGHT)
        currentPosition -= (TRACK_HEIGHT + TRACK_GAP)
        
        # Plot gene track legend
        if sector.name == circos.sectors[0].name:
            if args.legendLocation == "top":
                circos.text("Gene", r=geneTrack.r_center, color=colourPalette[trackIndex])
            elif args.legendLocation == "left":
                circos.text("Gene", r=geneTrack.r_center, deg=-90, color=colourPalette[trackIndex])
        handles.append(Patch(color=colourPalette[trackIndex], label="Gene"))
        trackIndex += 1
        
        # Plot each repeat track
        if args.repeatMaskerOut != None:
            for repeatType in args.repeatFamilies:
                # Plot repeat track data
                if args.binSize == -1:
                    repeatTrack = plot_repeat_features(seqid2features, sector, colourPalette[trackIndex], currentPosition,
                                                       repeatType, args.sectorOutline, TRACK_HEIGHT)
                else:
                    repeatTrack = plot_repeat_density(seqid2repeats, sector, colourPalette[trackIndex], currentPosition,
                                                      repeatType, seqid2size[sector.name], args.binSize,
                                                      args.sectorOutline, TRACK_HEIGHT)
                currentPosition -= (TRACK_HEIGHT + TRACK_GAP)
                
                # Plot repeat track legend
                if sector.name == circos.sectors[0].name:
                    if args.legendLocation == "top":
                        circos.text(repeatType, r=repeatTrack.r_center, color=colourPalette[trackIndex])
                    elif args.legendLocation == "left":
                        circos.text(repeatType, r=repeatTrack.r_center, deg=-90, color=colourPalette[trackIndex])
                handles.append(Patch(color=colourPalette[trackIndex], label=repeatType))
                trackIndex += 1
    
    # Generate the plot
    fig = circos.plotfig()
    if args.legendLocation == "centre":
        _ = circos.ax.legend(
            handles=handles,
            bbox_to_anchor=(0.5, 0.5),
            loc="center",
            ncols=2,
        )
    
    fig.savefig(args.outputFileName)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
