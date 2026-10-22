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

def detect_annotation_type(fileName, NUM_LINES_TO_CHECK=100):
    '''
    Sniffs the contents of an annotation file to determine whether it is a GFF3
    or a RepeatMasker .out file.
    
    Parameters:
        fileName -- a string indicating the location of the annotation file
        NUM_LINES_TO_CHECK -- an integer indicating how many non-blank lines to
                              inspect before giving up
    Returns:
        fileType -- a string with value "gff3" or "repeatmasker"
    '''
    with open(fileName, "r") as fileIn:
        numChecked = 0
        for line in fileIn:
            l = line.rstrip("\r\n")
            if l.strip() == "":
                continue
            numChecked += 1
            if numChecked > NUM_LINES_TO_CHECK:
                break
            
            # Check for GFF3 indicators
            if l.startswith("##gff-version"):
                return "gff3"
            if l.startswith("#"):
                continue
            sl = l.split("\t")
            if len(sl) == 9 and sl[3].isdigit() and sl[4].isdigit():
                return "gff3"
            
            # Check for RepeatMasker indicators
            sl = l.strip().split()
            if sl[0] == "SW" and "score" in sl:
                return "repeatmasker"
            if len(sl) >= 14 and sl[0].isdigit() and sl[5].isdigit() and sl[6].isdigit():
                return "repeatmasker"
    raise ValueError(f"-a '{fileName}' could not be recognised as a GFF3 or RepeatMasker .out file!")

def validate_args(args):
    # Validate input files
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f"-f '{args.genomeFasta}' is not a file!")
    
    args.annotationTypes = []
    for annotationFile in args.annotationFiles:
        if not os.path.isfile(annotationFile):
            raise FileNotFoundError(f"-a '{annotationFile}' is not a file!")
        args.annotationTypes.append(detect_annotation_type(annotationFile))
    
    numGFF3 = args.annotationTypes.count("gff3")
    numRepeats = args.annotationTypes.count("repeatmasker")
    if numRepeats > 0 and args.repeatFamilies == []:
        raise ValueError(f"RepeatMasker file(s) provided to -a but no --families provided!")
    if numRepeats == 0 and args.repeatFamilies != []:
        print("WARNING: --families was provided but no RepeatMasker files were given to -a; " +
              "it will be ignored")
    
    # Validate that number of tracks is supported
    args.numTracks = numGFF3 + (numRepeats * len(args.repeatFamilies))
    if args.numTracks > len(colourPalette):
        raise ValueError(f"Number of tracks to plot ({args.numTracks}) exceeds the number " + 
                         f"of supported colours ({len(colourPalette)})")
    
    # Validate that legend labels line up with the tracks
    if args.legendLabels != [] and len(args.legendLabels) != args.numTracks:
        raise ValueError(f"Number of --labels ({len(args.legendLabels)}) does not match the number " +
                         f"of tracks to plot ({args.numTracks})")
    
    # Validate numeric values
    if args.binSize < -1:
        raise ValueError(f"--binSize '{args.binSize}' must be greater than or equal to -1")
    elif args.binSize == 0:
        raise ValueError(f"If setting --binSize, you must provide a positive integer")
    if args.binSize == -1 and args.scaleMode != "sector":
        print("WARNING: --scale only applies when --binSize is set; it will be ignored")
    if not (10 < args.totalHeight <= 85):
        raise ValueError(f"--totalHeight must be a value between 10 and 85 (lower or higher gives nonsensical output)")
    
    # Validate output file
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o output file '{args.outputFileName}' already exists!")
    args.outputFileName = os.path.abspath(args.outputFileName)
    
    if not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"-o parent directory '{os.path.dirname(args.outputFileName)}' does not exist!")

def parse_gff3(fileName):
    '''
    Parameters:
        fileName -- a string indicating the location of the GFF3 file
    Returns:
        seqid2features -- a dictionary with structure like:
                          {
                              seqid1: [ SeqFeature1, SeqFeature2, ... ],
                              ...
                          }
    '''
    parser = Gff(fileName)
    return parser.get_seqid2features(feature_type=None)

def parse_repeatmasker(fileName, repeatFamilies, seqid2size):
    '''
    Parameters:
        fileName -- a string indicating the location of the RepeatMasker .out file
        repeatFamilies -- a list of strings indicating the repeat classes and/or
                          families to retain
        seqid2size -- a dictionary of contig IDs to their lengths; repeats on
                      contigs not in this dictionary are skipped
    Returns:
        seqid2repeats -- a dictionary with structure like:
                         {
                             seqid1: [ SeqFeature1, SeqFeature2, ... ],
                             ...
                         }
    '''
    seqid2repeats = {}
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            # Skip header lines
            sl = line.strip().split()
            if sl == [] or sl[0] == "SW" or sl[0] == "score":
                continue
            
            # Parse relevant info
            seqid, start, end, repeatID, classFam = *sl[4:7], sl[9], sl[10]
            if "/" in classFam:
                repeatClass, repeatFam = classFam.split("/", 1)
            else:
                repeatClass, repeatFam = classFam, classFam
            if not seqid in seqid2size:
                continue
            
            # Store repeat info if requested
            if repeatClass in repeatFamilies or repeatFam in repeatFamilies:
                if seqid not in seqid2repeats:
                    seqid2repeats[seqid] = []
                
                # Format and store a SeqFeature
                feature = SeqFeature.SeqFeature(
                    location=SeqFeature.FeatureLocation(int(start) - 1, int(end)), # 1-based to 0-based
                    type="repeat",
                    id=repeatID,
                    qualifiers={"class": repeatClass, "family": repeatFam}
                )
                
                seqid2repeats[seqid].append(feature)
    return seqid2repeats

def build_track_specs(args, seqid2size):
    '''
    Parses each annotation file and produces a list of track specifications in
    the order they should be plotted (outermost to innermost).
    
    Parameters:
        args -- the argparse Namespace object after validate_args() has been run
        seqid2size -- a dictionary of contig IDs to their lengths
    Returns:
        trackSpecs -- a list of dictionaries with keys "type" ("gene" or "repeat"),
                      "features" (a seqid2features dictionary), "repeatType"
                      (a string, or None for gene tracks), "label", and "colour"
    '''
    # Only disambiguate labels with file names if needed
    numGFF3 = args.annotationTypes.count("gff3")
    numRepeats = args.annotationTypes.count("repeatmasker")
    
    trackSpecs = []
    for annotationFile, annotationType in zip(args.annotationFiles, args.annotationTypes):
        baseName = os.path.basename(annotationFile)
        if annotationType == "gff3":
            trackSpecs.append({
                "type": "gene",
                "features": parse_gff3(annotationFile),
                "repeatType": None,
                "label": "Gene" if numGFF3 == 1 else f"Gene ({baseName})"
            })
        else:
            seqid2repeats = parse_repeatmasker(annotationFile, args.repeatFamilies, seqid2size)
            for repeatType in args.repeatFamilies:
                trackSpecs.append({
                    "type": "repeat",
                    "features": seqid2repeats,
                    "repeatType": repeatType,
                    "label": repeatType if numRepeats == 1 else f"{repeatType} ({baseName})"
                })
    
    for trackIndex, spec in enumerate(trackSpecs):
        spec["colour"] = colourPalette[trackIndex]
        if args.legendLabels != []:
            spec["label"] = args.legendLabels[trackIndex]
    return trackSpecs

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
    
    x = np.arange(0, contigSize, adjustedBinSize)
    y = np.zeros(len(x), dtype=int)
    for feature in features:
        startBin = int(feature.location.start) // adjustedBinSize
        endBin = min((int(feature.location.end) - 1) // adjustedBinSize, len(y) - 1) # -1 since end is exclusive
        for binIndex in range(startBin, endBin+1):
            y[binIndex] += 1
    return x, y

def select_track_features(spec, seqid):
    '''
    Parameters:
        spec -- a track specification dictionary as produced by build_track_specs()
        seqid -- a string indicating the contig to obtain features for
    Returns:
        features -- a list of SeqFeature objects to plot on this track
    '''
    features = spec["features"].get(seqid, [])
    if spec["type"] == "gene":
        return [ feature for feature in features if feature.type == "gene" ]
    else:
        repeatType = spec["repeatType"]
        return [ feature for feature in features
                if feature.qualifiers["class"] == repeatType or feature.qualifiers["family"] == repeatType ]

def bin_track_features(trackSpecs, seqid2size, binSize):
    '''
    Bins the features of each track for each contig, storing the result in each
    track specification under the "bins" key as a dictionary of seqid: (x, y).
    
    Parameters:
        trackSpecs -- a list of track specification dictionaries as produced by
                      build_track_specs()
        seqid2size -- a dictionary of contig IDs to their lengths
        binSize -- an integer value indicating the size of the bins/windows
    '''
    for spec in trackSpecs:
        spec["bins"] = {
            seqid: bin_values(select_track_features(spec, seqid), contigSize, binSize)
            for seqid, contigSize in seqid2size.items()
        }

def set_scale_limits(trackSpecs, scaleMode):
    '''
    Sets the "vmax" key of each binned track specification so that the density
    scale maximum is shared across contigs ("track") or across all tracks and
    contigs ("global"). With "sector" mode, nothing is set and each track is
    scaled for each contig individually during plotting. The scale minimum is
    always 0.
    
    Parameters:
        trackSpecs -- a list of track specification dictionaries that have been
                      through bin_track_features()
        scaleMode -- a string with value "sector", "track", or "global"
    '''
    if scaleMode == "sector":
        return
    
    def get_max(specs):
        return max(1, max( int(y.max()) for spec in specs for _, y in spec["bins"].values() )) # avoid a zero-height scale
    
    if scaleMode == "global":
        vmax = get_max(trackSpecs)
        for spec in trackSpecs:
            spec["vmax"] = vmax
    else:
        for spec in trackSpecs:
            spec["vmax"] = get_max([spec])

def plot_track(spec, sector, currentPosition, binSize,
               showOutline, TRACK_HEIGHT):
    '''
    Plots a single track for a sector either as individual features (binSize == -1)
    or as a density bar chart.
    
    Parameters:
        spec -- a track specification dictionary as produced by build_track_specs()
        sector -- a pyCirclize Sector object
        currentPosition -- a number indicating the outer radius of this track
        binSize -- an integer value indicating the size of the bins/windows, or -1
        showOutline -- a boolean indicating whether to outline the track
        TRACK_HEIGHT -- a number indicating the radial height of the track
    Returns:
        track -- the pyCirclize Track object
    '''
    track = sector.add_track((currentPosition-TRACK_HEIGHT, currentPosition))
    if showOutline:
        track.axis()
    
    if binSize == -1:
        features = select_track_features(spec, sector.name)
        if features != []:
            track.genomic_features(features, fc=spec["colour"])
    else:
        x, y = spec["bins"][sector.name]
        vmax = spec["vmax"] if "vmax" in spec else max(1, max(y))
        track.bar(x, y, vmax=vmax, width=binSize, fc=spec["colour"], align="edge")
    return track

def main():
    usage = """%(prog)s produces a circos plot representation of a genome and its annotations.
    Annotations are provided as one or more GFF3 and/or RepeatMasker .out files to -a,
    with each file being plotted as track(s) in the order given, from outermost to innermost.
    GFF3 files produce one track of gene features, while RepeatMasker files produce one
    track for each value given to --families.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Required arguments
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-a", dest="annotationFiles",
                   required=True,
                   nargs="+",
                   help="""Specify the location of one or more GFF3 and/or RepeatMasker
                   .out files; tracks are plotted in the order given, from outermost
                   to innermost""")
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
    p.add_argument("--scale", dest="scaleMode",
                   required=False,
                   choices=["sector", "track", "global"],
                   help="""Optionally, when --binSize is set, specify how the density
                   scale maximum is determined (the minimum is always 0); 'sector'
                   scales each track for each contig individually, 'track' shares
                   the maximum across all contigs for each track, and 'global'
                   shares it across all tracks and contigs; default is 'sector'""",
                   default="sector")
    p.add_argument("--families", dest="repeatFamilies",
                   required=False,
                   nargs="+",
                   help="""If you are plotting RepeatMasker output, specify
                   the repeat classes and/or families to plot; these apply to
                   all RepeatMasker files given to -a""",
                   default = [])
    p.add_argument("--labels", dest="legendLabels",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify the legend label for each track in
                   plotting order; the number of labels must match the number of
                   tracks plotted (one per GFF3, and one per --families value for
                   each RepeatMasker file)""",
                   default = [])
    p.add_argument("--outline", dest="sectorOutline",
                   required=False,
                   action="store_true",
                   help="Provide this flag to outline each sector in the plot",
                   default=False)
    p.add_argument("--totalHeight", dest="totalHeight",
                   required=False,
                   type=int,
                   help="""Optionally, specify the total height (as a value from 10
                   to 85) indicating how much space the circle tracks should occupy;
                   the remaining portion of 100 is empty space at the circle centre;
                   default is 50""",
                   default=50)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    seqid2size = { record.id:len(record) for record in genomeRecords }
    
    # Raise error if no contigs are found in genome FASTA
    if seqid2size == {}:
        raise ValueError(f"No contigs found in genome FASTA '{args.genomeFasta}'; is it actually a FASTA file?")
    
    # Parse annotation files into ordered track specifications
    trackSpecs = build_track_specs(args, seqid2size)
    
    # Bin features and determine density scale limits
    if args.binSize != -1:
        bin_track_features(trackSpecs, seqid2size, args.binSize)
        set_scale_limits(trackSpecs, args.scaleMode)
    
    # Initialise circos instance
    space = 0 if len(seqid2size) == 1 else 2
    if args.legendLocation == "top":
        circos = Circos(seqid2size, space=space,
                        start=10 if len(trackSpecs) < 5 else 20,
                        end=350 if len(trackSpecs) < 5 else 340,
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
    TRACK_HEIGHT = int(args.totalHeight / len(trackSpecs))
    
    for sector in circos.sectors:
        currentPosition = 95
        
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
        
        # Plot each annotation track in order
        for spec in trackSpecs:
            track = plot_track(spec, sector, currentPosition, args.binSize,
                               args.sectorOutline, TRACK_HEIGHT)
            currentPosition -= (TRACK_HEIGHT + TRACK_GAP)
            
            # Plot track legend
            if sector.name == circos.sectors[0].name:
                if args.legendLocation == "top":
                    circos.text(spec["label"], r=track.r_center, color=spec["colour"])
                elif args.legendLocation == "left":
                    circos.text(spec["label"], r=track.r_center, deg=-90, color=spec["colour"])
    
    # Generate the plot
    fig = circos.plotfig()
    if args.legendLocation == "centre":
        handles = [ Patch(color=spec["colour"], label=spec["label"]) for spec in trackSpecs ]
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
