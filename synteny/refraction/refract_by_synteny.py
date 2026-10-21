#! python3
'''
refract_by_synteny.py

Projects ("refracts") the coordinates of a GFF3 annotation onto the
coordinate space of a different, broadly collinear genome using the synteny
predictions made by SyRI. The intention is to allow gene annotations from
multiple genomes to be plotted on a single shared karyotype (e.g., a circos
density track) such that syntenically equivalent regions approximately line up.

How positions are mapped:
    1) The SYNAL and INVAL alignments from SyRI form a collinear "backbone"
       of anchor segments between the input genome and the target genome.
       Positions inside an anchor are mapped linearly within that segment;
       positions inside an INVAL segment are mapped in reverse, and features
       anchored within an inversion have their strand flipped.
    2) Positions between anchors (unaligned sequence, HDRs, translocated or
       duplicated regions, indels) are linearly interpolated between the
       flanking anchors, stretching or compressing them to fit the gap in
       the target. Chromosome ends act as implicit anchors so that every
       feature lands within the target chromosome's length.
    3) Optionally, genes lying within translocated regions (TRANS / INVTR)
       can be moved to the position of their translocated counterpart
       instead of being interpolated.

All features descending from the same top-level feature (e.g., a gene and its
mRNA, exon and CDS children) share one mapping decision, so a gene model is
never split across mapping modes. Each output feature receives a
"refract_mode" attribute describing how it was placed (syn, inv, interp,
trans, invtr) to help identify interpolation artefacts in downstream plots.

Note that stretching changes feature lengths, so output models are only
suitable for positional analyses (e.g., density) and not for sequence
extraction; CDS phases are left untouched.
'''

import os, sys, math, argparse
from bisect import bisect_right

BACKBONE_TYPES = {"SYNAL": "syn", "INVAL": "inv"}
TRANSLOCATION_TYPES = {"TRANS": "trans", "INVTR": "invtr"}

# Define functions
def validate_args(args):
    for fileName in [args.gff3, args.targetLengths, args.syri]:
        if not os.path.isfile(fileName):
            raise FileNotFoundError(f"Input file '{fileName}' does not exist!")
    if args.queryLengths is not None and not os.path.isfile(args.queryLengths):
        raise FileNotFoundError(f"Input file '{args.queryLengths}' does not exist!")
    if args.outputFileName is not None:
        if os.path.isfile(args.outputFileName):
            raise FileExistsError(f"Output file '{args.outputFileName}' already exists!")
        outputDir = os.path.dirname(os.path.abspath(args.outputFileName))
        if not os.path.isdir(outputDir):
            raise FileNotFoundError(f"Output directory '{outputDir}' does not exist!")
    if args.minTransLength < 0:
        raise ValueError("--min_trans_length must be >= 0")
    if args.gapsBed is not None:
        if os.path.isfile(args.gapsBed):
            raise FileExistsError(f"Gaps BED file '{args.gapsBed}' already exists!")
        gapsDir = os.path.dirname(os.path.abspath(args.gapsBed))
        if not os.path.isdir(gapsDir):
            raise FileNotFoundError(f"Gaps BED directory '{gapsDir}' does not exist!")
    if args.minGapLength < 1:
        raise ValueError("--min_gap_length must be >= 1")

def parse_lengths(fileName):
    '''
    Parses a FASTA file or a two-column TSV (chrom<tab>length; a samtools .fai
    also works) into a dictionary of sequence lengths. The format is detected
    by whether the first non-blank line begins with '>'.
    '''
    lengths = {}
    with open(fileName, "r") as fileIn:
        firstLine = ""
        for firstLine in fileIn:
            if firstLine.strip() != "":
                break
        fileIn.seek(0)

        if firstLine.startswith(">"):
            seqID = None
            for line in fileIn:
                if line.startswith(">"):
                    seqID = line[1:].split()[0]
                    if seqID in lengths:
                        raise ValueError(f"Duplicate sequence ID '{seqID}' in '{fileName}'")
                    lengths[seqID] = 0
                elif seqID is not None:
                    lengths[seqID] += len(line.strip())
        else:
            for line in fileIn:
                if line.strip() == "" or line.startswith("#"):
                    continue
                sl = line.rstrip("\r\n").split("\t")
                if len(sl) < 2:
                    raise ValueError(f"Lengths file '{fileName}' has a line without two " +
                                     f"tab-separated columns: '{line.rstrip()}'")
                lengths[sl[0]] = int(sl[1])
    return lengths

class Segment:
    '''
    A pairing of a query (input genome) interval with a target interval.
    Query coordinates are always ascending; inverted segments map the query
    start to the target end.
    '''
    def __init__(self, qchr, qstart, qend, tchr, tstart, tend, mode, inverted):
        self.qchr, self.tchr = qchr, tchr
        self.qstart, self.qend = qstart, qend
        self.tstart, self.tend = tstart, tend
        self.mode = mode
        self.inverted = inverted

    def __len__(self):
        return self.qend - self.qstart + 1

    def map(self, pos):
        if self.qend == self.qstart:
            frac = 0.5
        else:
            frac = (pos - self.qstart) / (self.qend - self.qstart)
        frac = min(max(frac, 0.0), 1.0)
        if self.inverted:
            return round(self.tend - frac * (self.tend - self.tstart))
        return round(self.tstart + frac * (self.tend - self.tstart))

def parse_syri(fileName, inputSide):
    '''
    Streams a SyRI .out file, keeping only the lines relevant to synteny
    anchoring. SNPs, indels, NOTAL etc. are only used to record the maximum
    coordinate seen on each side, which is a fallback for chromosome lengths.

    Returns:
        backbone -- dict of {qchr: [Segment, ...]} from SYNAL/INVAL lines
        translocations -- dict of {qchr: [Segment, ...]} from TRANS/INVTR lines
        queryMax -- dict of {qchr: maxCoordinate} for the input genome
        targetMax -- dict of {tchr: maxCoordinate} for the target genome
    '''
    backbone, translocations = {}, {}
    queryMax, targetMax = {}, {}

    with open(fileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if len(sl) < 11:
                continue
            lchr, lstart, lend, rchr, rstart, rend, alignType = \
                sl[0], sl[1], sl[2], sl[5], sl[6], sl[7], sl[10]

            if inputSide == "left":
                qchr, qa, qb, tchr, ta, tb = lchr, lstart, lend, rchr, rstart, rend
            else:
                qchr, qa, qb, tchr, ta, tb = rchr, rstart, rend, lchr, lstart, lend

            # Track the extent of each chromosome on both sides
            if qchr != "-" and qa != "-":
                queryMax[qchr] = max(queryMax.get(qchr, 0), int(qa), int(qb))
            if tchr != "-" and ta != "-":
                targetMax[tchr] = max(targetMax.get(tchr, 0), int(ta), int(tb))

            if alignType in BACKBONE_TYPES:
                storage, mode = backbone, BACKBONE_TYPES[alignType]
            elif alignType in TRANSLOCATION_TYPES:
                storage, mode = translocations, TRANSLOCATION_TYPES[alignType]
            else:
                continue

            # SyRI reports inversions by giving the right side's coordinates in
            # descending order (INVAL) or by the annotation type itself (INVTR)
            qa, qb, ta, tb = int(qa), int(qb), int(ta), int(tb)
            inverted = (qa > qb) != (ta > tb) or mode in ("inv", "invtr")
            storage.setdefault(qchr, []).append(Segment(
                qchr, min(qa, qb), max(qa, qb),
                tchr, min(ta, tb), max(ta, tb),
                mode, inverted
            ))
    return backbone, translocations, queryMax, targetMax

def collinear_chain(segments):
    '''
    Selects the subset of segments that is collinear (ascending in both query
    and target coordinates) with the greatest total query length. SyRI's
    SYN/INV annotations should already be collinear, so this is a safeguard
    that typically discards nothing (or only tiny overlapping alignments).
    '''
    segments = sorted(segments, key=lambda s: (s.qstart, s.qend))
    if len(segments) == 0:
        return segments

    # Weighted longest increasing subsequence over target start
    score = [len(s) for s in segments]
    prev = [-1] * len(segments)
    for i, si in enumerate(segments):
        for j in range(i):
            sj = segments[j]
            if sj.tstart <= si.tstart and sj.qend <= si.qend and score[j] + len(si) > score[i]:
                score[i] = score[j] + len(si)
                prev[i] = j

    i = max(range(len(segments)), key=lambda k: score[k])
    chain = []
    while i != -1:
        chain.append(segments[i])
        i = prev[i]
    return chain[::-1]

class ChromosomeMapper:
    '''
    Maps positions on one input chromosome onto its target chromosome using
    a collinear backbone of anchor segments, interpolating linearly between them.
    '''
    def __init__(self, qchr, tchr, segments, queryLength, targetLength):
        self.qchr, self.tchr = qchr, tchr
        self.segments = collinear_chain(segments)
        self.queryLength, self.targetLength = queryLength, targetLength
        self.qstarts = [s.qstart for s in self.segments]

    def find_segment(self, pos):
        '''
        Returns (index of the last segment starting at or before pos, whether
        pos lies within that segment).
        '''
        i = bisect_right(self.qstarts, pos) - 1
        if i >= 0 and pos <= self.segments[i].qend:
            return i, True
        return i, False

    def locate(self, pos):
        '''
        Returns the Segment containing pos, or None if pos lies between anchors.
        '''
        i, inside = self.find_segment(pos)
        return self.segments[i] if inside else None

    def map(self, pos):
        i, inside = self.find_segment(pos)
        if inside:
            return self.segments[i].map(pos)

        # Interpolate between flanking anchors, with chromosome ends as fallback
        if i >= 0:
            leftQ, leftT = self.segments[i].qend, self.segments[i].tend
        else:
            leftQ, leftT = 0, 0
        if i + 1 < len(self.segments):
            rightQ, rightT = self.segments[i+1].qstart, self.segments[i+1].tstart
        else:
            rightQ, rightT = self.queryLength + 1, self.targetLength + 1
        rightT = max(rightT, leftT) # guards against tiny overlaps in the target

        if rightQ == leftQ:
            return leftT
        mapped = leftT + (pos - leftQ) * (rightT - leftT) / (rightQ - leftQ)
        return min(max(round(mapped), 1), self.targetLength)

    def gaps(self):
        '''
        Yields each interpolated gap between anchors (including the chromosome
        ends) as (qstart, qend, tstart, tend) in 1-based inclusive coordinates.
        A gap that is compressed to nothing in the target has tend == tstart - 1.
        '''
        leftQ, leftT = 0, 0
        for s in self.segments + [None]:
            if s is not None:
                rightQ, rightT = s.qstart, s.tstart
            else:
                rightQ, rightT = self.queryLength + 1, self.targetLength + 1
            rightT = max(rightT, leftT + 1)
            if rightQ - leftQ > 1:
                yield leftQ + 1, rightQ - 1, leftT + 1, rightT - 1
            if s is not None:
                leftQ, leftT = max(leftQ, s.qend), max(leftT, s.tend)

def build_mappers(backbone, targetLengths, queryLengths):
    '''
    Creates a ChromosomeMapper for each input chromosome. The target chromosome
    is the one receiving the greatest share of that chromosome's backbone
    alignment; backbone segments to any other chromosome are discarded.
    '''
    mappers = {}
    for qchr, segments in backbone.items():
        tchrLengths = {}
        for s in segments:
            tchrLengths[s.tchr] = tchrLengths.get(s.tchr, 0) + len(s)
        tchr = max(tchrLengths, key=tchrLengths.get)
        if tchr not in targetLengths:
            print(f"WARNING: target chromosome '{tchr}' (syntenic to input '{qchr}') is not " +
                  "in the target lengths file; its features will be treated as unmapped",
                  file=sys.stderr)
            continue
        mappers[qchr] = ChromosomeMapper(
            qchr, tchr,
            [s for s in segments if s.tchr == tchr],
            queryLengths[qchr], targetLengths[tchr]
        )
    return mappers

def find_translocation(translocations, qchr, pos, minLength, allowInterchromosomal, targetLengths):
    '''
    Returns the longest qualifying translocation Segment containing pos, or None.
    '''
    best = None
    for s in translocations.get(qchr, []):
        if not (s.qstart <= pos <= s.qend):
            continue
        if len(s) < minLength or s.tchr not in targetLengths:
            continue
        if s.tchr != qchr and not allowInterchromosomal:
            continue
        if best is None or len(s) > len(best):
            best = s
    return best

def parse_gff3(fileName):
    '''
    Reads a GFF3 into a list of lines; feature lines are split into columns,
    and comment lines are retained as strings. Also returns a dict linking
    each feature ID to its (first) Parent ID.
    '''
    lines, parents = [], {}
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            line = line.rstrip("\r\n")
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or line.strip() == "":
                lines.append(line)
                continue
            sl = line.split("\t")
            if len(sl) != 9:
                raise ValueError(f"GFF3 line does not have 9 columns: '{line}'")
            attributes = parse_attributes(sl[8])
            if "ID" in attributes and "Parent" in attributes:
                parents[attributes["ID"]] = attributes["Parent"].split(",")[0]
            lines.append(sl)
    return lines, parents

def parse_attributes(attributeString):
    attributes = {}
    for field in attributeString.strip().strip(";").split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attributes[key.strip()] = value
    return attributes

def find_root(featureID, parents):
    seen = set()
    while featureID in parents and featureID not in seen:
        seen.add(featureID)
        featureID = parents[featureID]
    return featureID

def refract(gff3Lines, parents, mappers, translocations, targetLengths,
            translocationMode, minTransLength, keepUnmapped):
    '''
    Decides how each top-level feature should be mapped, then applies that
    decision to it and all its descendants.

    Returns:
        outputLines -- list of strings to write
        modeCounts -- dict of {mode: count of top-level features}
    '''
    # Resolve each feature line to its root (top-level feature) ID
    rootOf = []
    for i, line in enumerate(gff3Lines):
        if isinstance(line, str):
            rootOf.append(None)
            continue
        attributes = parse_attributes(line[8])
        if "ID" in attributes:
            rootOf.append(find_root(attributes["ID"], parents))
        elif "Parent" in attributes:
            rootOf.append(find_root(attributes["Parent"].split(",")[0], parents))
        else:
            rootOf.append(f"__line{i}")

    # Determine the span of each root (uses all descendants in case the root is absent)
    rootSpans = {}
    for line, root in zip(gff3Lines, rootOf):
        if root is None:
            continue
        chrom, start, end = line[0], int(line[3]), int(line[4])
        if root not in rootSpans:
            rootSpans[root] = [chrom, start, end]
        else:
            rootSpans[root][1] = min(rootSpans[root][1], start)
            rootSpans[root][2] = max(rootSpans[root][2], end)

    # Decide the mapping for each root: (mode, tchr, positionFunction, inverted)
    decisions = {}
    modeCounts = {}
    for root, (chrom, start, end) in rootSpans.items():
        midpoint = (start + end) // 2
        mapper = mappers.get(chrom)
        decision = None

        if mapper is not None:
            segment = mapper.locate(midpoint)
            if segment is not None:
                decision = (segment.mode, mapper.tchr, mapper.map, segment.inverted)
        if decision is None and translocationMode != "ignore":
            segment = find_translocation(translocations, chrom, midpoint, minTransLength,
                                         translocationMode == "all", targetLengths)
            if segment is not None:
                decision = (segment.mode, segment.tchr, segment.map, segment.inverted)
        if decision is None and mapper is not None:
            decision = ("interp", mapper.tchr, mapper.map, False)
        if decision is None:
            decision = ("unmapped", chrom, None, False)

        decisions[root] = decision
        modeCounts[decision[0]] = modeCounts.get(decision[0], 0) + 1

    # Apply decisions
    outputLines = []
    for line, root in zip(gff3Lines, rootOf):
        if root is None:
            if not line.startswith("##sequence-region"):
                outputLines.append(line)
            continue
        mode, tchr, mapFunction, inverted = decisions[root]
        if mode == "unmapped":
            if keepUnmapped:
                outputLines.append("\t".join(line[:8] + [add_attribute(line[8], "refract_mode", mode)]))
            continue

        newStart, newEnd = sorted([mapFunction(int(line[3])), mapFunction(int(line[4]))])
        strand = line[6]
        if inverted:
            strand = {"+": "-", "-": "+"}.get(strand, strand)
        outputLines.append("\t".join([
            tchr, line[1], line[2], str(newStart), str(newEnd),
            line[5], strand, line[7], add_attribute(line[8], "refract_mode", mode)
        ]))
    return outputLines, modeCounts

def write_gaps_bed(mappers, outputFileName, minGapLength):
    '''
    Writes the interpolated gaps in target coordinates as a BED file with
    two extra columns: the stretch ratio (target length / input length; <1 is
    compression, >1 is stretching) and its log2 (useful for colour scales).
    The name column gives the gap's original location in the input genome.
    '''
    with open(outputFileName, "w") as fileOut:
        fileOut.write("#chrom\tstart\tend\tinput_region\tstretch_ratio\tlog2_stretch_ratio\n")
        for qchr, mapper in sorted(mappers.items()):
            for qstart, qend, tstart, tend in mapper.gaps():
                queryGapLength = qend - qstart + 1
                if queryGapLength < minGapLength:
                    continue
                targetGapLength = tend - tstart + 1
                ratio = targetGapLength / queryGapLength
                log2Ratio = f"{math.log2(ratio):.3f}" if ratio > 0 else "-inf"
                fileOut.write(f"{mapper.tchr}\t{tstart - 1}\t{tend}\t{qchr}:{qstart}-{qend}\t" +
                              f"{ratio:.4f}\t{log2Ratio}\n")

def add_attribute(attributeString, key, value):
    attributeString = attributeString.rstrip().rstrip(";")
    return f"{attributeString};{key}={value}" if attributeString else f"{key}={value}"

def main():
    usage = """%(prog)s refracts (stretches/compresses) the coordinates of a GFF3
    annotation onto the coordinate space of a target genome, using SyRI synteny
    predictions as anchors and linearly interpolating positions between anchors.

    The input GFF3's genome must be one side of the SyRI comparison; specify which
    side with --input_side. SyRI's .out files list the reference (-r) genome on the
    left and the query (-q) genome on the right.

    Features in inverted syntenic regions (INV) are mapped in reverse with their
    strand flipped. Features in translocated regions are interpolated by default,
    or can be moved to their translocated position with --translocations.
    """
    p = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
    # Required
    p.add_argument("-g", dest="gff3", required=True,
                   help="GFF3 file to refract")
    p.add_argument("-t", dest="targetLengths", required=True,
                   help="Target genome lengths as a FASTA or two-column TSV (chrom<tab>length)")
    p.add_argument("-s", dest="syri", required=True,
                   help="SyRI .out file comparing the GFF3's genome and the target genome")
    p.add_argument("--input_side", dest="inputSide", required=True, choices=["left", "right"],
                   help="Which side of the SyRI file the GFF3's genome is on")
    # Optional
    p.add_argument("-o", dest="outputFileName", default=None,
                   help="Output GFF3 file name (default: stdout)")
    p.add_argument("-q", dest="queryLengths", default=None,
                   help="""Optionally, the GFF3 genome's lengths as a FASTA or two-column TSV;
                   if not given, these are inferred from the SyRI file""")
    p.add_argument("--translocations", dest="translocationMode", default="ignore",
                   choices=["ignore", "intra", "all"],
                   help="""How to handle features in TRANS/INVTR regions: 'ignore' interpolates
                   them like any other unanchored position, 'intra' moves them to their
                   translocated position if on the same chromosome, 'all' also allows moves
                   between chromosomes (default: ignore)""")
    p.add_argument("--min_trans_length", dest="minTransLength", type=int, default=10000,
                   help="""Minimum length (bp, input genome) of a translocation for it to be
                   used with --translocations intra/all (default: 10000)""")
    p.add_argument("--keep_unmapped", dest="keepUnmapped", action="store_true", default=False,
                   help="""Keep features on chromosomes with no synteny to the target, at their
                   original coordinates (default: drop them)""")
    p.add_argument("--gaps_bed", dest="gapsBed", default=None,
                   help="""Optionally, write a BED file of the interpolated gaps (in target
                   coordinates) with their stretch/compression ratios""")
    p.add_argument("--min_gap_length", dest="minGapLength", type=int, default=10000,
                   help="""Minimum length (bp, input genome) of a gap for it to be written
                   to --gaps_bed (default: 10000)""")
    args = p.parse_args()
    validate_args(args)

    # Parse inputs
    targetLengths = parse_lengths(args.targetLengths)
    backbone, translocations, queryMax, targetMax = parse_syri(args.syri, args.inputSide)
    gff3Lines, parents = parse_gff3(args.gff3)

    # Check the target lengths are consistent with SyRI
    for tchr, maxCoord in targetMax.items():
        if tchr in targetLengths and maxCoord > targetLengths[tchr]:
            print(f"WARNING: SyRI reports coordinates up to {maxCoord} on target '{tchr}' which " +
                  f"exceeds its given length of {targetLengths[tchr]}; check --input_side and " +
                  "the target lengths file", file=sys.stderr)

    # Establish input genome lengths
    if args.queryLengths is not None:
        queryLengths = parse_lengths(args.queryLengths)
    else:
        queryLengths = dict(queryMax)
    for line in gff3Lines:
        if not isinstance(line, str) and line[0] in queryLengths:
            queryLengths[line[0]] = max(queryLengths[line[0]], int(line[4]))

    # Build mappers and refract the annotation
    mappers = build_mappers(backbone, targetLengths, queryLengths)
    outputLines, modeCounts = refract(gff3Lines, parents, mappers, translocations, targetLengths,
                                      args.translocationMode, args.minTransLength, args.keepUnmapped)

    # Write output with sequence-region pragmas for the target genome
    fileOut = open(args.outputFileName, "w") if args.outputFileName else sys.stdout
    try:
        headerDone = False
        for line in outputLines:
            if not headerDone and not line.startswith("##gff-version"):
                if line is outputLines[0]:
                    fileOut.write("##gff-version 3\n")
                for tchr in sorted({m.tchr for m in mappers.values()}):
                    fileOut.write(f"##sequence-region {tchr} 1 {targetLengths[tchr]}\n")
                headerDone = True
            fileOut.write(line + "\n")
    finally:
        if fileOut is not sys.stdout:
            fileOut.close()
    if args.gapsBed is not None:
        write_gaps_bed(mappers, args.gapsBed, args.minGapLength)

    # Report a summary
    print("Top-level features by refract_mode:", file=sys.stderr)
    for mode in ["syn", "inv", "interp", "trans", "invtr", "unmapped"]:
        if mode in modeCounts:
            print(f"\t{mode}\t{modeCounts[mode]}", file=sys.stderr)
    for qchr, mapper in sorted(mappers.items()):
        print(f"\t{qchr} -> {mapper.tchr}: {len(mapper.segments)} anchors " +
              f"(input length {mapper.queryLength}, target length {mapper.targetLength})",
              file=sys.stderr)

    print("Program completed successfully!", file=sys.stderr)

if __name__ == "__main__":
    main()
