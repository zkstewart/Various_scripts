#! python3
# phase_contigs.py
# Script to receive set of query contigs (e.g., from hifiasm *.p_ctg.fasta) and
# a reference assembly, and phase those contigs into partitions that can be
# separately assembled through RagTag.
# Note that a significant amount of the code development here was accomplished
# with Microsoft Copilot

import os, argparse, sys, shutil, re, subprocess, platform, gzip, codecs
import numpy as np

from intervaltree import IntervalTree
from Bio import SeqIO
from contextlib import contextmanager
from itertools import combinations, product

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import RANSACRegressor

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_SeqIO, ZS_AlignIO

# Define functions
def validate_args(args):
    def _not_found_error(program):
        raise FileNotFoundError(f"{program} not discoverable in your system PATH and was not specified as an argument.")
    def _specified_wrong_error(program, path):
        raise FileNotFoundError(f"{program} was not found at the indicated location '{path}'")
    
    # Validate input file locations
    for inputFile in args.inputGenome:
        if not os.path.isfile(inputFile):
            raise FileNotFoundError(f"Input genome file '{inputFile}' not found.")
    args.inputGenome = [ os.path.abspath(inputFile) for inputFile in args.inputGenome ]
    
    if not os.path.isfile(args.referenceGenome):
        raise FileNotFoundError(f"Reference genome file '{args.referenceGenome}' not found.")
    args.referenceGenome = os.path.abspath(args.referenceGenome)
    
    # Validate program discoverability
    if args.minimap2 == None:
        args.minimap2 = shutil.which("minimap2")
        if args.minimap2 == None:
            _not_found_error("minimap2")
    else:
        if not os.path.isfile(args.minimap2):
            _specified_wrong_error("minimap2", args.minimap2)
    
    if args.samtools == None:
        args.samtools = shutil.which("samtools")
        if args.samtools == None:
            _not_found_error("samtools")
    else:
        if not os.path.isfile(args.samtools):
            _specified_wrong_error("samtools", args.samtools)
    
    if args.ragtag == None:
        args.ragtag = shutil.which("ragtag.py")
        if args.ragtag == None:
            _not_found_error("ragtag.py")
    else:
        if not os.path.isfile(args.ragtag):
            _specified_wrong_error("ragtag.py", args.ragtag)
    
    # Validate output file location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"# Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("# But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"# Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        f.close()
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            f.close()
            return "utf-16"
        except UnicodeDecodeError:
            print(f"'{fileName}' is neither utf-8 nor utf-16 encoded; please convert to one of these formats.")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

def run_ragtag(inputFile, referenceFile, outputDir, ragtagPath, threads=1):
    '''
    Runs RagTag to scaffold the input file using the reference file as a guide.
    
    Parameters:
        inputFile -- a string indicating the location of the file to be scaffolded
        referenceFile -- a string indicating the location of the reference file
        outputDir -- a string indicating the location to write ragtag outputs
        ragtagPath -- a string indicating the location of the ragtag.py file
        threads -- (OPTIONAL) an integer indicating the number of threads when running ragtag
    '''
    # Format ragtag command
    cmd = [
        ragtagPath, "scaffold", "-t", str(threads), "-o", outputDir, "-r",
        referenceFile, inputFile
    ]
    
    # Run ragtag
    if platform.system() != "Windows":
        run_ragtag = subprocess.Popen(" ".join(cmd), shell = True,
                                      stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    else:
        run_ragtag = subprocess.Popen(cmd, shell = True,
                                      stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    ragtagout, ragtagerr = run_ragtag.communicate()
    if not "INFO: Finished running" in ragtagerr.decode("utf-8"):
        raise Exception('ragtag error text below\n' + ragtagerr.decode("utf-8"))

def find_longest_seq(fastaFile):
    '''
    Parse a FASTA file (without loading it into memory) to find the longest sequence.
    
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file
    Returns:
        longestSeqID -- a string indicating the ID of the longest sequence
        numSeqs -- an integer indicating the number of sequences in the FASTA file
    '''
    numSeqs = 0
    with open(fastaFile, "r") as fastaFile:
        records = SeqIO.parse(fastaFile, "fasta")
        longestSeq = [None, 0] # [seq, length]
        for record in records:
            numSeqs += 1
            length = len(record)
            if length > longestSeq[1]:
                longestSeq[0] = record.id
                longestSeq[1] = length
    return longestSeq[0], numSeqs

def parse_fasta_identifiers(fastaFile):
    '''
    Parse a FASTA file and return a set of all sequence identifiers
    
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file
    Returns:
        seqIDs -- a set containing strings for each sequence identifier found
                  within the FASTA file
    '''
    seqIDs = set()
    with open(fastaFile, "r") as fastaFile:
        records = SeqIO.parse(fastaFile, "fasta")
        for record in records:
            seqIDs.add(record.id)
    return seqIDs

def find_overlaps(intervals):
    if len(intervals) < 2:
        return intervals.copy()
    
    indexed = sorted(
        enumerate(intervals),
        key=lambda x: x[1][0]
    )
    
    overlapping = set()
    
    prev_idx, (prev_start, prev_end) = indexed[0]
    
    for curr_idx, (curr_start, curr_end) in indexed[1:]:
        if curr_start < prev_end:  # overlap
            overlapping.add(prev_idx)
            overlapping.add(curr_idx)
            prev_end = max(prev_end, curr_end)
        else:
            prev_idx = curr_idx
            prev_end = curr_end
    
    return sorted(overlapping)

def line_fit(qPoint, tPoint, alnLen):
    '''
    Runs outlier-aware line fitting to obtain the straight line fit
    through the query and alignment data points
    
    Parameters:
        qPoint -- the x-axis points relating to query alignments
        tPoint -- the y-axis points relating to target alignments
        alnLen -- the length of each aligned region
    '''
    # Identify outliers with RANSAC regression fitting
    ransac = RANSACRegressor(
        estimator=LinearRegression(),
        random_state=42
    )
    ransac.fit(qPoint.reshape(-1, 1), tPoint, sample_weight=alnLen)
    inliers = ransac.inlier_mask_
    
    # Filter data down to inliers only
    qFilt = qPoint[inliers]
    tFilt = tPoint[inliers]
    alnFilt = alnLen[inliers]
    
    # Run standard regression for model attributes
    model = LinearRegression()
    model.fit(qFilt.reshape(-1, 1), tFilt, sample_weight=alnFilt)
    
    slope = model.coef_[0]
    intercept = model.intercept_
    
    return slope, intercept

def overlap_bp(a_start, a_end, b_start, b_end):
    '''
    Calculate overlap between two intervals, assuming they are sorted.
    '''
    return max(
        0,
        min(a_end, b_end) - max(a_start, b_start)
    )

def interval_coverage(intervals):
    '''
    Calculate unique bp covered by a set of intervals.
    '''
    if len(intervals) == 0:
        return 0

    coords = sorted(
        [(s, e) for s, e, _ in intervals],
        key=lambda x: x[0]
    )

    merged = [list(coords[0])]

    for s, e in coords[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(
                merged[-1][1],
                e
            )
        else:
            merged.append([s, e])

    return sum(
        e - s
        for s, e in merged
    )

def total_overlap(intervals):
    '''
    Sum pairwise overlap within a haplotype.
    '''
    overlap = 0

    for (s1, e1, _), (s2, e2, _) in combinations(intervals, 2):
        overlap += overlap_bp(
            s1, e1,
            s2, e2
        )

    return overlap

def hap_score(intervals, overlap_penalty=3.0):
    '''
    Higher is better.
    '''
    coverage = interval_coverage(intervals)
    overlap = total_overlap(intervals)

    return (
        coverage
        - overlap_penalty * overlap
    )

def partition_score(hap1, hap2, overlap_penalty=3.0, balance_penalty=0.25):
    '''
    Global score for one partition. Higher is better.
    
    When both haplotypes give equal coverage of the reference scaffold,
    balance_cost approaches zero and no penalisation occurs.
    '''
    cov1 = interval_coverage(hap1)
    cov2 = interval_coverage(hap2)

    balance_cost = abs(
        cov1 - cov2
    )

    return (
        hap_score(hap1, overlap_penalty) + \
        hap_score(hap2, overlap_penalty) - \
        balance_penalty * balance_cost
    )

def optimise_haplotypes(intervals, overlap_penalty=3.0, balance_penalty=0.25):
    '''
    Exact diploid partition search.

    Parameters
    ----------
    intervals : list
        [(start, end, contigID), ...]

    Returns
    -------
    hap1
    hap2
    score
    '''

    n = len(intervals)

    best_score = float("-inf")
    best_hap1 = None
    best_hap2 = None
    best_assignment = None

    # Brute force each contig combination for optimal score
    for assignment in product(
        [0, 1],
        repeat=n
    ):

        # Skip products where all contigs are assigned to a single group
        if all(x == 0 for x in assignment):
            continue

        if all(x == 1 for x in assignment):
            continue
        
        # Separate contigs according to their potential grouping
        hap1 = [
            intervals[i]
            for i in range(n)
            if assignment[i] == 0
        ]

        hap2 = [
            intervals[i]
            for i in range(n)
            if assignment[i] == 1
        ]

        # Score the overall grouping
        score = partition_score(
            hap1,
            hap2,
            overlap_penalty,
            balance_penalty
        )

        # Hold onto the optimal (highest) score
        if score > best_score:
            best_score = score
            best_hap1 = hap1
            best_hap2 = hap2
            best_assignment = assignment

    return (
        best_hap1,
        best_hap2,
        best_assignment,
        best_score
    )

def print_partition_result(hap1, hap2):
    '''
    Provides a display of the phasing results for user verification.
    '''
    BUFFER = max([
        len(contig)
        for hap in (hap1, hap2)
        for _, _, contig in hap
    ]) + 2
    
    TOTAL = BUFFER + 1 + 11 + 11 # 11 since it is " :>10", +1 for unknown reason
    
    print("Phasing groups were defined as below:\n")
    print("Group 1")
    print("-" * TOTAL)
    for start, end, contig in sorted(hap1):
        print(
            "{contig:<{width}}".format(contig=contig, width=BUFFER),
            f" {int(start):>10}"
            f" {int(end):>10}"
        )
    
    print("\nGroup 2")
    print("-" * TOTAL)
    for start, end, contig in sorted(hap2):
        print(
            "{contig:<{width}}".format(contig=contig, width=BUFFER),
            f" {int(start):>10}"
            f" {int(end):>10}"
        )
    print()

def print_ragtag_result(ragtagAgpFile, groupNum):
    '''
    Provides a display akin to print_partition_result() to see
    how the phasing expectation compares to the actual scaffolding
    result.
    '''
    data = []
    with open(ragtagAgpFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            
            scaffID, scaffStart, scaffEnd, partNum, \
                componentType, componentID, componentStart, \
                componentEnd, strand = line.rstrip().split("\t")
            scaffID = scaffID.split("_RagTag")[0]
            
            # Skip gaps
            if strand == "align_genus":
                continue
            
            # Store inclusions
            data.append((scaffID, componentID, scaffStart, scaffEnd))
    
    BUFFER = max([
        len(componentID)
        for _, componentID, _, _ in data
    ]) + 2
    TOTAL = BUFFER + 1 + 11 + 11 # 11 since it is " :>10", +1 for unknown reason
    
    print(f"RagRag scaffolding for group {groupNum+1} was:")
    prevScaff = None
    for scaffID, componentID, scaffStart, scaffEnd in data:
        if scaffID != prevScaff:
            print(f"\n{scaffID}")
            print("-" * TOTAL)
            prevScaff = scaffID
        
        print(
            "{contig:<{width}}".format(contig=componentID, width=BUFFER),
            f" {int(scaffStart):>10}"
            f" {int(scaffEnd):>10}"
        )
    print()

## Main
def main():
    MULTILINE_LENGTH = 70
    
    # User input
    usage = """%(prog)s automates the scaffolding of a haplotype assembly using a reference assembly.
    It first ensures that haplotypes are equivalently assigned (e.g., hap1 is, to the best of the program's ability,
    consistently hap1 rather than hap2) and then uses ragtag to scaffold the haplotype assembly where
    necessary. Note 1: Contigs for hap1 and hap2 of the reference genome must be equivalently named.
    Note 2: Defaults for --minQueryAlign and --minAlignLen are based on paf2dotplot.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputGenome",
                   required=True,
                   nargs="+",
                   help="Input one or more genome FASTA files requiring scaffolding")
    p.add_argument("-r", dest="referenceGenome",
                   required=True,
                   help="Input FASTA file for use as a reference")
    p.add_argument("-c", dest="referenceContig",
                   required=True,
                   help="Specify the reference contig to be used")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (minimap2)
    p.add_argument("--preset", dest="preset",
                   required=False,
                   choices=["asm5", "asm10", "asm20"],
                   help="""Optionally, specify the preset to use for minimap2;
                   default == 'asm20'""",
                   default="asm20")
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify the number of threads to use for minimap2;
                   default == 1""",
                   default=1)
    # Opts (programs)
    p.add_argument("--minimap2", dest="minimap2",
                   required=False,
                   help="""Optionally, specify the minimap2 executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--samtools", dest="samtools",
                   required=False,
                   help="""Optionally, specify the samtools executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--ragtag", dest="ragtag",
                   required=False,
                   help="""Optionally, specify the ragtag.py file
                   if it is not discoverable in the path""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Set up the working directory structure
    referenceDir = os.path.join(args.outputDirectory, "reference")
    os.makedirs(referenceDir, exist_ok=True)
    
    inputDir = os.path.join(args.outputDirectory, "input")
    os.makedirs(inputDir, exist_ok=True)
    
    alignmentDir = os.path.join(args.outputDirectory, "alignment")
    os.makedirs(alignmentDir, exist_ok=True)
    
    # Extract the reference contig
    refContigLocation = os.path.join(referenceDir, "reference.fasta")
    if not os.path.exists(refContigLocation):
        with read_gz_file(args.referenceGenome) as fileIn:
            records = SeqIO.to_dict(SeqIO.parse(fileIn, "fasta"))
            contig = records[args.referenceContig]
            with open(refContigLocation, "w") as fileOut:
                fileOut.write(contig.format("fasta"))
    
    if not os.path.exists(f"{refContigLocation}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(refContigLocation, args.samtools)
    
    # Extract the input contigs
    inputContigsLocation = os.path.join(inputDir, "input.fasta")
    if not os.path.exists(inputContigsLocation):
        foundContigs = set()
        with open(inputContigsLocation, "w") as fileOut:
            for inputFile in args.inputGenome:
                with read_gz_file(inputFile) as fileIn:
                    records = SeqIO.parse(fileIn, "fasta")
                    for record in records:
                        contigID = record.id
                        if contigID in foundContigs:
                            raise ValueError(f"Found duplicated contig '{contigID}' in '{inputFile}'")
                        fileOut.write(record.format("fasta"))
                        foundContigs.add(contigID)
    
    if not os.path.exists(f"{inputContigsLocation}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(inputContigsLocation, args.samtools)
    
    # Run minimap2 of input files against the reference files
    minimapFileLocation = os.path.join(alignmentDir, "minimap2.paf")
    if not (os.path.exists(minimapFileLocation) and os.path.exists(minimapFileLocation + ".ok")):
        runner = ZS_AlignIO.Minimap2(inputContigsLocation, refContigLocation, args.preset, args.minimap2, args.threads)
        runner.minimap2(minimapFileLocation, force=True) # allow overwriting since the flag was not created
        open(minimapFileLocation + ".ok", "w").close()
    else:
        print(f"# Minimap2 alignment has already been performed; skipping.")
    
    # Parse minimap2 PAF files
    pafDict = {}
    qlenDict = {}
    with open(minimapFileLocation, "r") as fileIn:
        for line in fileIn:
            # Extract relevant data
            sl = line.rstrip("\r\n ").split("\t")
            qid, qlen, qstart, qend, strand, tid, tlen, tstart, tend, \
                numresidues, lenalign, mapq = sl[0:12]
            
            # Skip if the alignment doesn't meet Q score minimum
            if int(mapq) != 60:
                continue
            
            # Store the alignment
            pafDict.setdefault(qid, {"q": [], "t": []})
            pafDict[qid]["q"].append((int(qstart), int(qend)))
            pafDict[qid]["t"].append((int(tstart), int(tend)))
            
            qlenDict[qid] = int(qlen)
    tlen = int(tlen)
    
    # Retain only segments of the query that do not align redundantly
    for key in pafDict.keys():
        overlappingIndices = find_overlaps(pafDict[key]["q"])
        pafDict[key]["q"] = [
            interval
            for i, interval in enumerate(pafDict[key]["q"])
            if i not in overlappingIndices
        ]
        pafDict[key]["t"] = [
            interval
            for i, interval in enumerate(pafDict[key]["t"])
            if i not in overlappingIndices
        ]
    
    # Convert alignment intervals into scatterplot data for regression
    intervals = []
    for key in pafDict.keys():
        
        # Run line through data points
        qPoint = np.array([ (start + end) / 2 for start, end in pafDict[key]["q"] ])
        tPoint = np.array([ (start + end) / 2 for start, end in pafDict[key]["t"] ])
        alnLen = np.array([ end - start for start, end in pafDict[key]["q"] ])
        
        slope, intercept = line_fit(qPoint, tPoint, alnLen)
        
        # Derive the projected start and end points
        if slope > 0:
            start = intercept if intercept > 0 else 1 # if projection would be longer than the reference
            end = (slope * qlenDict[key]) + intercept
        else:
            start = (slope * qlenDict[key]) + intercept
            end = intercept if intercept < tlen else tlen # if projection would be longer than the reference
        
        # Store this contig's aligned region
        intervals.append((start, end, key))
    
    # Separate contigs into two haplotype groupings
    hap1, hap2, assignment, score = (
        optimise_haplotypes(
            intervals,
            overlap_penalty=3.0,
            balance_penalty=0.25
        )
    )
    
    # Let the user know what the result was
    print_partition_result(hap1, hap2)
    
    # Load input sequences prior to scaffolding
    with open(inputContigsLocation, "r") as fileIn:
        records = SeqIO.to_dict(SeqIO.parse(fileIn, "fasta"))
    
    # Generate scaffolds for each haplotype chromosome
    for i, hap in enumerate((hap1, hap2)):
        contigs = [ contig for _, _, contig in hap ]
        
        runDir = os.path.join(args.outputDirectory, f"group{i+1}")
        os.makedirs(runDir, exist_ok=True)
        
        if not os.path.isfile(runDir + ".ok"):
            # Produce the input FASTA
            groupFile = os.path.join(runDir, "group.fasta")
            with open(groupFile, "w") as fileOut:
                for contig in contigs:
                    fileOut.write(records[contig].format("fasta"))
            
            # Run RagTag to scaffold the sequences
            run_ragtag(groupFile, refContigLocation,
                       os.path.join(runDir),
                       args.ragtag, threads=args.threads)
            
            # Create a flag to indicate that RagTag was successful
            open(runDir + ".ok", "w").close()
            
        # Print RagTag scaffolding result for comparison to the phasing expectation
        ragtagAgpFile = os.path.join(runDir, f"ragtag.scaffold.agp")
        print_ragtag_result(ragtagAgpFile, i)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
