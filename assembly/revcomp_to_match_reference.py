#! python3
# revcomp_to_match_reference.py
# Script to take in a directory of FASTA files and reverse-complement them to match a reference genome.

import os, argparse, sys, shutil, re, subprocess, platform
from intervaltree import IntervalTree
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_SeqIO, ZS_AlignIO

# Define functions
def validate_args(args):
    def _not_found_error(program):
        raise FileNotFoundError(f"{program} not discoverable in your system PATH and was not specified as an argument.")
    def _specified_wrong_error(program, path):
        raise FileNotFoundError(f"{program} was not found at the indicated location '{path}'")
    
    # Locate and validate input files
    args.fastaFiles = []
    fastaPrefixes = set()
    for location in args.inputLocations:
        location = os.path.abspath(location)
        # Handle files
        if os.path.isfile(location):
            args.fastaFiles.append(location)
        # Handle directories
        elif os.path.isdir(location):
            foundAny = False
            for f in os.listdir(location):
                if f.endswith(args.fastaSuffix):
                    fastaPrefix = os.path.basename(location).rsplit(".", maxsplit=1)[0]
                    if fastaPrefix in fastaPrefixes:
                        raise ValueError(f"Duplicate FASTA prefix found: '{fastaPrefix}'")
                    
                    args.fastaFiles.append(os.path.join(location, f))
                    foundAny = True
            if not foundAny:
                raise FileNotFoundError(f"No FASTA files found in directory '{location}' ending with '{args.fastaSuffix}'")
        # Handle other cases
        else:
            raise FileNotFoundError(f"Input FASTA file or directory '{location}' not found!")
    
    # Validate reference genome file
    if not os.path.isfile(args.referenceGenome):
        raise FileNotFoundError(f"Reference genome file '{args.referenceGenome}' not found.")
    args.referenceGenome = os.path.abspath(args.referenceGenome)
    
    # Validate numeric arguments
    if args.minAlignLen < 0:
        raise ValueError(f"--minAlignLen must be greater than or equal to 0.")
    if args.threads < 1:
        raise ValueError(f"--threads must be greater than or equal to 1.")
    if args.minQ < 0:
        raise ValueError(f"--minQ must be greater than or equal to 0.")
    
    # Validate program discoverability
    if args.minimap2 == None:
        args.minimap2 = shutil.which("minimap2")
        if args.minimap2 == None:
            _not_found_error("minimap2")
    else:
        if not os.path.isfile(args.minimap2):
            _specified_wrong_error("minimap2", args.minimap2)
    
    # Validate output file location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

## Main
def main():
    MULTILINE_LENGTH = 70
    
    # User input
    usage = """%(prog)s automates the reverse-complementing of FASTA files to match a reference genome.
    It will take in a directory of FASTA files and reverse-complement them to match the reference genome.
    All FASTA files must contain the same contigs as the reference genome.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputLocations",
                   required=True,
                   nargs="+",
                   help="Input one or more file names and/or directories containing FASTA files")
    p.add_argument("-r", dest="referenceGenome",
                   required=True,
                   help="Input genome FASTA file for use as a reference")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (behavioural)
    p.add_argument("--minAlignLen", dest="minAlignLen",
                   required=False,
                   type=int,
                   help="""Optionally, specify the minimum alignment length to
                   consider a match; default == 10000""",
                   default=10000)
    p.add_argument("--minQ", dest="minQ",
                   required=False,
                   type=int,
                   help="""Optionally, specify the minimap map Q score for considering an alignment
                   to be a match; default == 1""",
                   default=1)
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
    p.add_argument("--fastaSuffix", dest="fastaSuffix",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more suffixes to use when scanning for FASTA
                   files in input directories; default == '.fasta'""",
                   default=[".fasta"])
    p.add_argument("--minimap2", dest="minimap2",
                   required=False,
                   help="""Optionally, specify the minimap2 executable file
                   if it is not discoverable in the path""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Validate that all files have the same sequence IDs as the reference genome
    refSeqIDs = set()
    with open(args.referenceGenome, "r") as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                # Extract the sequence ID from the header line
                seqID = line.split()[0][1:]
                refSeqIDs.add(seqID)
    
    for fastaFile in args.fastaFiles:
        # Extract the sequence IDs from the input FASTA file
        inputSeqIDs = set()
        with open(fastaFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith(">"):
                    # Extract the sequence ID from the header line
                    seqID = line.split()[0][1:]
                    inputSeqIDs.add(seqID)
        # Check if the input sequence IDs match the reference sequence IDs
        if inputSeqIDs != refSeqIDs:
            raise ValueError(f"Input FASTA file '{fastaFile}' does not have the same sequence IDs " + 
                             f"as the reference genome '{args.referenceGenome}'")
    
    # Set up the working directory structure
    alignmentsDir = os.path.join(args.outputDirectory, "alignments")
    os.makedirs(alignmentsDir, exist_ok=True)
    
    # Run minimap2 of input files against the reference files
    for fastaFile in args.fastaFiles:
        outputName = os.path.join(alignmentsDir, f"{os.path.basename(fastaFile)}.paf")
        flagName = os.path.join(alignmentsDir, f"{os.path.basename(fastaFile)}.is.ok.flag")
        if not os.path.exists(flagName) or not os.path.exists(outputName):
            runner = ZS_AlignIO.Minimap2(fastaFile, args.referenceGenome, args.preset, args.minimap2, args.threads)
            runner.minimap2(outputName, force=True) # allow overwriting since the flag was not created
            open(flagName, "w").close()
    
    # Parse minimap2 PAF files and determine if reverse-complementing is needed
    needsRevComp = {}
    for fastaFile in args.fastaFiles:
        fastaPrefix = os.path.basename(fastaFile).rsplit(".", maxsplit=1)[0]
        needsRevComp[fastaPrefix] = set()
        
        # Parse PAF file and tally alignment length for each contig by strand
        pafFile = os.path.join(alignmentsDir, f"{os.path.basename(fastaFile)}.paf")
        alignLen = { r: {"+": 0, "-": 0} for r in refSeqIDs }
        with open(pafFile, "r") as fileIn:
            for line in fileIn:
                # Extract relevant data
                sl = line.rstrip("\r\n ").split("\t")
                qid, qlen, qstart, qend, strand, tid, tlen, tstart, tend, \
                    numresidues, lenalign, mapq = sl[0:12]
                
                # Skip if the alignment is not between the same contig
                if qid != tid:
                    continue
                
                # Convert to integers
                qlen, qstart, qend, tlen, tstart, tend, numresidues, lenalign, mapq = \
                    map(int, [qlen, qstart, qend, tlen, tstart, tend, numresidues, lenalign, mapq])
                refContigs.add(tid)
                
                # Skip if the alignment doesn't meet length minimum
                if lenalign < args.minAlignLen:
                    continue
                
                # Skip if the alignment doesn't meet Q score minimum
                if mapq < args.minQ:
                    continue
                
                # Store the alignment length
                alignLen[qid][strand] += lenalign
        
        # Check if the alignment length is greater for the reverse strand for each contig
        for contig in refSeqIDs:
            if alignLen[contig]["+"] == 0 and alignLen[contig]["-"] == 0:
                print(f"Warning: No alignments found for contig '{contig}' in file '{fastaFile}'.")
                continue
            
            if alignLen[contig]["-"] > alignLen[contig]["+"]:
                needsRevComp[fastaPrefix].add(contig)
    
    # Write output files with reverse-complemented sequences where needed
    for fastaFile in args.fastaFiles:
        fastaPrefix = os.path.basename(fastaFile).rsplit(".", maxsplit=1)[0]
        outputFile = os.path.join(args.outputDirectory, f"{fastaPrefix}.fasta")
        
        with open(fastaFile, "r") as fileIn, open(outputFile, "w") as fileOut:
            records = SeqIO.parse(fileIn, "fasta")
            for record in records:
                contig = record.id
                inputSeq = record.seq.reverse_complement() if contig in needsRevComp[fastaPrefix] else record.seq
                # Write multiline output
                fileOut.write(f">{contig}\n")
                sequence = str(inputSeq)
                fileOut.write("\n".join(
                    [
                        sequence[i:i+MULTILINE_LENGTH]
                        for i in range(0, len(sequence), MULTILINE_LENGTH)
                    ]
                ) + "\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
