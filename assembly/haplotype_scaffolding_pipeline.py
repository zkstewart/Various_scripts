#! python3
# haplotype_scaffolding_pipeline.py
# Script to take in a hifiasm haplotype assembly and a reference haplotype assembly
# to scaffold the hifiasm assembly using the reference assembly as a guide.

import os, argparse, sys, shutil

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
    for referenceFile in args.referenceGenome:
        if not os.path.isfile(referenceFile):
            raise FileNotFoundError(f"Reference genome file '{referenceFile}' not found.")
    
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
    # User input
    usage = """%(prog)s
    
    Contigs for hap1 and hap2 of the reference genome must be equivalently named.
    Defaults for --minQueryAlign and --minAlignLen are based on paf2dotplot.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputGenome",
                   required=True,
                   nargs=2,
                   help="Input hap1 and hap2 genome FASTA files requiring scaffolding")
    p.add_argument("-r", dest="referenceGenome",
                   required=True,
                   nargs=2,
                   help="Input hap1 and hap2 genome FASTA files for use as a reference")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (behavioural)
    p.add_argument("--minQueryAlign", dest="minQueryAlign",
                   required=False,
                   type=int,
                   help="""Optionally, specify the minimum cumulative alignment length for
                   a query to be considered a match (i.e., the sum of all alignments
                   that meet at least --minAlignLen); default == 400000""",
                   default=400000)
    p.add_argument("--minAlignLen", dest="minAlignLen",
                   required=False,
                   type=int,
                   help="""Optionally, specify the minimum alignment length to
                   consider a match; default == 10000""",
                   default=10000)
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
    
    args = p.parse_args()
    validate_args(args)
    
    # Set up the working directory structure
    referenceDir = os.path.join(args.outputDirectory, "reference")
    os.makedirs(referenceDir, exist_ok=True)
    
    inputDir = os.path.join(args.outputDirectory, "input")
    os.makedirs(inputDir, exist_ok=True)
    
    hap1Dir = os.path.join(args.outputDirectory, "hap1")
    os.makedirs(hap1Dir, exist_ok=True)
    
    hap2Dir = os.path.join(args.outputDirectory, "hap2")
    os.makedirs(hap2Dir, exist_ok=True)
    
    # Symlink files to the working directory locations
    ref1File = os.path.join(referenceDir, "hap1.fasta")
    if not os.path.exists(ref1File):
        os.symlink(args.referenceGenome[0], ref1File)
    
    ref2File = os.path.join(referenceDir, "hap2.fasta")
    if not os.path.exists(ref2File):
        os.symlink(args.referenceGenome[1], ref2File)
    
    input1File = os.path.join(inputDir, "hap1.fasta")
    if not os.path.exists(input1File):
        os.symlink(args.inputGenome[0], input1File)
    
    input2File = os.path.join(inputDir, "hap2.fasta")
    if not os.path.exists(input2File):
        os.symlink(args.inputGenome[1], input2File)
    
    # Index the reference files
    if not os.path.exists(f"{ref1File}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(ref1File, args.samtools)
    
    if not os.path.exists(f"{ref2File}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(ref2File, args.samtools)
    
    # Run minimap2 of input files against the reference files
    minimap2FlagName = os.path.join(inputDir, "minimap2_was_successful.flag")
    if not os.path.exists(minimap2FlagName):
        for queryIndex, queryFile in enumerate([input1File, input2File]):
            for refIndex, refFile in enumerate([ref1File, ref2File]):
                outputFileName = os.path.join(inputDir, f"i{queryIndex+1}_vs_r{refIndex+1}.paf")
                runner = ZS_AlignIO.Minimap2(queryFile, refFile, args.preset, args.minimap2, args.threads)
                runner.minimap2(outputFileName, force=True) # allow overwriting since the flag was not created
        open(minimap2FlagName, "w").close()
    else:
        print(f"Minimap2 alignment has already been performed; skipping.")
    
    # Parse minimap2 PAF files
    pafDict = {}
    hap1Contigs = set()
    hap2Contigs = set()
    refContigs = set()
    for refIndex in range(2):
        pafDict[refIndex+1] = {}
        for queryIndex in range(2):
            pafFile = os.path.join(inputDir, f"i{queryIndex+1}_vs_r{refIndex+1}.paf")
            with open(pafFile, "r") as fileIn:
                for line in fileIn:
                    # Extract relevant data
                    sl = line.rstrip("\r\n ").split("\t")
                    qid, qlen, qstart, qend, strand, tid, tlen, tstart, tend, \
                        numresidues, lenalign, mapq = sl[0:12]
                    
                    # Convert to integers
                    qlen, qstart, qend, tlen, tstart, tend, numresidues, lenalign, mapq = \
                        map(int, [qlen, qstart, qend, tlen, tstart, tend, numresidues, lenalign, mapq])
                    if queryIndex == 0:
                        hap1Contigs.add(qid)
                    else:
                        hap2Contigs.add(qid)
                    refContigs.add(tid)
                    
                    # Store the alignment length
                    pafDict[refIndex+1].setdefault(tid, {})
                    pafDict[refIndex+1][tid].setdefault(qid, 0)
                    pafDict[refIndex+1][tid][qid] += lenalign
    
    # 
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
