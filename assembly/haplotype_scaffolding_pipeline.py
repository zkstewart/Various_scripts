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
    if not os.path.isfile(args.inputGenome):
        raise FileNotFoundError(f"Input genome file '{args.inputGenome}' not found.")
    if not os.path.isfile(args.referenceGenome):
        raise FileNotFoundError(f"Reference genome file '{args.referenceGenome}' not found.")
    
    # Validate program discoverability
    if args.mmseqs2 == None:
        args.mmseqs2 = shutil.which("mmseqs2")
        if args.mmseqs2 == None:
            _not_found_error("mmseqs2")
    else:
        if not os.path.isfile(args.mmseqs2):
            _specified_wrong_error("mmseqs2", args.mmseqs2)
    
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
    usage = """%(prog)s ...
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputGenome",
                   required=True,
                   help="Input genome FASTA requiring scaffolding")
    p.add_argument("-r", dest="referenceGenome",
                   required=True,
                   help="Input genome FASTA for use as a reference")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (behavioural)
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
    p.add_argument("--mmseqs2", dest="mmseqs2",
                   required=False,
                   help="""Optionally, specify the mmseqs2 executable file
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
    
    hap1Dir = os.path.join(args.outputDirectory, "hap1")
    os.makedirs(hap1Dir, exist_ok=True)
    
    hap2Dir = os.path.join(args.outputDirectory, "hap2")
    os.makedirs(hap2Dir, exist_ok=True)
    
    # Symlink files to the working directory locations
    referenceFile = os.path.join(referenceDir, "reference.fasta")
    if not os.path.exists(referenceFile):
        os.symlink(args.referenceGenome, referenceFile)
    
    hap1File = os.path.join(hap1Dir, "hap1.fasta")
    if not os.path.exists(hap1File):
        os.symlink(args.inputGenome, hap1File)
    
    hap2File = os.path.join(hap2Dir, "hap2.fasta")
    if not os.path.exists(hap2File):
        os.symlink(args.inputGenome, hap2File)
    
    # Index the reference file
    if not os.path.exists(f"{referenceFile}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(referenceFile, args.samtools)
    
    # Run minimap2 of hap1/2 files against the reference file
    for hapDir, queryFile in zip([hap1Dir, hap2Dir], [hap1File, hap2File]):
        minimap2FlagName = os.path.join(hapDir, "minimap2_was_successful.flag")
        if not os.path.exists(minimap2FlagName):
            outputfile = os.path.join(hapDir, "minimap2.paf")
            runner = ZS_AlignIO.Minimap2(queryFile, referenceFile, args.preset, args.mmseqs2, args.threads)
            runner.minimap2(outputfile, force=True) # allow overwriting since the flag was not created
            open(minimap2FlagName, "w").close()
        else:
            print(f"Minimap2 alignment has already been performed for {queryFile}; skipping.")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
