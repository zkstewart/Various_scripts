#! python 3
# mafftAlign.py
# A script to receive one or more FASTA files, and perform MAFFT
# alignment for each file.

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from Function_packages import ZS_Utility, ZS_AlignIO

def validate_args(args):
    # Validate input files
    args.fileNames = []
    for value in args.input:
        value = os.path.abspath(value)
        
        if os.path.isfile(value):
            args.fileNames.append(value)
        elif os.path.isdir(value):
            for file in os.listdir(value):
                if any([ file.endswith(ext) for ext in args.fastaExtensions ]):
                    args.fileNames.append(os.path.join(value, file))
        else:
            print(f"I am unable to locate the file or directory ({value})")
            print("Make sure you've typed the file name or location correctly and try again.")
            quit()
    if args.fileNames == []:
        print("No FASTA files with the indicated extensions were found in any input locations.")
        print("Make sure you've typed the file name or location correctly and try again.")
        quit()
    
    # Validate MAFFT arguments
    if args.mafft is None:
        args.mafft = ZS_Utility.wsl_which("mafft")
        if args.mafft is None:
            print(f"ERROR: 'mafft' not discoverable in your system PATH and was not specified as an argument.")
            quit()
    else:
        if not os.path.isfile(args.mafft):
            print(f"ERROR: 'mafft' was not found at the location indicated ('{args.mafft}')")
            quit()
    if args.maxiterate < 0:
        print(f"ERROR: --maxiterate must be a non-negative integer.")
        quit()
    if args.threads < 1:
        print(f"ERROR: --threads must be an integer greater than or equal to 1.")
        quit()
    
    # Handle file output
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def main():
    usage = """%(prog)s will receive one or more input FASTA files, specified as individual file
    paths or as directories to scan through for files with the indicated suffix(es). It will then
    perform a multiple sequence alignment using MAFFT for each file, and write the output to the
    specified output directory. Some notes for use:
    
    1) --codons will align by codon, by first translating sequences into their +ve strand, frame 0
    protein sequences. The input sequences should thus be CDS' that match this assumption.
    2) --maxiterate can be set to a very high value, but the MAFFT documentation suggests that most
    improvement occurs within the first few iterations, so take that into consideration (e.g.,
    maybe set it to 5 or 10 if you're in a hurry, or 100+ if it needs to be perfection).
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="input",
                   required=True,
                   nargs="+",
                   help="Specify the location of FASTA files or directories to scan for FASTA files")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory for aligned FASTA file(s)")
    # Opts (MAFFT)
    p.add_argument("--mafft", dest="mafft",
                   required=False,
                   help="""Optionally, specify the mafft executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--alg", dest="algorithm",
                   required=False,
                   choices=["auto", "einsi", "linsi", "ginsi", "fftns1", "fftns2", "fftnsi"],
                   help="""Optionally, specify the MAFFT algorithm to use for alignment;
                   default == 'auto'""",
                   default=None)
    p.add_argument("--maxiterate", dest="maxiterate",
                   required=False,
                   type=int,
                   help="""Optionally, specify the --maxiterate value to use for alignment;
                   default == 0""",
                   default=0)
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify the number of threads to use for alignment;
                   default == 1""",
                   default=1)
    # Opts (program behaviour)
    p.add_argument("--codons", dest="codons",
                   required=False,
                   action="store_true",
                   help="""Optionally, specify this flag if your input sequences are nucleotides
                   that you would like to align by codon.""",
                   default=False)
    p.add_argument("--suffix", dest="suffix",
                   required=False,
                   help="""Optionally, specify a suffix to add to each FASTA file name prior to the
                   file extension to indicate that it has been aligned; default == '', which
                   means output files will be named identically to each input.""",
                   default="")
    p.add_argument("--fastaExtensions", dest="fastaExtensions",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more file extensions to
                   use for identifying FASTA files in directories; default == several
                   commonly used extensions.""",
                   default=[".fasta", ".fa", ".fas", ".fna",
                            ".faa", ".cds", ".nucl", ".trans"])
    
    args = p.parse_args()
    validate_args(args)
    
    aligner = ZS_AlignIO.MAFFT(args.mafft, algorithm=args.algorithm,
                               thread=args.threads, maxiterate=args.maxiterate)
    
    for fastaFile in args.fileNames:
        # Prepare output file name
        basename, fileExtension = os.path.splitext(fastaFile)
        basePrefix = os.path.basename(basename)
        outName = os.path.join(args.outputDirectory, basePrefix + args.suffix + fileExtension)
        
        # Skip if file already exists
        if os.path.isfile(outName):
            print(f"Output file '{outName}' already exists; skipping this file.")
            continue
        
        # Perform MAFFT alignment
        if args.codons:
            numSeqs = sum([1 for _ in SeqIO.parse(fastaFile, "fasta")])
            
            resultFASTA_obj = aligner.align_as_protein(
                fastaFile,
                strands=[1 for _ in range(numSeqs)],
                frames=[0 for _ in range(numSeqs)]
            )
        else:
            resultFASTA_obj = aligner.align(fastaFile)
        
        # Write output file
        resultFASTA_obj.write(outName, asAligned=True)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
