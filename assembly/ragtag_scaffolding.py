#! python3
# ragtag_scaffolding.py
# Script to facilitate easier use of RagTag to scaffold contigs
# particularly for tricky cases that haplotype_scaffolding_pipeline
# did not solve

import os, argparse, sys, shutil, subprocess, platform
from Bio import SeqIO

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
    
    # Validate numeric arguments
    if args.threads < 1:
        raise ValueError(f"--threads must be greater than or equal to 1.")
    
    # Validate program discoverability
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
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

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

## Main
def main():
    MULTILINE_LENGTH = 70
    
    # User input
    usage = """%(prog)s automates the scaffolding of two or more input contigs
    against a specified reference genome using RagTag.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputGenome",
                   required=True,
                   nargs="+",
                   help="Input one or more genome FASTA files with contigs requiring scaffolding")
    p.add_argument("-is", dest="inputSequences",
                   required=True,
                   nargs="+",
                   help="Input one or more contig IDs to grab from the input FASTAs and scaffold")
    p.add_argument("-r", dest="referenceGenome",
                   required=True,
                   help="Input genome FASTA file for use as a reference")
    p.add_argument("-rs", dest="referenceSequence",
                   required=True,
                   help="Input contig to use during scaffolding")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (ragtag)
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify the number of threads to use for ragtag;
                   default == 1""",
                   default=1)
    # Opts (programs)
    p.add_argument("--ragtag", dest="ragtag",
                   required=False,
                   help="""Optionally, specify the ragtag.py file
                   if it is not discoverable in the path""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Load in sequences as a FastaCollection object
    inputRecords = ZS_SeqIO.FastaCollection(args.inputGenome)
    refRecords = ZS_SeqIO.FastaCollection([args.referenceGenome])
    
    # Generate scaffolds for each haplotype chromosome
    outputFileName = os.path.join(args.outputDirectory, f"{args.referenceSequence}.fasta")
    outputFlagName = os.path.join(args.outputDirectory, f"{args.referenceSequence}.is.ok.flag")
    if not os.path.exists(outputFlagName):
        # Create a working directory for the files
        chrDir = os.path.join(args.outputDirectory, args.referenceSequence)
        os.makedirs(chrDir, exist_ok=True)
        
        # Write the input sequences to a file
        "Use a flag here to allow for hacky resume behaviours"
        rawSequencesFile = os.path.join(chrDir, f"{args.referenceSequence}.input.fasta")
        rawSequencesFlag = os.path.join(chrDir, f"{args.referenceSequence}.input.fasta.is.ok.flag")
        if not os.path.exists(rawSequencesFlag):
            with open(rawSequencesFile, "w") as fileOut:
                for inputContig in args.inputSequences:
                    inputSeq = inputRecords[inputContig]
                    fileOut.write(f">{inputContig}\n{str(inputSeq)}\n")
            open(rawSequencesFlag, "w").close()
        
        # Write the reference sequence to a file
        refSequenceFile = os.path.join(chrDir, f"{args.referenceSequence}.ref.fasta")
        refSequenceFlag = os.path.join(chrDir, f"{args.referenceSequence}.ref.fasta.is.ok.flag")
        if not os.path.exists(refSequenceFlag):
            with open(refSequenceFile, "w") as fileOut:
                refSeq = refRecords[args.referenceSequence]
                fileOut.write(f">{args.referenceSequence}\n{str(refSeq)}\n")
            open(refSequenceFlag, "w").close()
        
        # Run ragtag to scaffold the sequences
        run_ragtag(rawSequencesFile, refSequenceFile, os.path.join(chrDir, "ragtag_output"),
                    args.ragtag, threads=1)
        
        # Rewrite the output file to the haplotype directory
        ragtagOutputFile = os.path.join(chrDir, "ragtag_output", f"ragtag.scaffold.fasta")
        seqID, numIDs = find_longest_seq(ragtagOutputFile) # detect if ragtag couldn't scaffold them all together
        if numIDs > 1:
            print(f"{numIDs} IDs were found in the ragtag output '{ragtagOutputFile}'; " + 
                    f"we will just use the longest sequence '{seqID}'")
        
        with open(ragtagOutputFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
            records = SeqIO.parse(fileIn, "fasta")
            for record in records:
                if record.id == seqID:
                    fileOut.write(f">{args.referenceSequence}\n")
                    fileOut.write(str(record.seq) + "\n")
        
        # Create a flag to indicate that the file was created
        open(outputFlagName, "w").close()
    else:
        print(f"contig '{args.referenceSequence}' for '{os.path.basename(args.outputDirectory)}' already exists; skipping.")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
