#! python3
# run_signalp.py
# Script to take a FASTA file input and run SignalP
# on it for signal peptide predictions

import os, argparse, sys
from threading import Thread

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Function_packages import ZS_SignalPIO, ZS_SeqIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.fastaFile):
        print(f'I am unable to locate the input FASTA file ({args.fastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate program location
    if not os.path.isfile(args.signalpExe):
        print(f'I am unable to locate the SignalP exe file ({args.signalpExe})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric parameters
    if args.threads < 1:
        print("threads must be a positive integer")
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

class SignalPThread(Thread):
    '''
    This provides a modified Thread which allows for the output of signalP
    to be stored locally within this object. Necessary since we can't otherwise
    get return values from threads.
    '''
    def __init__(self, fastaFile, signalpExe, organism, cygwinDir=None):
        Thread.__init__(self)
        
        self.fastaFile = fastaFile
        self.signalpExe = signalpExe
        self.organism = organism
        self.cygwinDir = cygwinDir
        self.output = None
    
    def run(self):
        # Create and configure signalP handler object
        sigp = ZS_SignalPIO.SignalP(self.fastaFile, self.signalpExe, self.cygwinDir)
        sigp.organism = self.organism
        sigp.clean = True
        
        # Run prediction
        self.output = sigp.signalp()

def merge_dictionaries(dict1, dict2):
    '''
    This function assumes the two dictionaries have no keys in common. If
    there is, whatever is in dict1 will be kept.
    '''
    dict1.update(dict2)
    return dict1

## Main
def main():
    # User input
    usage = """%(prog)s reads in a FASTA file and generates an output TSV indicating
    the signal peptide predictions made. It's assumed your FASTA contains protein
    sequences; nucleotide input may lead to unpredictable outcomes.
    
    This script relies on the ZS_SignalPIO package, which at this time attempts to
    handle SignalP versions 4, 5, and 6. However, version 4 is non-functional, and
    version 6 is only functional on Linux (not WSL!).
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="fastaFile",
                   required=True,
                   help="Input FASTA file containing protein sequences.")
    p.add_argument("-s", dest="signalpExe",
                   required=True,
                   help="Input full path to the signalP executable")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for signalP predictions")
    # Optional
    p.add_argument("--org", dest="organism",
                   required=False,
                   choices=["euk", "eukarya", "gram+", "gram-", "arch", "other"],
                   help="""Optionally, specify which organism type to search for
                   (default == "euk"; make sure to use the right term for your
                   version of signalP)""",
                   default="euk")
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Specify how many threads to operate with; default == 1""",
                   default=1)
    p.add_argument("--cygwinDir", dest="cygwinDir",
                   required=False,
                   help="""If you intend on running signalP 4 on Windows, you must provide
                   this value; WSL does not work, hence we must run it through Cygwin""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Handle chunking for multithreading
    if args.threads > 1:
        # Parse in FASTA as a ZS_SeqIO.FASTA object
        FASTA_obj = ZS_SeqIO.FASTA(args.fastaFile)
        
        # Generate temporary file prefix
        outputLocation = os.path.dirname(args.outputFileName)
        tmpFilePrefix = ZS_SeqIO.Conversion.get_hash_for_input_sequences(FASTA_obj)
        outputFilePrefix = os.path.join(outputLocation, tmpFilePrefix)
        
        # Chunk FASTA into separate files for multi-threading
        fastaFiles = FASTA_obj.write_as_chunks(outputFilePrefix, args.threads)
    else:
        # Use file as-is with a single thread
        fastaFiles = [args.fastaFile]
    
    # Run SignalP in multiple threads
    processing = []
    for i in range(args.threads):
        fastaFile = fastaFiles[i]
        workerThread = SignalPThread(fastaFile, args.signalpExe, args.organism, args.cygwinDir)
        processing.append(workerThread)
        workerThread.start()
    
    # Gather results
    sigpResultsDict = {}
    for workerThread in processing:
        # Wait for thread to end ...
        workerThread.join()
                
        # ... then merge dictionary outputs together
        sigpResultsDict.update(workerThread.output)
    
    # Create output file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#sequence_id\tstart\tend\n")
        # Content lines
        for seqid, coords in sigpResultsDict.items():
            start, end = coords
            fileOut.write(f"{seqid}\t{start}\t{end}\n")
    
    # Clean up chunked FASTA files if necessary
    if args.threads > 1:
        for fastaFile in fastaFiles:
            os.unlink(fastaFile)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
