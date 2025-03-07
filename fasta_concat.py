#! python3
# fasta_concat
# Simple script to concatenate FASTA files together
# by specifying individual files and/or directories
# containing FASTA files.

import argparse, os
from Bio import SeqIO

def ends_with_suffix(file, suffixes):
    return any([file.endswith(s) for s in suffixes])

def validate_args(args):
    # Validate input file locations depending on type of input
    args.files = []
    for inputLocation in args.inputLocations:
        inputLocation = os.path.abspath(inputLocation)
        
        # Handle directories
        if os.path.isdir(inputLocation):
            foundFiles = False
            for f in os.listdir(inputLocation):
                file = os.path.join(inputLocation, f)
                if os.path.isfile(file) and ends_with_suffix(file, args.fileSuffixes):
                    args.files.append(file)
                    foundFiles = True
            if not foundFiles:
                raise FileNotFoundError(f"'{inputLocation}' does not contain any files.")
        
        # Handle files
        elif os.path.isfile(inputLocation): # don't need to check suffix here
            args.files.append(inputLocation)
        else:
            raise FileNotFoundError(f"'{inputLocation}' is not a valid file or directory.")
    
    # Validate numeric arguments
    if args.minimumLength < 0:
        raise ValueError(f"--minimum '{args.minimumLength}' must be an integer >= 0")
    if args.multilineSize < 0:
        raise ValueError(f"--multiline '{args.multilineSize}' must be an integer >= 0")
    
    # Validate output file
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o output file '{args.outputFileName}' already exists!")
    args.outputFileName = os.path.abspath(args.outputFileName)
    
    if not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"-o parent directory '{os.path.dirname(args.outputFileName)}' does not exist!")

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s concatenates FASTA files together, ensuring
    that proper newline spacing is maintained. It offers some filtering
    ability on sequence length, and can also multiline the output.
    """
    # Required arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputLocations",
                   required=True,
                   nargs="+",
                   help="Location of files or directories containing files to concatenate")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write concatenated outputs to")
    # Optional arguments
    p.add_argument("--minimum", dest="minimumLength",
                   required=False,
                   type=int,
                   help="""Optionally, specify the minimum length of a sequence to be
                   included in the output; default == 0 for no minimum""",
                   default=0)
    p.add_argument("--multiline", dest="multilineSize",
                   required=False,
                   type=int,
                   help="""Optionally, specify the size of each line in the output;
                   default == 0 for no multilining""",
                   default=0)
    p.add_argument("--suffixes", dest="fileSuffixes",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more file endings you want
                   to include in the concatenation; default == 
                   ['.fasta', '.fas', '.fa', '.fna', '.faa', '.aa', '.prot', '.nucl', '.trans', '.mrna']""",
                   default=['.fasta', '.fas', '.fa', '.fna', '.faa', '.aa', '.prot', '.nucl', '.trans', '.mrna'])
    
    args = p.parse_args()
    validate_args(args)
    
    # Produce output file
    with open(args.outputFileName, "w") as fileOut:
        for f in args.files:
            for record in SeqIO.parse(f, "fasta"):
                sequence = str(record.seq)
                
                # Skip sequences that are too short
                if len(sequence) < args.minimumLength:
                    continue
                
                # Format sequence as multiline if requested
                if args.multilineSize > 0:
                    sequence = "\n".join([ sequence[i:i+args.multilineSize] for i in range(0, len(sequence), args.multilineSize) ])
                
                # Write sequence to output file
                fileOut.write(f">{record.description}\n{sequence}\n")
    
    # Notify user of successful completion
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
