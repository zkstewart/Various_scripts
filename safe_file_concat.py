#! python3
# safe_file_concat
# Self-explanatory script to concatenate files without having to worry about
# newline endings and such.

import argparse, os

# Define functions for later use
## Validate arguments
def validate_args(args):
    # Validate the input file locations depending on type of input
    for inputLocation in args.inputLocation:
        # Check that specified value is a path or file
        if not os.path.isdir(inputLocation) and not os.path.isfile(inputLocation):
            print('The provided value "' + inputLocation + '" is not a directory or file; fix your input and try again.')
            quit()
    # Handle file overwrites
    if os.path.isdir(args.outputFileName):
        print(args.outputFileName + ' is a directory. Delete/move/rename this directory and run the program again.')
        quit()
    elif os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()
    elif os.path.dirname(args.outputFileName) != '':
        if not os.path.isdir(os.path.dirname(args.outputFileName)):
            print(args.outputFileName + ' is to be created in a non-existent directory. Make this directory or rename your output file and run the program again.')
            quit()

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s is a simple program to concatenate files together, ensuring
    that proper newline spacing is maintained. It is especially useful for FASTA files
    where some have newline endings and others do not.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", "-input", nargs="+",
                   required=True,
                   dest="inputLocation",
                   help="Input director(y/ies) and/or file location(s)")
    p.add_argument("-o", "-output",
                   required=True,
                   dest="outputFileName",
                   help="Output file name.")
    
    args = p.parse_args()
    validate_args(args)
    
    # Find files depending on how inputLocation was specified
    input_files = []
    for inputLocation in args.inputLocation:
        # Handle dir values
        if os.path.isdir(inputLocation):
            for f in os.listdir(inputLocation):
                f_path = os.path.join(inputLocation, f)
                if os.path.isfile(f_path):
                    input_files.append(os.path.abspath(f_path))
        # Handle file values [validate_args has already proven that if it isn't a dir, it's a file]
        else:
            input_files.append(os.path.abspath(inputLocation))
    
    # Produce output file
    with open(args.outputFileName, 'w') as file_out:
        for f in input_files:
            with open(f, 'r') as input_file:
                for line in input_file:
                    l = line.rstrip("\r\n")
                    if l != "":
                        file_out.write(l + "\n") # write line as posix standard
    
    # Notify user of successful completion
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
