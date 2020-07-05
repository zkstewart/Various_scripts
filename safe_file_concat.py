#! python3
# safe_file_concat
# Self-explanatory script to concatenate files without having to worry about
# newline endings and such.

import argparse, os

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure no None arguments exist
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument was not specified. Fix this and try again.')
                        quit()
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

##### USER INPUT SECTION
usage = """%(prog)s is a simple program to concatenate files together, ensuring
that proper newline spacing is maintained. It is especially useful for FASTA files
where some have newline endings and others do not.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="inputLocation", nargs="+",
               help="Input director(y/ies) and/or file location(s)")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name.")

args = p.parse_args()
validate_args(args)

# Find files depending on how inputLocation was specified
input_files = []
for inputLocation in args.inputLocation:
        if os.path.isdir(inputLocation):
                # Scan through files and detect our files of interest
                for f in os.listdir(inputLocation):
                        input_files.append(os.path.join(inputLocation, f))
        else:
                # Add the absolute path of the file to our list
                input_files.append(os.path.abspath(inputLocation))

# Produce output file
with open(args.outputFileName, 'w') as file_out:
        for f in input_files:
                with open(f, 'r') as input_file:
                        for line in input_file:
                                # Write line as posix standard
                                file_out.write(line.rstrip('\r\n') + '\n')
