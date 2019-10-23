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
        if len(args.inputLocation) == 1:
                # Check that specified value is a path
                if not os.path.isdir(args.inputLocation[0]):
                        print('One value was provided for -i, which means you should have provided a directory containing FASTA files.')
                        print('The provided value "' + args.inputLocation[0] + '" is not a directory; either it does not exist or it is a file; fix your input and try again.')
                        quit()
        else:
                # Check that the specified values are files
                for file in args.inputLocation:
                        if not os.path.isfile(file):
                                print('Multiple values were provided for -i, which means you should have provided the location of individual FASTA files.')
                                print('The provided value "' + file + '" is not a file; either it does not exist or it is a directory; fix your input and try again.')
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
               help="Input a single directory or multiple file locations")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name.")

args = p.parse_args()
validate_args(args)

# Find files depending on how inputLocation was specified
input_files = []
if len(args.inputLocation) == 1:
        # Scan through files and detect our files of interest
        for file in os.listdir(args.inputLocation[0]):
                input_files.append(os.path.join(args.inputLocation[0], file))
else:
        # Make sure that the provided files all exist
        for file in args.inputLocation:
                if not os.path.isfile(file):
                        print('Input file "' + file + '" either does not exist or is not a file.')
                        print('Make sure you spelled the file name/location correctly and try again.')
                        quit()
                input_files.append(os.path.abspath(file))

# Produce output file
with open(args.outputFileName, 'w') as file_out:
        for f in input_files:
                with open(f, 'r') as input_file:
                        for line in input_file:
                                if not line.endswith('\n'):
                                        file_out.write(line + '\n')
                                else:
                                        file_out.write(line)
