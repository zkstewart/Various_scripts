#! python3
# orthofinder_group_containing_species
# Script to produce an abbreviated Orthogroups.csv file which contains only
# orthogroups which contain the species of interest

# Load packages
import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure all arguments are specified
        for key, value in vars(args).items():
                if value == None or value == []:
                        print(key + ' argument was not specified; fix your input and try again.')
                        quit()
        # Validate Orthogroup file location
        if not os.path.isfile(args.orthogroups):
                print('I am unable to locate the Orthogroups file (' + args.orthogroup + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

### USER INPUT
usage = """%(prog)s will parse the Orthogroups.csv file output by Orthofinder and 
produce a modified output file with only orthogroups which contain sequences from
the specified species of interest (SOI). Specification of SOI is in the form of
providing the name of the column(s) you want; if specifying more than one column,
separate each value with a comma e.g., -s species1,species2
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-c", "-csv", dest="orthogroups",
               help="Specify the orthogroup file")
p.add_argument("-p", "-output", dest="outputFileName",
               help="Output file name")
p.add_argument("-s", "-soi", dest="soi", nargs="+",
               help="Column name(s) of the species of interest (SOI)")

args = p.parse_args()
validate_args(args)

# Parse Orthofinder file
header = None
with open(args.orthogroups, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if header == None:
                        soiIndices = []
                        header = line.rstrip('\r\n').strip('"').split('\t')[1:]         # [1:] gets rid of the blank space at the start of the file
                        for soi in args.soi:
                                if soi not in header:
                                        print('I could not find the species of interest (SOI) you specified. These are the column headers I found:')
                                        print(header)
                                        print('And this is what you provided:')
                                        print(soi)
                                        print('Make sure your SOI is within the list above and try again.')
                                        quit()
                                else:
                                        soiIndices.append(header.index(soi))
                        fileOut.write(line)
                else:
                        sl = line.rstrip('\r\n').split('\t')
                        for index in soiIndices:
                                if sl[index+1] != '':           # Add +1 since these lines have the OG# at the start, whereas we got rid of the blank space in the header that corresponds to this column
                                        fileOut.write(line)
