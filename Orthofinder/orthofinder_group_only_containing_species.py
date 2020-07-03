#! python3
# orthofinder_group_only_containing_species
# Script to produce an abbreviated Orthogroups.csv file which optionally
# 1) returns groups ONLY containing the SOI, or 2) culls groups ONLY containing
# the SOI

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
produce a modified output file with 1) orthogroups which ONLY contain sequences from
the specified species of interest (SOI) or 2) all orthogroups OTHER THAN those resulting
from method (1); this means that the output will contain orthogroups where at least one
other species is present other than the SOI. Specify the behaviour (-b) a 1 or 2 to
correspond to the aforementioned scenarios.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-t", "--tsv", dest="orthogroups",
               help="Specify the orthogroup file")
p.add_argument("-o", "--output", dest="outputFileName",
               help="Output file name")
p.add_argument("-s", "--soi", dest="soi",
               help="Column name of the species of interest (SOI)")
p.add_argument("-b", "--behaviour", dest="behaviour", choices=["1", "2"],
               help="Specify program behaviour for orthogroups with the SOI")

args = p.parse_args()
validate_args(args)

# Parse Orthofinder file
header = None
with open(args.orthogroups, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if header == None:
                        soiIndex = None
                        header = line.rstrip('\r\n').strip('"').split('\t')[1:]         # [1:] gets rid of the blank space at the start of the file
                        if args.soi not in header:
                                print('I could not find the species of interest (SOI) you specified. These are the column headers I found:')
                                print(header)
                                print('And this is what you provided:')
                                print(args.soi)
                                print('Make sure your SOI is within the list above and try again.')
                                quit()
                        else:
                                soiIndex = header.index(args.soi) + 1 # Add +1 since the lines below have the OG# at the start, whereas we got rid of the blank space in the header here
                        fileOut.write(line)
                else:
                        sl = line.rstrip('\r\n').split('\t')
                        # Behaviour handling
                        behaviourCondition = True
                        for i in range(1, len(sl)):
                                # If the SOI column contains sequences
                                if i == soiIndex and sl[i] != "":
                                        continue
                                # If the SOI column does not contain sequences
                                elif i == soiIndex and sl[i] == "":
                                        behaviourCondition = False
                                        break
                                # If the column is other than the SOI and contains sequences
                                elif i != soiIndex and sl[i] != "":
                                        behaviourCondition = False
                                        break
                                # If the sequence is other than the SOI and does not contain sequences
                                else:
                                        continue
                        if args.behaviour == "1" and behaviourCondition == True:
                                fileOut.write(line)
                        elif args.behaviour == "1" and behaviourCondition == False:
                                continue
                        elif args.behaviour == "2" and behaviourCondition == True:
                                continue
                        elif args.behaviour == "2" and behaviourCondition == False:
                                fileOut.write(line)

# Done!
print('Program completed successfully!')
