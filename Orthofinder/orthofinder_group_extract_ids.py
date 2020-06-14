#! python3
# orthofinder_group_extract_ids
# Script to extract sequence IDs corresponding to groups for
# each respective file

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure all arguments are specified
        for key, value in vars(args).items():
                if key != "soi" and key != "outputFileName":
                        if value == None or value == []:
                                print(key + ' argument was not specified; fix your input and try again.')
                                quit()
                if key == "outputFileName" and args.one_file:
                        if value == None or value == []:
                                print(key + ' argument was not specified; fix your input and try again.')
                                quit()
        # Validate Orthogroup file location
        if not os.path.isfile(args.orthogroups):
                print('I am unable to locate the Orthogroups file (' + args.orthogroups + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Ensure that the output location is sensible
        if args.outputFileName != None:
                if os.path.exists(args.outputFileName):
                        print('The specified output file name "' + args.outputFileName + '" already exists; fix your input and try again.')
                        quit()

## OrthoFinder-related
def parse_orthogroups_csv(orthogroupsFile):
        # Setup
        header = None
        orthoDict = {}
        # Main function
        with open(orthogroupsFile , 'r') as orthoFile:
                for line in orthoFile:
                        if header == None:
                                header = line.rstrip('\r\n').strip('"').split('\t')[1:]                 # [1:] gets rid of the blank space at the start of the file
                                # Establish subdicts for later indexing
                                for i in range(len(header)):
                                        orthoDict[header[i].strip('"')] = {}
                        else:
                                sl = line.rstrip('\r\n').strip('"').replace('""', '"').split('\t')      # For some reason sequence IDs with " in them have it replaced with ""
                                for i in range(1, len(sl)):
                                        sl[i] = sl[i].strip('"').split(', ')    # OrthoFinder separates sequence IDs like so; this works out since OrthoFinder removes commas from the IDs (which we need to handle during the parse)
                                        if sl[i] == ['']:
                                                sl[i] = []
                                        # Index species by orthogroup
                                        if i == 1:
                                                orthoDict[sl[0]] = {}
                                        orthoDict[sl[0]][header[i-1]] = sl[i]   # i-1 since sl has the orthogroup ID as its first value
                                # Index orthogroup by species
                                for i in range(len(header)):
                                        orthoDict[header[i]][sl[0]] = sl[i+1]   # i+1 since the first value in sl is the orthogroup ID
        return orthoDict, header

#### USER INPUT SECTION
usage = """%(prog)s is intended to receive an Orthogroups.tsv file which has
been pre-processed by orthofinder_group_containing_species.py. It returns
text files containing the sequence members for each group. Optionally, only groups
that contain one or more species of interest (SOI) may be returned.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-c", "-csv", dest="orthogroups",
               help="Specify the orthogroup file")
p.add_argument("-s", "-soi", dest="soi", nargs="+",
               help="Column name(s) of the species of interest (SOI)")
p.add_argument("--one_file", dest="one_file", action="store_true", default=False,
               help="Specify if you do not want files to be created for each group")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Specify output file name (only relevant if --one_file is set)")

args = p.parse_args()
validate_args(args)

# Parse orthogroup file
orthoDict, header = parse_orthogroups_csv(args.orthogroups)

# Check SOIs for sensibility
if args.soi != []:
        for soi in args.soi:
                if soi not in header:
                        print(soi + " is not a recognised species")
                        print("Valid SOIs are:")
                        print(header)
                        quit()

# Main output loop
for key, value in orthoDict.items():
        if key in header: # This lets us skip the alternative indexing system values in orthoDict
                continue
        # Continue condition: Check for SOI(s) if relevant
        if args.soi != []:
                skip = True
                for soi in args.soi:
                        if value[soi] != []:
                                skip = False
                                break
                if skip == True:
                        continue
        # Make a single list for sequence IDs
        seqids = []
        for subkey, subvalue in value.items():
                seqids += subvalue
        # Make list writeable to file
        seqids = '\n'.join(seqids)
        # Write output
        if args.one_file:
                with open(args.outputFileName, "a") as fileOut:
                        fileOut.write(seqids + '\n')
        else:
                with open(key + ".txt", "w") as fileOut:
                        fileOut.write(seqids + '\n')

# Done!
print('Program completed successfully!')
