#! python3
# orthofinder_group_statistics
# Script to calculate some basic group statistics from the Orthogroups.csv file
# output by OrthoFinder. These statistics currently are limited to the similarity
# in copy number within orthogroups.

# Load packages
import os, argparse
from decimal import Decimal

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

### USER INPUT
usage = """%(prog)s will parse the Orthogroups.csv file output by Orthofinder and 
produce a short tab-delimited text file indicating the number of species who have
orthogroup occurrence counts (or copy numbers) greater than, less than, or equal to a specified
species of interest (SOI). It is expected that the SOI should correspond to a column
name in the input file. Note that, for the final section under the commented (#) line
the leftmost species is treated as the "original" and the rightmost is the "new" which
should be how percentage increases or decreases as interpreted. 
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-c", "-csv", dest="orthogroups",
               help="Specify the orthogroup file")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name")
p.add_argument("-s", "-soi", dest="soi", type=str,
               help="Column name of the species of interest (SOI) to use as a reference")

args = p.parse_args()
validate_args(args)

# Parse orthogroup file
orthoDict, header = parse_orthogroups_csv(args.orthogroups)

# Count the copy number present in each orthogroup for later manipulation
orthoDictNum = {}
for key, value in orthoDict.items():
        if key in header:       # This lets us skip the alternative indexing system values in orthoDict
                continue
        orthoDictNum[key] = {}
        # Count values
        for subkey, subvalue in value.items():
                orthoDictNum[key][subkey] = len(subvalue)

# Compare orthogroup(s) to assess similarity in copy number
orthoStats = {}
for key, value in orthoDictNum.items():
        # Extract values for handling
        enumeratedValue = list(enumerate(value.items()))
        # Instantiate a dictionary structure for storing results
        if orthoStats == {}:
                for x in range(len(enumeratedValue)):
                        for y in range(len(enumeratedValue)):
                                if x == y:
                                        continue
                                # Extract values for handling
                                xval = enumeratedValue[x]
                                yval = enumeratedValue[y]
                                # Equal to count
                                orthoStats[xval[1][0] + '_=_' + yval[1][0]] = 0
                                # Greater than count
                                orthoStats[xval[1][0] + '_>_' + yval[1][0]] = 0
                                # Less than count
                                orthoStats[xval[1][0] + '_<_' + yval[1][0]] = 0
        # Tally each comparison in our dictionary structure
        for x in range(len(enumeratedValue)):
                for y in range(len(enumeratedValue)):
                        if x == y:
                                continue
                        # Extract values for handling
                        xval = enumeratedValue[x]
                        yval = enumeratedValue[y]
                        # Compare group sizes and tally results
                        if xval[1][1] == yval[1][1]:
                                orthoStats[xval[1][0] + '_=_' + yval[1][0]] += 1
                        elif xval[1][1] > yval[1][1]:
                                orthoStats[xval[1][0] + '_>_' + yval[1][0]] += 1
                        else:
                                orthoStats[xval[1][0] + '_<_' + yval[1][0]] += 1

# Compute additional proportion values
orthoProp = {}
for key, value in orthoStats.items():
        proportion = value / len(orthoDictNum)
        proportion = round(Decimal(value / len(orthoDictNum)) * 100, 2)
        orthoProp[key] = float(proportion)

# Compare group(s) to reference SOI to assess which one is more/less similar
orthoComps = []
for symbol in ['<', '>', '=']:
        compKeys = []
        compValues = []
        # Extract values that are $symbol the SOI's value
        for key, value in orthoStats.items():
                if '_' + symbol + '_' + args.soi in key:
                        compKeys.append(key.split('_' + symbol + '_')[0])       # This gives us just the column name that is $symbol the SOI's value
                        compValues.append(value)
        # Compare these extracted values
        for x in range(len(compKeys)):
                for y in range(len(compKeys)):
                        if x == y:
                                continue
                        xminusy = compValues[x] - compValues[y]
                        # Calculate percentage increase/decrease
                        if xminusy != 0:
                                percent = abs(((compValues[x] - compValues[y])/compValues[y]) * 100)
                        else:
                                percent = 0.00
                        percent = round(Decimal(percent), 2)
                        # Store as a string summary
                        if xminusy > 0:
                                orthoComps.append(compKeys[x] + '_vs_' + compKeys[y] + '\t' + symbol + args.soi + '\t' + str(xminusy) + '\t' + str(percent) + '% increase')
                        elif xminusy < 0:
                                orthoComps.append(compKeys[x] + '_vs_' + compKeys[y] + '\t' + symbol + args.soi + '\t'  + str(xminusy) + '\t' + str(percent) + '% decrease')
                        else:
                                orthoComps.append(compKeys[x] + '_vs_' + compKeys[y] + '\t' + symbol + args.soi + '\t'  + str(xminusy) + '\tno change')

# Tabulate results in human-readable format
with open(args.outputFileName, 'w') as fileOut:
        # Header information
        fileOut.write('#%(prog)s output statistics file for: ' + os.path.abspath(args.orthogroups) + '\n')
        fileOut.write('#total orthogroups = ' + str(len(orthoDictNum)) + '\n')
        # Write counts details
        fileOut.write('#column_comparison\tnumber_of_groups\t%_of_total\n')
        for key, value in orthoStats.items():
                count = value
                proportion = orthoProp[key]
                fileOut.write(key + '\t' + str(count) + '\t' + str(proportion) + '\n')
        # Write group comparison details
        fileOut.write('#column_comparison\tcomparison_direction\tnumeric_difference\tpercentage_difference\n')
        for line in orthoComps:
                fileOut.write(line + '\n')

# Done!
print('Program completed successfully!')
