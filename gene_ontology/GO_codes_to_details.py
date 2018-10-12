#! python3

import os, argparse
from goatools import obo_parser

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputTextFile):
                print('I am unable to locate the input text file (' + args.inputTextFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.oboFile):
                print('I am unable to locate the go-basic.obo file (' + args.oboFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
                quit()

def text_go_parse(textFile):
        # Setup
        codeDict = {}
        # Main function
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        line = line.rstrip('\r\n')
                        if line == '.':
                                continue
                        else:
                                codes = line.split('; ')
                                for code in codes:
                                        codeDict[code] = ''
        return codeDict

#### USER INPUT SECTION
usage = """This program will read in a text file listing GO codes separated with '; '
characters alongside the go-basic.obo file and produces an output text file listing the
names for these ontologies. '.' can be used to mark lines with no GO results, and this
will be carried over into the output file.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="inputTextFile",
                   help="Input text file name listing GO codes (with '.' acting as a marker for no codes).")
p.add_argument("-g", "-goobo", dest="oboFile",
                   help="Input go-basic.obo file.")
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output text file name.")

args = p.parse_args()
validate_args(args)

# Parse the input text file to identify which GO codes are listed
codeDict = text_go_parse(args.inputTextFile)

# Parse .obo file
go = obo_parser.GODag(args.oboFile)

# Hardcoded GO replacements for problem entries
obsoletedGOs = ['GO:0097034', 'GO:0097033']   # I can't find any suitable replacements for these terms
replacedGOs = {'GO:0004871': 'GO:0007165', 'GO:0004702': 'GO:0007165', 'GO:0005057': 'GO:0007165', # These replacements were made manually as the keys are not present within the go-basic.obo file downloaded 24/05/2018
               'GO:0042993': 'GO:0042307', 'GO:0042991': 'GO:0006606', 'GO:0044376': 'GO:0031503', # The idmapping_selected.tab file did contain these keys; file version was dated 25/04/2018
               'GO:1990022': 'GO:0031503', 'GO:1904721': 'GO:1903895', 'GO:0042992': 'GO:0042308', # Additional modifications were made to this code 01/10/2018
               'GO:0030819': 'GO:0030816', 'GO:0051436': 'GO:0000278', 'GO:0051442': 'GO:0051321',
               'GO:0006987': 'GO:0036498', 'GO:0001007': 'GO:0006359', 'GO:0004716': 'GO:0023014',
               'GO:0001191': 'GO:0001085', 'GO:0001076': 'GO:0001085', 'GO:0000989': 'GO:0008134',
               'GO:0001190': 'GO:0001085', 'GO:0000991': 'GO:0140110', 'GO:0001129': 'GO:0051123',
               'GO:2001275': 'GO:1900078', 'GO:0031659': 'GO:0045737'} 

# Produce output text file of GO descriptions
with open(args.inputTextFile, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                line = line.rstrip('\r\n')
                # Handle blank lines
                if line == '.':
                        fileOut.write('.\n')
                # Handle GO code lines
                else:
                        # Get the names of GO entries
                        codes = line.split('; ')
                        names = []
                        for code in codes:
                                names.append(go[code].name)
                        # Format and write to file
                        fileOut.write('; '.join(names) + '\n')

# Done!
print('Program completed successfully!')
