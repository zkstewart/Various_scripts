#! python3
# domtblout2go

# Python script to receive a parsed domtblout file and, alongside
# the pfam2go data file, create gene ontology annotations for gene
# models.

import argparse, os

# Define functions for later use

## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.parsedDomtblout):
                print('I am unable to locate the parsed domtblout file (' + args.parsedDomtblout + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.pfam2goFile):
                print('I am unable to locate the pfam2go file (' + args.pfam2goFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

## Input parsing
def read_parsed_domtblout(parsedDomtblout):
        domDict = {}
        with open(parsedDomtblout, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('#') or line == '' or line == '\n' or line == '\r\n':
                                continue
                        # Parse line and extract domains associated with each sequence
                        sl = line.rstrip('\r\n').split()
                        pid = sl[0]
                        domains = []
                        for i in range(1, len(sl[1:])):
                            d = sl[i].split(",")[0].strip("'[]")
                            domains.append(d)
                        # Store details in dictionary
                        domDict[pid] = domains
        return domDict

def parse_pfam2go(pfam2goFile):
        pfam2goDict = {}
        with open(pfam2goFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('!') or line == '' or line == '\n' or line == '\r\n':
                                continue
                        # Parse line and extract relevant details
                        line = line.rstrip("\r\n")
                        pfamAbbrev = line.split(" > ")[0].split(" ")[1]
                        goAccession = line.split(" ; ")[1]
                        # Store details in dictionary
                        if pfamAbbrev not in pfam2goDict:
                                pfam2goDict[pfamAbbrev] = [goAccession]
                        else:
                                pfam2goDict[pfamAbbrev].append(goAccession)
        return pfam2goDict

## Association
def associate_parsed_domtblout_and_pfam2go(domDict, pfam2goDict):
    associatedDict = {}
    for sequenceID, domains in domDict.items():
        # Perform the association
        goTerms = []
        for domain in domains:
            if domain in pfam2goDict:
                goTerms += pfam2goDict[domain]
        goTerms = list(set(goTerms)) # remove redundancy
        # Format for text output
        if goTerms == []:
            goTerms = "."
        else:
            goTerms = "; ".join(goTerms)
        # Store in output dictionary
        associatedDict[sequenceID] = goTerms
    return associatedDict

#### USER INPUT SECTION
usage = """%(prog)s reads a parsed HMMER domtblout file as well as the pfam2go annotations
data file. It produces gene ontology annotations for each gene model that are associated with
the gene's pfam domain annotations.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="parsedDomtblout", required=True,
                   help="Input domtblout HMMER3 domtblout file.")
p.add_argument("-p", "-pfam2goFile", dest="pfam2goFile", required=True,
                   help="Input pfam2go file.")
p.add_argument("-o", "-output", dest="outputFileName", required=True,
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Read in the parsed hmmer domblout file
domDict = read_parsed_domtblout(args.parsedDomtblout)

# Read in the pfam2go file
pfam2goDict = parse_pfam2go(args.pfam2goFile)

# Associate domDict and pfam2goDict
associatedDict = associate_parsed_domtblout_and_pfam2go(domDict, pfam2goDict)

# Generate output
with open(args.outputFileName, 'w') as fileOut:
        for sequenceID, goTerms in associatedDict.items():
                fileOut.write("{0}\t{1}\n".format(sequenceID, goTerms))

# All done!
print('Program completed successfully!')
