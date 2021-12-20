#! python3
# domtblout2go

# Python script to receive a parsed domtblout file and, alongside
# the pfam2go data file, create gene ontology annotations for gene
# models. It will also expand the terms to include parent terms
# if a go-basic.obo file is provided.

import argparse, os
from goatools import obo_parser

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
        if args.oboFile != None:
            if not os.path.isfile(args.oboFile):
                print('I am unable to locate the GO .obo file (' + args.oboFile + ')')
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

def parse_pfam2go(pfam2goFile, oboFile):
        # Optional processing of GO parent expansion
        if oboFile != None:
                go = obo_parser.GODag(oboFile)
                # Update annotations file
                obsoletedGOs = ['GO:0055114'] # go-basic.obo file was downloaded 16-12-21
                replacedGOs = {'GO:0140603': 'GO:0016887', 'GO:0036425': 'GO:0036424', 'GO:0005671': 'GO:0140671',
                               'GO:2000574': 'GO:0140659', 'GO:0102132': 'GO:0004316', 'GO:0102131': 'GO:0004316',
                               'GO:0009405': 'GO:0052031', 'GO:0015002': 'GO:0016491', 'GO:0052331': 'GO:0044179',
                               'GO:0000186': 'GO:0000165', 'GO:2000575': 'GO:0140661'} # Modifications were made 20-12-21
        
        # Normal processing
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
                        # Optional expansion to include parent terms
                        if oboFile != None:
                                # Handle GO obsoletions
                                skipThis = False
                                for entry in obsoletedGOs:
                                        if entry == goAccession:
                                                print('Deleted ' + entry)
                                                skipThis = True
                                if skipThis == True:
                                        continue # completely skip this line in pfam2go file
                                # Handle GO replacements
                                for key, value in replacedGOs.items():
                                        if key == goAccession:
                                                print('Replaced ' + key)
                                                goAccession = value
                                # Populate ancestors of GO terms
                                if goAccession not in go:
                                        print('GO term needs replacement/obsoletion! == ' + goAccession)
                                else:
                                        goAccession = list(go[goAccession].get_all_parents())
                        # Store details in dictionary
                        if pfamAbbrev not in pfam2goDict:
                                if type(goAccession).__name__ == "list":
                                        pfam2goDict[pfamAbbrev] = goAccession
                                else:
                                        pfam2goDict[pfamAbbrev] = [goAccession]
                        else:
                                if type(goAccession).__name__ == "list":
                                        pfam2goDict[pfamAbbrev] += goAccession
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
p.add_argument("-g", "-goObo", dest="oboFile",
                   help="Optionally input a go-basic.obo file to expand GOs to include parent terms (not required).")
p.add_argument("-o", "-output", dest="outputFileName", required=True,
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Read in the parsed hmmer domblout file
domDict = read_parsed_domtblout(args.parsedDomtblout)

# Read in the pfam2go file
pfam2goDict = parse_pfam2go(args.pfam2goFile, args.oboFile)

# Associate domDict and pfam2goDict
associatedDict = associate_parsed_domtblout_and_pfam2go(domDict, pfam2goDict)

# Generate output
with open(args.outputFileName, 'w') as fileOut:
        for sequenceID, goTerms in associatedDict.items():
                fileOut.write("{0}\t{1}\n".format(sequenceID, goTerms))

# All done!
print('Program completed successfully!')
