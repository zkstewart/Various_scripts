#! python3
# reconstitute_domains.py
# Script to assess whether a set of domain HMMs can be used to reconstitute the groupings of
# sequences used to create the HMMs.

import argparse, os, re
from Bio import SeqIO

# Various functions for program operations

## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.domtbloutParse):
                print('I am unable to locate the domtblout parse file (' + args.domtbloutParse + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate the input file locations depending on type of input
        if not os.path.isdir(args.referenceFastaDir):
                print('The provided reference fasta dir value "' + args.referenceFastaDir + '" does not point to a directory.')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
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

def parse_parsed_domtblout(fileLocation):
        sequenceDict = {}
        domainDict = {}
        with open(fileLocation, "r") as fileIn:
                for line in fileIn:
                        # Clean data
                        sl = line.rstrip("\r\n").split("\t")
                        seqID = sl[0]
                        domains = []
                        for i in range(1, len(sl)):
                                domains.append(eval(sl[i])) # Data is stored in a list-like format already within parsed domtblout
                        # Associate to sequenceDict
                        sequenceDict[seqID] = domains
                        # Associate to domainDict
                        for domain in domains:
                                if domain[0] not in domainDict:
                                        domainDict[domain[0]] = [seqID]
                                else:
                                        domainDict[domain[0]].append(seqID)
        return sequenceDict, domainDict

def locate_files_in_directory(directory):
        files = []
        for item in os.listdir(directory):
                item_path = os.path.abspath(os.path.join(directory, item))
                if os.path.isfile(item_path):
                        files.append(item_path)
        return files

def sequence_ownership_from_multiple_fastas(fastaFileList):
        sequenceOwnershipDict = {}
        for fastaFile in fastaFileList:
                with open(fastaFile, "r") as fileIn:
                        fastaFileBasename = os.path.basename(fastaFile)
                        records = SeqIO.parse(fileIn, "fasta")
                        sequenceOwnershipDict[fastaFileBasename] = {"ids": [], "descriptions": []}
                        for record in records:
                                sequenceOwnershipDict[fastaFileBasename]["ids"].append(record.id)
                                sequenceOwnershipDict[fastaFileBasename]["descriptions"].append(record.description)
        return sequenceOwnershipDict

def domaindict_sequenceownershipdict_tabulate(domainDict, sequenceOwnershipDict):
        # Get a clean list of keys
        domainKeys = []
        domainCleanPairs = {}
        for key in domainDict.keys():
                if key.startswith(".\\"):
                        cleanKey = key[2:]
                domainKeys.append(cleanKey)
                domainCleanPairs[cleanKey] = key
        
        ownershipKeys = []
        ownershipCleanPairs = {}
        for key in sequenceOwnershipDict.keys():
                suffix = key.rsplit(".", maxsplit=1)[-1]
                if suffix.startswith("fa"):
                        cleanKey = key.rsplit(".", maxsplit=1)[0]
                ownershipKeys.append(cleanKey)
                ownershipCleanPairs[cleanKey] = key
        
        # Validate that the two dicts pair with each other
        assert len(domainKeys) == len(ownershipKeys)
        for key in domainKeys:
                assert key in ownershipKeys
        
        # Tabulate an output to facilitate comparison of groups
        tableGroups = []
        for key in domainKeys:
                # Re-obtain the pre-cleaned key
                domainKey = domainCleanPairs[key]
                ownershipKey = ownershipCleanPairs[key]
                # Get the values for each group
                domainGroup = domainDict[domainKey]
                ownershipGroupId = sequenceOwnershipDict[ownershipKey]["ids"]
                ownershipGroupDescription = sequenceOwnershipDict[ownershipKey]["descriptions"]
                # Clean the values for each group
                '''
                This is only necessary since the files I'm working with currently
                have issues with their sequence IDs
                '''
                pipeCdsRegex = re.compile(r"\|CDS\d{1,3}.+$")
                for i in range(len(ownershipGroupId)):
                        ownershipGroupId[i] = ownershipGroupId[i].rsplit("|", maxsplit=1)[0]
                        suffixToRemove = pipeCdsRegex.findall(ownershipGroupId[i])
                        if suffixToRemove != []:
                                ownershipGroupId[i] = ownershipGroupId[i][:len(suffixToRemove[0])]
                for i in range(len(ownershipGroupDescription)):
                        ownershipGroupDescription[i] = ownershipGroupDescription[i].rsplit("|", maxsplit=1)[0]
                        suffixToRemove = pipeCdsRegex.findall(ownershipGroupDescription[i])
                        if suffixToRemove != []:
                                ownershipGroupDescription[i] = ownershipGroupDescription[i][:len(suffixToRemove[0])]
                for i in range(len(domainGroup)):
                        domainGroup[i] = domainGroup[i].rsplit("|", maxsplit=1)[0]
                        suffixToRemove = pipeCdsRegex.findall(domainGroup[i])
                        if suffixToRemove != []:
                                domainGroup[i] = domainGroup[i][:len(suffixToRemove[0])]
                ownershipGroupId = set(ownershipGroupId)
                ownershipGroupDescription = set(ownershipGroupDescription)
                domainGroup = set(domainGroup)
                # Pair our ownership values up with domain values
                if ownershipGroupId == ownershipGroupDescription:
                        ownershipGroup = ownershipGroupId
                else:
                        idsIntersection = domainGroup.intersection(ownershipGroupId)
                        descriptionsIntersection = domainGroup.intersection(ownershipGroupDescription)
                        if len(idsIntersection) > len(descriptionsIntersection):
                                ownershipGroup = ownershipGroupId
                        else:
                                ownershipGroup = ownershipGroupDescription
                assert len(ownershipGroup) != 0
                # Compare groups to find values in common and values unique to either group
                common = list(domainGroup.intersection(ownershipGroup))
                uniqueToDomain = list(domainGroup.difference(ownershipGroup))
                uniqueToOwnership = list(ownershipGroup.difference(domainGroup))
                common.sort()
                uniqueToDomain.sort()
                uniqueToOwnership.sort()
                # Tabulate for the table grouping
                tableGroups.append([key, common + ["... unique ..."] + uniqueToDomain, common + ["... unique ..."] + uniqueToOwnership])
        



# Main call
def main():
        #### USER INPUT SECTION
        usage = """Wrapper script to perform MMseqs2 search. Provide the arguments below."""
        # Required
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-dp", "--domtblout_parse", dest="domtbloutParse", type = str,
                          help="Specify the parsed domtblout file")
        p.add_argument("-r", "--reference_fasta_dir", dest="referenceFastaDir", type = str,
                          help="Specify the directory containing FASTA files used to generate HMMs")
        p.add_argument("-o", "--output_file", dest="outputFileName", type = str,
                          help="Specify the output file name")
        
        args = p.parse_args()
        ## HARD-CODED TESTING
        args.domtbloutParse = r"F:\toxins_annot\analysis\reconstitute_test_database\refs_toxins_domains.domtblout_parse"
        domtbloutParse = r"F:\toxins_annot\analysis\reconstitute_test_database\refs_toxin_domains.domtblout_parse"
        args.referenceFastaDir = r"F:\toxins_annot\alns\Manual_domain_curation_MSAs"
        referenceFastaDir = r"F:\toxins_annot\alns\Manual_domain_curation_MSAs"
        args.outputFileName = r"F:\toxins_annot\analysis\reconstitute_test_database\test_reconstitute.txt"
        outputFileName = r"F:\toxins_annot\analysis\reconstitute_test_database\test_reconstitute.txt"
        #validate_args(args)
        ## END

        # Parse parsed domtblout
        #sequenceDict, domainDict = parse_parsed_domtblout(args.domtbloutParse)
        sequenceDict, domainDict = parse_parsed_domtblout(domtbloutParse)

        # Parse reference directory FASTAs
        referenceFastas = locate_files_in_directory(referenceFastaDir)
        sequenceOwnershipDict = sequence_ownership_from_multiple_fastas(referenceFastas)

        # Format group comparison output file

        # Done!
        print('Program completed successfully!')

if __name__ == '__main__':
        main()
