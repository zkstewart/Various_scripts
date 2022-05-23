#! python3
# get_cafe_family_ids_containing_sequences.py
# Script to pull out the IDs of families that
# expanded/contracted which contain one or more of any
# sequences identified in one or more FASTA input files

import os, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.cafeResultsFileName):
        print('I am unable to locate the Base_family_results.txt file (' + args.cafeResultsFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.orthogroupsFileName):
        print('I am unable to locate the Orthogroups.tsv file (' + args.orthogroupsFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for value in args.sequenceFilesOrDirs:
        if not os.path.isdir(value) and not os.path.isfile(value):
            print('The specified directory or file does not exist (' + value + ')')
            print('Make sure you\'ve typed the location correctly and try again.')
            quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location ("{0}")'.format(args.outputFileName))
        print('Make sure you specify a unique file name and try again.')
        quit()

def get_significant_cafe_families(baseFamilyResultsFileName):
    '''
    Parses the Base_family_results.txt file, produced by CAFE 5, and returns
    a list of all the families that are detected as significantly expanded/contracted
    in one or more lineage to a P-value threshold of 0.05.
    
    Parameters:
        baseFamilyResultsFileName -- a string indicating the location of the Base_family_results.txt file.
    Returns:
        significantFamilies -- a list containing strings indicating the gene families (/orthogroups) that
                               expanded or contracted to a significant P-value (==0.05).
    '''
    significantFamilies = []
    with open(baseFamilyResultsFileName, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            sl = line.rstrip("\r\n").split("\t")
            
            # Hold onto familiy ID if it's significant
            if sl[2] == "y": # "y" is 0.05 significant, otherwise it's "n" in CAFE 5
                significantFamilies.append(sl[0])
    return significantFamilies

def get_orthogroup_tsv_sequences(tsvFileName):
    '''
    Parses the Orthogroups.tsv file, produced by OrthoFinder, and returns a dictionary
    relating orthogroup IDs to the sequences that are part of that familie. There is no
    segregation of different species' IDs in this; it is just a flat list.
    
    Parameters:
        tsvFileName -- a string indicating the location of the Orthogroups.tsv file.
    Returns:
        orthoSeqsDict -- a dictionary with structure:
                                {
                                    orthogroupA: [
                                        seqID_1, seqID_2, ...
                                        ],
                                    orthogroupB: [ ... ],
                                    ...
                                }
    '''
    orthoSeqsDict = {}
    with open(tsvFileName, "r") as fileIn:
        for line in fileIn:
            if line.startswith("Orthogroup"):
                continue
            sl = line.rstrip("\r\n").split("\t")
            
            # Get a flat list of tab-separated column values
            orthoSeqsDict[sl[0]] = [seqID for famGroup in sl[1:] for seqID in famGroup.split(", ") if seqID != ""] 
    return orthoSeqsDict

def hierarchical_parse_of_fasta_files(fastaFileNamesOrDirs):
    '''
    Parameters:
        fastaFileNamesOrDirs -- a list containing strings pointing to an individual FASTA
                                file, or to a directory that contains FASTA files. Subdirectories
                                will be ignored i.e., it does not parse recursively.
    Returns:
        hierarchicalFastaIDsDict -- a dictionary with structure like:
                                        {
                                            fileStringA: [
                                                seqID_1, seqID_2, ...
                                            ],
                                            dirString: {
                                                fileStringB: [
                                                    seqID_3, seqID_4, ...
                                                ],
                                                fileStringC: [
                                                    ...
                                                ], ...
                                            },
                                            ...
                                        }
    '''
    hierarchicalFastaIDsDict = {}
    for value in fastaFileNamesOrDirs:
        # Handle directory or FASTA file
        if os.path.isdir(value):
            fastaFiles = [os.path.join(value, subFasta) for subFasta in os.listdir(value) if os.path.isfile(os.path.join(value, subFasta))]
            hierarchicalFastaIDsDict[os.path.basename(value)] = {} # this will give us the hierarchical structure to our dictionary
            isDir = True
        elif os.path.isfile(value):
            fastaFiles = [value]
            fileBase = os.path.basename(value.rsplit(".", maxsplit=1)[0]) # remove the file .suffix
            hierarchicalFastaIDsDict[fileBase] = [] # this will be a flat result associated to just the parent file name key (sans suffix)
            isDir = False
        else:
            raise Exception("Hierarchical FASTA parse failed because '{0}' is neither file nor directory".format(value))
        
        # Parse any FASTA file(s) found
        for file in fastaFiles:
            fileBase = os.path.basename(file.rsplit(".", maxsplit=1)[0]) # remove the file .suffix (this is the same as the above fileBase declaration if .isfile() == True)
            with open(file, "r") as fileIn:
                    for id, seq in SimpleFastaParser(fileIn):
                        if isDir:
                            hierarchicalFastaIDsDict[os.path.basename(value)].setdefault(fileBase, [])
                            hierarchicalFastaIDsDict[os.path.basename(value)][fileBase].append(id)
                        else:
                            hierarchicalFastaIDsDict[fileBase].append(id)
    return hierarchicalFastaIDsDict

def get_families_containing_a_sequence(orthoSeqsDict, familiesList, hierarchicalFastaIDsDict):
    '''
    Parameters:
        orthoSeqsDict -- a dictionary with the below structure, where the seqID values
                         should map to values in hierarchicalFastaIDsDict:
                                {
                                    orthogroupA: [
                                        seqID_1, seqID_2, ...
                                        ],
                                    orthogroupB: [ ... ],
                                    ...
                                }
        familiesList -- a list containing strings representing families identified
                        as significant for expansion/contraction by CAFE. Strings
                        should map to keys in orthoSeqsDict.
        hierarchicalFastaIDsDict -- a dictionary with structure like:
                                        {
                                            fileStringA: [
                                                seqID_1, seqID_2, ...
                                            ],
                                            dirString: {
                                                fileStringB: [
                                                    seqID_3, seqID_4, ...
                                                ],
                                                fileStringC: [
                                                    ...
                                                ], ...
                                            },
                                            ...
                                        }
    Returns:
        familiesWithSequence -- a list containing string values indicating the
                                families that contain a FASTA sequence, and where
                                that sequence comes from in the hierarchy of
                                hierarchicalFastaIDsDict i.e., something like
                                [
                                    "orthogroupA (\t) fileStringA, fileStringB, ...",
                                    "orthogroupB (\t) dirString (\t) fileStringC, ...",
                                    ...
                                ]
    '''
    hierarchicalResults = {} # we'll store results here temporarily
    for parentKey, parentValue in hierarchicalFastaIDsDict.items():
        # Split into separate code path to handle file and directory hierarchies
        if type(parentValue).__name__ == "dict":
            # Dig into the next level of the dictionary hierarchy
            for subkey, idsList in parentValue.items():
                # See if we can relate any sequence IDs from the three input parameters
                for seqID in idsList:
                    found = None
                    for og, ogSeqsList in orthoSeqsDict.items():
                        if seqID in ogSeqsList and og in familiesList:
                            found = og
                            break
                    if found != None:
                        hierarchicalResults.setdefault(og, {})
                        hierarchicalResults[og].setdefault(parentKey, [])
                        if subkey not in hierarchicalResults[og][parentKey]:
                            hierarchicalResults[og][parentKey].append(subkey)
        else:
            # See if we can relate any sequence IDs from the three input parameters
            for seqID in parentValue:
                found = None
                for og, ogSeqsList in orthoSeqsDict.items():
                    if seqID in ogSeqsList and og in familiesList:
                        found = og
                        break
                if found != None:
                    hierarchicalResults.setdefault(og, {})
                    hierarchicalResults[og].setdefault(parentKey, None) # none since we don't have hierarchy
            
    # Provide a flattened, human readable output list
    familiesWithSequence = []
    for family in familiesList:
        if family in hierarchicalResults:
            result = hierarchicalResults[family]
            # Coerce results into more sensible outputs
            individualResults = []
            directoryResults = []
            for key, value in result.items():
                # Handle individual FASTA file results
                if value == None:
                    individualResults.append(key)
                # Handle directory-level results
                else:
                    directoryResults.append("{0}\t{1}".format(key, ", ".join(value)))
            individualResults = ", ".join(individualResults)
            # Store outputs
            if directoryResults != []:
                for dResult in directoryResults:
                    familiesWithSequence.append("{0}\t{1}".format(family, dResult))
            if individualResults != "":
                familiesWithSequence.append("{0}\t{1}".format(family, individualResults))
        else:
            pass # do nothing
    return familiesWithSequence

def main():
    # User input
    usage = """%(prog)s reads in the Base_family_results.txt file produced by CAFE,
    the Orthogroups.tsv file produced by OrthoFinder, and one or more files / directories
    containing FASTA files. With this, it will output any expanded/contracted families
    that contained one or more of the sequences from within the FASTA file(s), noting the
    FASTA file name and the sequence itself.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-r", dest="cafeResultsFileName",
        help="Input Base_family_results.txt file name and location")
    p.add_argument("-g", dest="orthogroupsFileName",
        help="Input Orthogroups.tsv file name and location")
    p.add_argument("-s", dest="sequenceFilesOrDirs", nargs="+",
        help="Specify one or more files or directories that ONLY contain FASTA files")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name where tabular output will be written")
    args = p.parse_args()
    validate_args(args)
    
    # Parse CAFE results file
    significantFamilies = get_significant_cafe_families(args.cafeResultsFileName)
    
    # Parse OrthoFinder TSV file
    orthoSeqsDict = get_orthogroup_tsv_sequences(args.orthogroupsFileName)
    
    # Parse all locateable FASTA files
    hierarchicalFastaIDsDict = hierarchical_parse_of_fasta_files(args.sequenceFilesOrDirs)
    
    # Locate CAFE families that contain a sequence identified in a parsed FASTA
    familiesWithSequence = get_families_containing_a_sequence(orthoSeqsDict, significantFamilies, hierarchicalFastaIDsDict)
    
    # Produce output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(familiesWithSequence))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
