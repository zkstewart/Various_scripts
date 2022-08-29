#! python3
# prep_qtlseq.py
# Script to make it easier to format commands to run QTLseq
# for bulk segregant analysis

import os, argparse, re

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isdir(args.bamDirectory):
        print('I am unable to locate the BAM directory (' + args.bamDirectory + ')')
        print('Make sure you\'ve typed the location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print('I am unable to locate the metadata file (' + args.metadataFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.referenceFastaFile):
        print('I am unable to locate the reference FASTA file (' + args.referenceFastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def parse_metadata_file(metadataFile, sampleIdCol, bulkCol):
    '''
    This function will parse a metadata file for QTL-seq or MutMap file,
    as a TSV or CSV, and produce a dictionary associating sample IDs
    to their bulk segregant analysis grouping.
    
    Parameters:
        metadataFile -- a string indicating the file location of the metadata
                        with TSV or CSV formatting.
        sampleIdCol -- a string indicating the header value for the column containing
                       sample identifiers.
        bulkCol -- a string indicating the header value for the column containing
                   bulk segregant groups for comparison.
    Returns:
        sampleDict -- a dictionary with structure like:
                      {
                          'sampleid1': 'p1',
                          'sampleid2': 'p1_c1',
                          ...
                      }
    '''
    header = None
    splitChar = None
    sampleDict = {}
    
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            # Check if parsing CSV or TSV
            if splitChar is None:
                if "\t" in line:
                    splitChar = "\t"
                else:
                    splitChar = ","
            
            # Handle header line
            sl = line.rstrip("\r\n ").split(splitChar)
            if header == None:
                header = sl
                
                # Raise errors if columns don't exist
                if sampleIdCol not in header:
                    raise Exception("Sample ID column does not exist")
                if bulkCol not in header:
                    raise Exception("Genotype column(s) do not exist")
                
                # Get column indices
                sampleIndex = header.index(sampleIdCol)
                bulkIndex = header.index(bulkCol)
            # Handle content lines
            else:
                sampleID = sl[sampleIndex]
                bulk = sl[bulkIndex]
                assert sampleID not in sampleDict, "Sample ID '{0}' isn't unique!".format(sampleDict)
                
                sampleDict[sampleID] = bulk
                
    return sampleDict

def process_qtlseq_comparisons(sampleDict, blankChar="."):
    '''
    Receives the output sampleDict from parse_metadata_file() and
    parses it so we can have a new dictionary indicating the bulks
    for comparison.
    
    Parameters:
        sampleDict -- a dictionary with structure like:
                      {
                          'sampleid1': 'p1',
                          'sampleid2': 'p1_c1',
                          ...
                      }
    Returns:
        comparisonDict -- a dictionary with structure like:
                          {
                              comparisonNum:
                              {
                                  'parent': parentSampleID,
                                  'bulk1': [bulk1SampleID1, ...],
                                  'bulk2': [bulk2SampleID1, ...]
                              },
                              ...
                          }
    '''
    QTLSEQ_BULK_FORMAT = re.compile(r"p\d{1,2}(_c[12])?$")
    
    # Check and valiate the number of comparisons to perform
    comparisons = set()
    for value in sampleDict.values():
        if value == blankChar:
            continue
        assert QTLSEQ_BULK_FORMAT.search(value) is not None, \
            f"'{value}' is not properly formatted"
        
        parent = int(value.split("_")[0][1:]) # remove _ suffix, remove the p prefix
        comparisons.add(parent)
    comparisons = list(comparisons)
    comparisons.sort()
    assert 1 in comparisons, \
        "Comparison groupings must count from 1"
    assert all([comparisons[i] == comparisons[i-1]+1 for i in range(1, len(comparisons))]), \
        "Comparison groups must be sequential (no skipped numbers)"
    
    # Parse comparisons into dictionary format
    comparisonDict = {num: {"parent": None, "bulk1": [], "bulk2": []} for num in comparisons}
    for num in comparisons:
        for sampleid, bulkGroup in sampleDict.items():
            if f"p{num}" in bulkGroup:
                # Handle parents
                if f"p{num}" == bulkGroup:
                    assert comparisonDict[num]["parent"] is None, \
                        f"More than one parent not allowed for comparison {num}"
                    comparisonDict[num]["parent"] = sampleid
                # Handle bulks
                else:
                    if "c1" in bulkGroup:
                        comparisonDict[num]["bulk1"].append(sampleid)
                    else:
                        comparisonDict[num]["bulk2"].append(sampleid)
    
    # Validate that all comparisons have a parent and one or more bulked samples
    for num in comparisons:
        assert comparisonDict[num]["parent"] is not None, \
            f"Comparison '{num}' lacks a parent"
        assert comparisonDict[num]["bulk1"] is not [], \
            f"Comparison '{num}' lacks samples in bulk 1 group"
        assert comparisonDict[num]["bulk2"] is not [], \
            f"Comparison '{num}' lacks samples in bulk 2 group"
    
    return comparisonDict

def get_comparisonDict_filenames(comparisonDict, bamDirectory, bamSuffix):
    '''
    Receives a comparisonDict as created by process_qtlseq_comparisons() or
    process_mutmap_comparisons() and replaces the sampleID values with
    paths to the relevant BAM files, after having validated that they exist.
    
    This modification occurs in-place, so the original dictionary will change!
    
    Parameters:
        comparisonDict -- a dictionary with structure like:
                          {
                              comparisonNum:
                              {
                                  'parent': parentSampleID,
                                  'bulk1': [bulk1SampleID1, ...],
                                  'bulk2': [bulk2SampleID1, ...]
                              },
                              ...
                          }
        bamDirectory -- a string indicating the location of the BAM files.
        bamSuffix -- a string indicating the suffix following sample ID values
                     which creates the appropriate BAM file name.
    '''
    for groupsDict in comparisonDict.values():
        for key, value in groupsDict.items():
            # Handle parent group
            if key == "parent":
                fileName = os.path.join(bamDirectory, f"{value}{bamSuffix}")
                assert os.path.isfile(fileName), \
                    f"'{fileName}' does not exist or cannot be found at the path provided"
                groupsDict[key] = fileName
            
            # Handle bulk groups
            else:
                for i in range(len(value)):
                    fileName = os.path.join(bamDirectory, f"{value[i]}{bamSuffix}")
                    assert os.path.isfile(fileName), \
                        f"'{fileName}' does not exist or cannot be found at the path provided"
                    groupsDict[key][i] = fileName

def create_qtlseq_shell_script(comparisonDict, referenceFastaFile, outputFileName):
    '''
    Receives a comparisonDict as modified by get_comparisonDict_filenames()
    and formats a shell script that can be qsub-ed to run QTL-seq on the
    appropriate bulk groupings.
    
    Parameters:
        comparisonDict -- a dictionary with structure like:
                          {
                              comparisonNum:
                              {
                                  'parent': parentSampleBAM,
                                  'bulk1': [bulk1SampleBAM1, ...],
                                  'bulk2': [bulk2SampleBAM1, ...]
                              },
                              ...
                          }
        referenceFastaFile -- a string indicating the location of the reference FASTA
                              file that mapping was performed against.
        outputFileName -- a string indicating the location to write the shell script to.
    '''
    # Specify hard-coded script features
    JOBNAME = "qtlseq"
    WALLTIME = "120:00:00"
    MEM = "30G"
    
    # Get qtlseq commands
    cmds = []
    for comparisonNum, groupsDict in comparisonDict.items():
        # Derive relevant parameters from comparisonDict
        n1 = len(groupsDict["bulk1"])
        n2 = len(groupsDict["bulk2"])
        
        # Format and store command
        cmd = """qtlseq -r {0} \\
    -p {1} \\
{2}
{3}
    -n1 {4} \\
    -n2 {5} \\
    -o qtlseq_output_comparison{6}
""".format(
    referenceFastaFile, groupsDict["parent"],
    "\n".join(["    -b1 {0} \\".format(b) for b in groupsDict["bulk1"]]),
    "\n".join(["    -b2 {0} \\".format(b) for b in groupsDict["bulk2"]]),
    n1, n2, comparisonNum
    )
        cmds.append(cmd)
    cmds = "".join(cmds)
    
    # Setup the script's contents
    scriptText = f"""#!/bin/bash -l
#PBS -N {JOBNAME}
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus=1\n
cd ${{PBS_O_WORKDIR}}\n
{cmds}"""

    # Write to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def main():
    # User input
    usage = """%(prog)s enables easier running of QTL-seq by receiving a metadata file
    containing rows and columns which will facilitate the creation of a qsub script to
    run QTL-seq with.
    
    The metadata file is expected to contain two columns; one indicating the prefixes
    to sample BAM files, and another indicating one or more comparisons to perform.
    A single parent sample (p#) should be specified, with two child bulks (p#_c1, p#_c2)
    specified. The # should be 1 for the first comparison, 2 for the next, etc.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="bamDirectory",
                   required=True,
                   help="Input directory where BAM files reside")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Input Metadata TSV (or CSV) file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the file name to write qsub script to")
    p.add_argument("-r", dest="referenceFastaFile",
                   required=True,
                   help="Specify the file name of the reference FASTA that mapping was done with")
    # Opts
    p.add_argument("--sampleid", dest="sampleIdCol",
                   required=False,
                   help="""Column name where sample ID is located; this should provide the
                   prefix to the BAM files""",
                   default="sampleid")
    p.add_argument("--qtlseq", dest="qtlseqCol",
                   required=False,
                   help="Column name where QTL-seq comparison details can be found",
                   default="qtlseq")
    p.add_argument("--bamSuffix", dest="bamSuffix",
                   required=False,
                   help="Suffix which follows every sampleid value to indicate a BAM file",
                   default=".sorted.bam")
    args = p.parse_args()
    validate_args(args)
    
    # Parse metadata file columns
    sampleDict = parse_metadata_file(args.metadataFile, args.sampleIdCol, args.qtlseqCol)
    
    # Process for QTL-seq comparisons
    comparisonDict = process_qtlseq_comparisons(sampleDict)
    
    # Process BAM file names and validate that they exist
    get_comparisonDict_filenames(comparisonDict, args.bamDirectory, args.bamSuffix)
    
    # Format and write QTL-seq script
    create_qtlseq_shell_script(comparisonDict, args.referenceFastaFile, args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
