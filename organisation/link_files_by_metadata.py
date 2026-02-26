#!/usr/bin/env python3
# link_files_by_metadata.py
# Script to make easy the process of symlinking files for the purpose of
# renaming them, according to the contents of a metadata file. The script
# handles the logic of identifying paired/single files and validating that
# files which should exist, do actually exist.

import os, argparse
import pandas as pd

def validate_args(args):
    # Validate input location(s)
    args.inputLocations = [ os.path.abspath(x) for x in args.inputLocations ]
    for location in args.inputLocations:
        if not os.path.exists(location):
            raise FileNotFoundError(f"File or directory given to -i ('{location}') does not exist!")
    
    args.metadataFile = os.path.abspath(args.metadataFile)
    if not os.path.isfile(args.metadataFile):
        raise FileNotFoundError(f"Value given to -m ('{args.metadataFile}') is not a file!")
    
    # Validate that argument values are sensible
    if args.fileSuffix == "":
        raise ValueError(f"-s value cannot be blank!")
    if args.newSuffix == None:
        args.newSuffix = args.fileSuffix
    
    # Validate output location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        parentDir = os.path.dirname(args.outputDirectory)
        if os.path.isdir(parentDir):
            os.makedirs(args.outputDirectory, exist_ok=True)
            print(f"# Created output location '{args.outputDirectory}' as part of argument validation")
        else:
            raise ValueError(f"-o location '{args.outputDirectory}' does not exist and its parent dir also does not exist; "
                             f"this script won't create the -o location under these circumstances; " +
                             "please fix this issue manually before trying again.")

def handle_file_name(file, fileSuffix, found):
    '''
    Helper function for locate_files() to set values in its
    found dictionary.
    '''
    filePrefix = os.path.basename(file).rsplit(fileSuffix, maxsplit=1)[0].rstrip(". ")
    
    if filePrefix.endswith("_1") or filePrefix.endswith("_2"):
        pairedPrefix = filePrefix.rsplit("_", maxsplit=1)[0]
        found.setdefault(pairedPrefix, {1: None, 2: None})
    
    if filePrefix.endswith("_1"):
        if not isinstance(found[pairedPrefix], dict):
            raise KeyError(f"File '{file}' ends with paired sequencing suffix ('_1') but has a prefix that already " +
                           f"matches a seemingly unpaired file i.e., '{found[pairedPrefix]}'; duplicates cannot be tolerated")
        
        if found[pairedPrefix][1] != None:
            raise KeyError(f"File '{file}' appears to be the forward ('_1') read for the prefix '{pairedPrefix}' " + 
                           f"but another matching file has already been found previously i.e., '{found[pairedPrefix][1]}'; " +
                           "duplicates cannot be tolerated")
        found[pairedPrefix][1] = file
    elif filePrefix.endswith("_2"):
        if not isinstance(found[pairedPrefix], dict):
            raise KeyError(f"File '{file}' ends with paired sequencing suffix ('_2') but has a prefix that already " +
                           f"matches a seemingly unpaired file i.e., '{found[pairedPrefix]}'; duplicates cannot be tolerated")
        
        if found[pairedPrefix][2] != None:
            raise KeyError(f"File '{file}' appears to be the reverse ('_2') read for the prefix '{pairedPrefix}' " + 
                           f"but another matching file has already been found previously i.e., '{found[pairedPrefix][2]}'; " +
                           "duplicates cannot be tolerated")
        found[pairedPrefix][2] = file
    else:
        if filePrefix in found:
            raise KeyError(f"File '{file}' lacks an expected paired sequencing suffix and so I expect it to occur as " +
                           f"a single file. However, the prefix of this file ('{filePrefix}') was already found previously; " +
                           "duplicates cannot be tolerated")
        
        found[filePrefix] = file

def locate_files(inputLocations, fileSuffix):
    '''
    Returns:
        found -- a dictionary with structure like:
                 {
                     'origFilePrefix1': 'singleEndFile1',
                     'origFilePrefix2': { 1: 'pairedForwardFile1', 2: 'pairedReverseFile1' },
                     ...
                 }
    '''
    found = {}
    
    # Locate files
    for location in inputLocations:
        location = os.path.abspath(location)
        if os.path.isfile(location):
            if location.endswith(fileSuffix):
                handle_file_name(location, fileSuffix, found)
        elif os.path.isdir(location):
            for file in os.listdir(location):
                file = os.path.join(location, file)
                if file.endswith(fileSuffix):
                    handle_file_name(file, fileSuffix, found)
        else:
            raise TypeError(f"Unable to verify that '{location}' was a file or a directory; " +
                            "make sure this location exists and/or is entered correctly")
    
    # Run post-processing checks for data validity
    for filePrefix, value in found.items():
        if isinstance(value, dict):
            if value[1] == None or value[2] == None:
                problemIndex, goodIndex = (1, 2) if value[1] == None else (2, 1)
                readType = ["forward", "reverse"]
                
                raise KeyError(f"An expected pair of files with prefix '{filePrefix}' is lacking the {readType[problemIndex-1]} " +
                               f"read as no file with the expected '_{problemIndex}' suffix was found to match the " +
                               f"'{value[goodIndex]}' file.")
    
    return found

def parse_metadata(metadataFile, columnNumbers=[1, 2], hasHeader=False):
    '''
    Returns:
        metaDict -- a dictionary with structure like:
                    {
                        'origFilePrefix1': 'newFilePrefix1',
                        'origFilePrefix2': 'newFilePrefix2',
                        ...
                    }
    '''
    columnNumbers = [ x-1 for x in columnNumbers ] # make numbers 0-based after user 1-based input
    foundRight = set()
    metaDict = {}
    
    # Parse metadata file into dict
    if metadataFile.endswith(".xlsx"):
        df = pd.read_excel(metadataFile, header=None if not hasHeader else 0)
        if (len(df.columns) - 1) < max(columnNumbers):
            raise ValueError(f"Excel file has fewer columns than --columns would suggest i.e., there should be at least " + 
                             f"{max(columnNumbers)+1} columns")
        
        for _, row in df.iterrows():
            left = row.iloc[columnNumbers[0]]
            right = row.iloc[columnNumbers[1]]
            if left in metaDict:
                    raise ValueError(f"Non-unique value '{left}' found in left column of metadata file")
            if right in metaDict:
                raise ValueError(f"Non-unique value '{right}' found in right column of metadata file; " +
                                 "this may cause conflicts in symlinked files")
            
            metaDict[left] = right
            foundRight.add(right)
    
    elif metadataFile.endswith(".tsv") or metadataFile.endswith(".csv"):
        delim = "\t" if metadataFile.endswith(".tsv") else ","
        with open(metadataFile, "r") as fileIn:
            skipHeader = hasHeader
            for line in fileIn:
                if skipHeader:
                    skipHeader = False
                    continue
                
                sl = line.rstrip().split(delim)
                if (len(sl) - 1) < max(columnNumbers):
                    raise ValueError(f"Line '{line.rstrip()}' when split by the delimiter '{delim}' has fewer columns " +
                                     "than --columns would suggest i.e., there should be at least " + 
                                     f"{max(columnNumbers)+1} columns")
                left = sl[columnNumbers[0]]
                right = sl[columnNumbers[1]]
                if left in metaDict:
                    raise ValueError(f"Non-unique value '{left}' found in left column of metadata file")
                if right in metaDict:
                    raise ValueError(f"Non-unique value '{right}' found in right column of metadata file; " +
                                     "this may cause conflicts in symlinked files")
                
                metaDict[left] = right
                foundRight.add(right)
    else:
        raise ValueError(f"-m file '{metadataFile}' must end in .xlsx or .tsv or .csv to facilitate parsing")
    
    # Run post-processing checks for data validity
    if metaDict == {}:
        raise ValueError(f"Parsing of -m file '{metadataFile}' resulted in no metadata associations; is this file empty?")
    
    return metaDict

def match_files_to_metadata(metaDict, foundFiles):
    '''
    Updates the foundFiles dict to match up with the metaDict keys where minor disagreements might
    exist, but where those disagreements to not prevent a 1-to-1 match between metadata and file names.
    
    Returns:
        resolvedFoundFiles -- a dictionary with structure equivalent to the original foundFiles but with
                              key values updated where relevant
    '''
    metaSet = set(metaDict.keys())
    foundSet = set(foundFiles.keys())
    
    # Identify common overlaps which do not need resolution
    foundInCommon = metaSet.intersection(foundSet)
    resolvedFoundFiles = { key : metaDict[key] for key in foundFiles } # to be further populated with keys which are matched to the metaDict keys
    
    # Identify file names which explainably differ from metadata labels [through a common base suffix]
    metaDiff = metaSet.difference(foundSet)
    foundDiff = foundSet.difference(metaSet)
    if len(foundDiff) != 0:
        # Drop our list of found differences to only those which could conceivably be rescued
        foundDiff = list(set([ x for x in foundDiff for y in metaDiff if x.startswith(y) ])) # make it a list to ensure consistent iteration order
        
        # Find a common suffix among samples which could conceivably be rescued by removal of this suffix
        commonSuffix = os.path.commonprefix([ x[::-1] for x in foundDiff ])[::-1]
        foundSansSuffix = [ x.rsplit(commonSuffix, maxsplit=1)[0] for x in foundDiff ]
        
        # Rescue samples with differences if possible
        if all([ x in metaDiff for x in foundSansSuffix ]) and ( not any([ x in foundSet or x in resolvedFoundFiles for x in foundSansSuffix ]) ):
            # i.e., if all of these are among our 'not found in metadata' samples and they are not found in the original file names
            for prefix, sansSuffix in zip(foundDiff, foundSansSuffix):
                resolvedFoundFiles[sansSuffix] = foundFiles[prefix]
                foundSet.remove(prefix)
                if sansSuffix not in foundSet:
                    foundSet.add(sansSuffix)
                else:
                    raise ValueError(f"In the process of rescuing '{prefix}' by considering it instead to be '{sansSuffix}', " +
                                     "we ended up with sample name duplication. You should probably make sure your file name prefixes " +
                                     "are maximally descriptive to prevent this rescuing process from being needed here.")
    
    # Raise error for any remaining unresolved differences
    metaDiff = metaSet.difference(foundSet) # foundSet was updated so metaDiff needs to be recalculated
    if len(metaDiff) != 0:
        formattedDiff = ", ".join(sorted(metaDiff))
        raise ValueError(f"Metadata file indicates file prefixes that were not found in any input location(s); " +
                         f"these include '{formattedDiff}'. Maybe check that the --header flag is correctly set for your metadata file, " + 
                         "or that your -s suffix value is appropriate.")
    
    # Drop any found files which are not noted in our metadata
    foundDiff = foundSet.difference(metaSet)
    if len(foundDiff) != 0:
        resolvedFoundFiles = { key:value for key, value in resolvedFoundFiles.items() if key not in foundDiff }
    
    return resolvedFoundFiles

def main():
    # User input
    usage = """%(prog)s handles the symlinking of files with generic/lab sample names into more informative
    plant/accession identifiers. Where multiple files exist with a common prefix, this script assumes they
    are paired sequencing reads and should have a '_1' or '_2' immediately before the indicated file suffix.
    The input metadata file type is inferred by file suffix of '.xlsx' or '.tsv' or '.csv'; presence of headers
    and column numbers to obtain data from have assumed defaults but can be explicitly set using the
    --header or --columns arguments.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputLocations",
                   required=True,
                   nargs="+",
                   help="Input one or more file names and/or directories containing files for symlinking")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Input metadata file with two columns indicating current file prefix and intended output prefix")
    p.add_argument("-s", dest="fileSuffix",
                   required=True,
                   help="""Suffix which uniquely identifies all relevant read files
                   e.g., '.fq.gz' for sequencing reads or '.sorted.bam' for mapped data""")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write symlinks")
    # Opts (minimap2)
    p.add_argument("--header", dest="hasHeader",
                   required=False,
                   action="store_true",
                   help="Metadata file has a header line",
                   default=False)
    p.add_argument("--columns", dest="columnNumbers",
                   required=False,
                   type=int,
                   nargs=2,
                   help="""Column numbers for left (original) and right (new) pairs;
                   use 1-based numbering (first column == 1); default=='1 2'""",
                   default=[1, 2])
    p.add_argument("--newSuffix", dest="newSuffix",
                   required=False,
                   help="""Optionally specify a new suffix to append to symlinked files; default
                   is to use whatever -s has""",
                   default=None)
    args = p.parse_args()
    validate_args(args)
    
    # Parse the metadata file
    metaDict = parse_metadata(args.metadataFile, args.columnNumbers, args.hasHeader)
    
    # Locate file(s) in the specified location(s)
    foundFiles = locate_files(args.inputLocations, args.fileSuffix)
    
    # Tolerantly match metadata and files where minor differences might occur
    "For example, a suffix of 'S1' might be given in metadata but all files might look more like 'S1.trimmed'"
    foundFiles = match_files_to_metadata(metaDict, foundFiles)
    
    # Link each file
    for origPrefix, location in foundFiles.items():
        newPrefix = metaDict[origPrefix]
        if isinstance(location, str):
            original = [ location ]
            new = [ os.path.join(args.outputDirectory, f"{newPrefix}{args.newSuffix}") ]
        else:
            original = [ location[1], location[2] ]
            new = [
                os.path.join(args.outputDirectory, f"{newPrefix}_1{args.newSuffix}"),
                os.path.join(args.outputDirectory, f"{newPrefix}_2{args.newSuffix}")
            ]
        
        for o, n in zip(original, new):
            if not os.path.exists(n):
                os.symlink(o, n)
            else:
                print(f"# Symlink for '{o}' already exists at intended output '{n}'")

    # Print completion flag if we reach this point
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
