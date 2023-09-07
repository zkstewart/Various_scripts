#! python3
# basespace_duplicate_symlinker.py

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.downloadDirectory):
        print('I am unable to locate the base directory (' + args.downloadDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output location
    if os.path.isdir(args.outputDirectory):
        if len(os.listdir(args.outputDirectory)) != 0:
            print(f"Output location '{args.outputDirectory}' already exists and contains files")
            print("I don't want to risk making a mess of this. Choose a new directory or empty this one, and try again.")
            quit()
    elif os.path.exists(args.outputDirectory):
        print(f"Output location '{args.outputDirectory}' already exists and is not a directory")
        print("Make sure to specify an empty directory, or one that can be created.")
        quit()
    else:
        try:
            os.mkdir(args.outputDirectory)
            print(f"Created '{args.outputDirectory}' as part of argument validation")
        except:
            print(f"Unable to create '{args.outputDirectory}'; do intermediate directories exist?")
            print("Please check")

## Main
def main():
    # User input
    usage = """%(prog)s receives the directory where a basespace download has occurred
    which contains multiple failed FASTQ generation runs and at least one successful one.
    It will attempt to locate the successful run based on file size, and create symlinks
    to these files in the current directory.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="downloadDirectory",
                   required=True,
                   help="Input base directory where basespace download occurred")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory to write symlinks to")
    args = p.parse_args()
    validate_args(args)
    
    # Parse all subdirs looking for duplicates, and get the size of any containing files
    sampleDirs = {}
    for subLocation in os.listdir(args.downloadDirectory):
        if "_ds." in subLocation: # this string should mark a FASTQ generation output folder
            # Set up the dict
            sampleID = subLocation.split("_ds.")[0]
            fullSubLocation = os.path.join(args.downloadDirectory, subLocation)
            
            # Skip over any directories which don't fit expectations
            if not all([f.endswith(".fastq.gz") for f in os.listdir(fullSubLocation)]):
                continue
            
            sampleDirs.setdefault(sampleID, {})
            sampleDirs[sampleID][fullSubLocation] = []
            
            # Store file sizes
            for fastqFile in os.listdir(fullSubLocation):
                fullFastqFile = os.path.join(fullSubLocation, fastqFile)
                fileSize = os.path.getsize(fullFastqFile)
                
                sampleDirs[sampleID][fullSubLocation].append([fastqFile, fileSize])
    
    # Locate the best subdir when duplicates exist
    chosenSampleDirs = []
    bestSizes = []
    for sampleID, locationDict in sampleDirs.items():
        # Validate that all FASTQ files are named the same
        fileNames = None
        for fullSubLocation, fastqAndSize in locationDict.items():
            fastqFiles = [ fastq for fastq, size in fastqAndSize ]
            if fileNames == None:
                fileNames = set(fastqFiles)
            else:
                assert fileNames == set(fastqFiles), f"{fileNames} ... {set(fastqFiles)}" #"FASTQ files are not identical in same named subdir"
        
        # Pick out the dataset with the largest file size
        largest = [None, 0]
        for fullSubLocation, fastqFiles in locationDict.items():
            fastqSize = sum([ size for fastq, size in fastqAndSize ])
            largest = largest if largest[1] > fastqSize else [fullSubLocation, fastqSize]
        assert largest[0] != None, "Found no files in FASTQ fullSubLocation?"
        
        # Add this to our list
        chosenSampleDirs.append(largest[0])
        bestSizes.append(largest[1])
    
    # Symlink the file contents of our chosen locations
    for fullSubLocation in chosenSampleDirs:
        for fastqFile in os.listdir(fullSubLocation):
            fullFastqFile = os.path.join(fullSubLocation, fastqFile)
            os.symlink(fullFastqFile, os.path.join(args.outputDirectory, fastqFile))
    
    # Print some QC information for human to make sure this worked
    numWithDuplicates = sum([ 1 if len(locationDict) > 1 else 0 for _, locationDict in sampleDirs.items() ])
    
    print("# basespace_duplicate_symlink report")
    print(f"# {len(chosenSampleDirs)} unique samples were identified")
    print(f"# {numWithDuplicates} samples were found as duplicated")
    print(f"# Smallest combined size for a sample = {min(bestSizes)}")
    print(f"# Largest combined size for a sample = {max(bestSizes)}")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
