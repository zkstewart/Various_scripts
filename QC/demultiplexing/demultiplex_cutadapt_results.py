#! python 3
# demultiplex_cutadapt_results.py
# A script to receive a directory resulting from cutadapt's demultiplexing
# and _properly_ demultiplex the reads into separate files based on the
# sample name.

import os, argparse

def validate_args(args):
    # Validate input data
    if not os.path.isdir(args.fastqDirectory):
        print(f'I am unable to locate the cutadapt FASTQ results directory ({args.fastqDirectory})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.demultiplexingMetadata):
        print(f'I am unable to locate the TSV metadata file ({args.demultiplexingMetadata})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        if os.listdir(args.outputDirectory) != []:
            print("WARNING: I won't overwrite any existing files, so beware that if a previous run had issues, " +
                "you should have deleted or moved those files first.")
    else:
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def parse_cutadapt_metadata(metadataFile):
    '''
    Parse the metadata file and return a dictionary for use with associating cutadapt results
    to their intended sample identifier.
    
    Parameters:
        metadataFile -- a string indicating the location of a metadata TSV to load in; it is expected
                        to contain three columns and a header which WILL BE IGNORED; the columns
                        should be 'Sample', 'Forward adapter', and 'Reverse adapter'.
    Returns:
        metaDict -- a dictionary with structure like:
                    {
                        "forward_adapter1": {
                            "reverse_adapter1": "Sample1",
                            "reverse_adapter2": "Sample2",
                            ...
                        },
                        "forward_adapter2": {
                            "reverse_adapter1": "Sample3"
                        },
                        ...
                    }
    '''
    metaDict = {}
    with open(metadataFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            if firstLine:
                firstLine = False
            else:
                sl = line.strip("\r\n ").split("\t")
                if len(sl) != 3:
                    print(f"Error: Metadata file '{metadataFile}' does not have three columns; please check the file and try again.")
                    quit()
                
                sample, forward, reverse = sl
                metaDict.setdefault(forward, {})
                metaDict[forward][reverse] = sample
    
    # Index adapters in reverse
    reverseDict = {}
    for forward, _revDict in metaDict.items():
        for reverse, sample in _revDict.items():
            if reverse in metaDict:
                assert forward not in metaDict[reverse]
            reverseDict.setdefault(reverse, {})
            reverseDict[reverse][forward] = sample
    
    # Merge dicts together
    metaDict.update(reverseDict)
    
    return metaDict

def write_merged_files(componentFileNames, outputFileName):
    '''
    Helper function to write a series of files to a single output file.
    
    Parameters:
        componentFileNames -- a list of strings indicating the file names to merge together
        outputFileName -- a string indicating the name of the output file to write to
    '''
    with open(outputFileName, "w") as fileOut:
        for f in componentFileNames:
            if os.path.isfile(f):
                with open(f, "r") as fileIn:
                    for line in fileIn:
                        fileOut.write(line)

def main():
    usage = """%(prog)s will receive a directory containing cutadapt's demultiplexing output.
    It is expected that the directory contains a series of FASTQ files, each of which
    has a name with the format {{name1}}-{{name2}}.1.fastq and {{name1}}-{{name2}}.2.fastq for forward
    and reverse reads; this format should have been specified to cutadapt, and it results in the file name
    having the two adapters separated by a hyphen. Alongside a metadata TSV file indicating the forward
    and reverse adapters associated with each sample, this script will identify the reads associated with
    each sample and write them with the appropriate sample name to the output directory.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="fastqDirectory",
                   required=True,
                   help="Specify the location of the cutadapt results directory")
    p.add_argument("-m", dest="demultiplexingMetadata",
                   required=True,
                   help="""Specify a TSV file containing the metadata for demultiplexing; it is
                   expected to contain three columns: 'Sample', 'Forward adapter', and 'Reverse adapter';
                   headers are expected and the first line will be ignored.""")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory for demultiplexed FASTQ files")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the metadata file
    metaDict = parse_cutadapt_metadata(args.demultiplexingMetadata)
    sampleIDs = set([ sampleID for subdict in metaDict.values() for sampleID in subdict.values() ])
    
    # Loop through samples and write their demultiplexed reads to file
    foundSamples = set()
    knownFileSizes = 0
    for sampleID in sampleIDs:
        # Locate the forward and reverse adapters for the sample
        adapters = [
            [ forward, reverse ]
            for forward, _revDict in metaDict.items()
            for reverse, sample in _revDict.items()
            if sample == sampleID
        ]
        
        # Format the forward and reverse paired file names expected
        fwdFiles = [ os.path.join(args.fastqDirectory, f"{forward}-{reverse}.1.fastq") for forward, reverse in adapters ]
        revFiles = [ os.path.join(args.fastqDirectory, f"{forward}-{reverse}.2.fastq") for forward, reverse in adapters ]
        
        # Check that any pairs occur together
        foundSomething = False
        for fwdFile, revFile in zip(fwdFiles, revFiles):
            if os.path.isfile(fwdFile):
                foundSomething = True
                knownFileSizes += os.stat(fwdFile).st_size / (1024 * 1024)
                
                # Make sure the reverse pair also occurs
                assert os.path.isfile(revFile), \
                    (f"Error: The forward file '{fwdFile}' exists, but the reverse file '{revFile}' does not; " + 
                     "please check that your input files are correct and try again.")
            else:
                # Make sure the reverse pair is also absent
                assert not os.path.isfile(revFile), \
                    (f"Error: The reverse file '{revFile}' exists, but the forward file '{fwdFile}' does not; " +
                     "please check that your input files are correct and try again.")
        
        # If no pairs occurred, continue now
        if foundSomething == False:
            continue
        foundSamples.add(sampleID) # track that we found this sample
        
        # Otherwise, write the file(s) to the output directory
        newFwdFile = os.path.join(args.outputDirectory, f"{sampleID}.1.fastq")
        if os.path.isfile(newFwdFile):
            print(f"'{newFwdFile}' already exists; skipping...")
        else:
            write_merged_files(fwdFiles, newFwdFile)
        
        newRevFile = os.path.join(args.outputDirectory, f"{sampleID}.2.fastq")
        if os.path.isfile(newRevFile):
            print(f"'{newRevFile}' already exists; skipping...")
        else:
            write_merged_files(revFiles, newRevFile)
    
    # Calculate the file size for all demultiplexed outputs
    totalFileSize = 0
    for file in os.listdir(args.fastqDirectory):
        if file.endswith(".1.fastq"):
            totalFileSize += os.stat(os.path.join(args.fastqDirectory, file)).st_size / (1024 * 1024)
    
    # Print statistics then end program
    notFoundSamples = sampleIDs.difference(foundSamples)
    
    print("# demultiplex_cutadapt_results output statistics")
    print(f"# Processing of files located at: '{args.fastqDirectory}'")
    print(f"# > Of {len(sampleIDs)} samples noted in your metadata file, I identified {len(foundSamples)} and wrote them to your output directory")
    if len(notFoundSamples) != 0:
        print(f"# > The following samples were not found in the cutadapt output: {', '.join(notFoundSamples)}")
    print(f"# > Approximately {(knownFileSizes / totalFileSize) * 100}% of your reads were associated with known samples")
    print(f"# > The remaining reads are considered 'unknown' by cutadapt")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
