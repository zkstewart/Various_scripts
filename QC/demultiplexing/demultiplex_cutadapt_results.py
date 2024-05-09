#! python 3
# demultiplex_cutadapt_results.py
# A script to receive a directory resulting from cutadapt's demultiplexing
# and _properly_ demultiplex the reads into separate files based on the
# sample name.

import os, argparse, shutil

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
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
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
    return metaDict

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
    
    # Loop through cutadapt files to locate and copy reads files to the output dir
    foundSamples = set()
    numFiles = 0
    unknownSamples = 0
    for file in os.listdir(args.fastqDirectory):
        if file.endswith(".fastq"):
            numFiles += 1
            
            # Get the forward and reverse adapters
            adapters = file.split(".")[0].split("-")
            if len(adapters) != 2:
                print(f"Error: The filename '{file}' does not have the expected forward-reverse adapter format; " +
                      "please check the file and try again.""")
                quit()
            forward, reverse = adapters
            
            # Get the sample name
            sample = metaDict.get(forward, {}).get(reverse, None)
            if sample is None:
                unknownSamples += 1
                continue
            foundSamples.add(sample)
            
            # Identify whether it's the forward (.1) or reverse (.2) file
            if file.endswith(".1.fastq"):
                directionNum = 1
            elif file.endswith(".2.fastq"):
                directionNum = 2
            else:
                print(f"Error: The filename '{file}' does not have the expected direction indicator (.1.fastq or .2.fastq); " +
                      "please check the file and try again.")
                quit()
            
            # Copy the file to the output directory
            newFileName = os.path.join(args.outputDirectory, f"{sample}.{directionNum}.fastq")
            if os.path.isfile(newFileName):
                print(f"'{newFileName}' already exists; skipping...")
            else:
                shutil.copyfile(os.path.join(args.fastqDirectory, file), newFileName)
    
    # Print statistics then end program
    notFoundSamples = sampleIDs.difference(foundSamples)
    print("# demultiplex_cutadapt_results output statistics")
    print(f"# Processing of files located at: '{args.fastqDirectory}'")
    print(f"# > Of {len(sampleIDs)} samples noted in your metadata file, I identified {len(foundSamples)} and copied them to your output directory")
    if len(notFoundSamples) != 0:
        print(f"# > The following samples were not found in the cutadapt output: {', '.join(notFoundSamples)}")
    print(f"# > Additionally, of {numFiles} cutadapt files, {len(foundSamples)} had recognised forward-reverse adapter pairings")
    print(f"# > This leaves {unknownSamples} or {(unknownSamples / numFiles) * 100}% of demultiplexed files as being unknown")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
