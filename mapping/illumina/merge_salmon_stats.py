#! python3
# merge_star_stats.py
# A simple python program to combine STAR output ReadsPerGene.out.tab
# files into a table that's amenable to my way of doing DGE analysis

import os, argparse

# Define classes
class DGEQuantCollection():
    '''
    Adapted from the BINge class with some modification to ideally make it quicker
    and less memory intensive.
    '''
    def __init__(self):
        self.samples = [] # sample names
        self.quant = []
        self.length = []
        self.tpm = []
        
        self.numTranscripts = None # for checking quant file compatibility
    
    def parse_quant_file(self, quantFile, sample):
        # Validate input parameters
        assert os.path.isfile(quantFile), \
            f"Cannot parse '{quantFile}' as file does not exist!"
        assert isinstance(sample, str) and not sample in self.samples, \
            f"Sample name {sample} must be a string value that uniquely identifies this sample!"
        
        # Parse the file                
        with open(quantFile, "r") as fileIn:
            # Check header line
            sl = next(fileIn).rstrip("\r\n ").split("\t")
            assert sl == ["Name", "Length", "EffectiveLength", "TPM", "NumReads"], \
                "Quant file appears to lack the expected header? Cannot parse."
            
            # Handle content lines
            lineIndex = 0
            for line in fileIn:
                name, _, effectiveLength, tpm, numReads  = line.rstrip("\r\n ").split("\t") # might error here if file format is bad
                if self.numTranscripts == None:
                    self.quant.append([name, numReads])
                    self.length.append([name, effectiveLength])
                    self.tpm.append([name, tpm])
                else:
                    try:
                        self.quant[lineIndex].append(numReads)
                        self.length[lineIndex].append(effectiveLength)
                        self.tpm[lineIndex].append(tpm)
                    except IndexError:
                        print(f"File '{quantFile}' is longer than any we've seen previously")
                        print("This has led to an IndexError which indicates incompatibility " +
                              "in your files")
                        print("Make sure you've specified Salmon files mapped to the same " +
                              "FASTA reference; program exiting now.")
                        quit()
                
                lineIndex += 1
        
        # Check compatibility of quant file with existing ones
        if self.numTranscripts == None:
            self.numTranscripts = lineIndex
        else:
            assert self.numTranscripts == lineIndex, \
                "Quant files are incompatible since transcript numbers differ!"
        
        # Store relevant parameters now that parsing has completed successfully
        self.samples.append(sample)

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.salmonOutputsDir):
        print('I am unable to locate the parent directory where salmon subdirectories are (' + args.salmonOutputsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    args.outputFileNames = []
    for suffix in [".length", ".counts", ".abundance"]:
        outputFileName = args.outputPrefix.rstrip("._ ") + suffix
        if os.path.isfile(outputFileName):
            print(f'File already exists at output location ({outputFileName})')
            print('Make sure you specify a unique file name and try again.')
            quit()
        args.outputFileNames.append(outputFileName)

def write_dgequant_output(sampleList, dataList, outputFileName):
    with open(outputFileName, "w") as fileOut:
        fileOut.write("\t{0}\n".format("\t".join(sampleList)))
        for dataRow in dataList:
            fileOut.write("{0}\n".format("\t".join(dataRow)))

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing subdirectories that house
    salmon output files, notably the quant.sf file, and writes several files
    as output which can be loaded by DESeq2 for DGE analysis.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="salmonOutputsDir",
                   required=True,
                   help="Input directory containing subdirectories where salmon was run")
    p.add_argument("-o", dest="outputPrefix",
                   required=True,
                   help="Output file prefix for the DGE-amenable outputs")
    ## Optional
    p.add_argument("--split", dest="splitString",
                   required=False,
                   help="""Optionally, if your salmon output directories contain a string
                   you want to split at (getting anything that comes beforehand as your
                   sample ID), specify it here.""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Get quant.sf files list
    quantFiles = []
    for fileOrDir in os.listdir(args.salmonOutputsDir):
        if os.path.isdir(os.path.join(args.salmonOutputsDir, fileOrDir)):
            for subFileOrDir in os.listdir(os.path.join(args.salmonOutputsDir, fileOrDir)):
                if subFileOrDir == "quant.sf":
                    quantFiles.append(os.path.join(
                        args.salmonOutputsDir, fileOrDir, subFileOrDir
                    ))
    
    # Load in quant.sf files
    dgeQuantCollection = DGEQuantCollection()
    for quantFile in quantFiles:
        sample = os.path.basename(os.path.dirname(quantFile)) # get the dir quant.sf is contained within
        if args.splitString != None:
            sample = sample.split(args.splitString)[0]
        
        dgeQuantCollection.parse_quant_file(quantFile, sample)
    
    # Write outputs
    lengthFile, countFile, abundanceFile = args.outputFileNames
    write_dgequant_output(dgeQuantCollection.samples, dgeQuantCollection.length, lengthFile)
    write_dgequant_output(dgeQuantCollection.samples, dgeQuantCollection.quant, countFile)
    write_dgequant_output(dgeQuantCollection.samples, dgeQuantCollection.tpm, abundanceFile)
    
    # Alert user to program success with some QC statistics
    print("# merge_salmon_stats report")
    print(f"# {len(quantFiles)} salmon outputs were identified")
    print(f"# {dgeQuantCollection.numTranscripts} genes were mapped against")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
