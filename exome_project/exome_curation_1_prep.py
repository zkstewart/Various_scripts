#! python3
# exome_curation_prep.py
# Program to enable manual curation of exome sequencing
# to occur for the Oz Mammals Genomics initiative as part
# of Matthew Phillips and Andrew Baker (et. al.'s) group.

import sys, argparse, os, platform
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO, ZS_Utility

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print('I am unable to locate the metadata file (' + args.metadataFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.liftoverDirs != []:
        for loDir in args.liftoverDirs:
            if not os.path.isdir(loDir):
                print('I am unable to locate the directory where the liftover files are (' + loDir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
    
    # Validate MAFFT arguments
    if args.mafft is None:
        args.mafft = ZS_Utility.wsl_which("mafft")
        if args.mafft is None:
            print(f"ERROR: 'mafft' not discoverable in your system PATH and was not specified as an argument.")
            quit()
    else:
        if not os.path.isfile(args.mafft):
            print(f"ERROR: 'mafft' was not found at the location indicated ('{args.mafft}')")
            quit()
    
    # Handle file output
    if os.path.isdir(args.outputDir):
        if os.listdir(args.outputDir) != []:
            print(args.outputDir + ' already contains files. Specify a new location or move any existing files elsewhere.')
            quit()
    else:
        try:
            os.mkdir(args.outputDir)
            print("Created '{0}' directory as part of argument validation".format(args.outputDir))
        except:
            print("Wasn't able to create '{0}' directory; does '{1}' actually exist?".format(args.outputDir, os.path.dirname(args.outputDir)))

def get_dasyurid_metadata_dict(metadataFile):
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ").split("\t")
            metadataDict[l[0]] = l[1]
    return metadataDict

def set_alts(FASTA_obj, metadataDict):
    '''
    This function will receive a single FASTA object and perform
    Oz Mammals-specific processing to get an alt ID appropriately
    set for each FastASeq object.
    '''
    for FastASeq_obj in FASTA_obj:
        description = FastASeq_obj.description
        if " " in description:
            metadataID = description.split(" ")[1].rsplit("_", maxsplit=1)[0] # Get just the middle part sans _S### suffix
        else:
            metadataID = description
        altID = metadataDict[metadataID]
        FastASeq_obj.alt = altID

def add_missing_seqs(FASTA_obj, sequenceIDs):
    '''
    This function will receive a single FASTA object and perform
    Oz Mammals-specific processing to add missing sequences into
    the FASTA. The sequenceIDs object should be a list derived
    from the metadataDict's keys.
    '''
    # Figure out which IDs are missing
    altIDs = [FastASeq_obj.alt for FastASeq_obj in FASTA_obj]
    missingIDs = [id for id in sequenceIDs if id not in altIDs]
    
    # Add dummy sequences to FASTA object
    mockSequence = "-" * len(FASTA_obj[0].gap_seq) # Mock up a fully-gapped sequence
    for mID in missingIDs:
        dummyFastASeq_obj = ZS_SeqIO.FastASeq(mID, alt=mID, gapSeq = mockSequence)
        FASTA_obj.insert(0, dummyFastASeq_obj) # Just insert at index = 0, we'll sort things later

def main():
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is process these files to have more informative IDs,
    to ensure every file contains all species (even if they're just blank lines), and to keep
    them in a consistent ordering.
    
    Note: This should be step 1 in the Oz Mammals project!
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir",
                   required=True,
                   help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Specify the metadata file location")
    p.add_argument("--mafft", dest="mafft",
                   required=False,
                   help="""Optionally, specify the mafft executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("-o", dest="outputDir",
                   required=False,
                   help="Output directory location (default == \"1_prep\")",
                   default="1_prep")
    # Opts
    p.add_argument("-l", dest="liftoverDirs",
                   required=False,
                   nargs="+",
                   help="Optionally, specify one or more directories where exome liftover FASTAs can be found",
                   default=[])
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Parse metadata file
    metadataDict = get_dasyurid_metadata_dict(args.metadataFile)
    
    # Load FASTA files
    fastaObjs = []
    mergeStats = [0 for _ in range(len(args.liftoverDirs))]
    for file in files:
        # Locate the base alignment file
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        
        # Opt: Add in liftover FASTA file contents
        if args.liftoverDirs != []:
            
            # Locate relevant liftover file(s)
            '''
            This is ugly, but we have files with .fa and .fasta suffix that we need
            to know are matched. So we need to handle just the prefixes, then reattach
            the suffix later
            '''
            suffixlessBaseFile = os.path.basename(file).rsplit(".", maxsplit=1)[0] 
            loFiles = []
            for i in range(len(args.liftoverDirs)):
                loDir = args.liftoverDirs[i]
                suffixlessLoFiles = [x.rsplit(".", maxsplit=1)[0] for x in os.listdir(loDir)] # agnostic to file suffix
                loFileSuffixes = [x.rsplit(".", maxsplit=1)[1] for x in os.listdir(loDir)] # paired suffixes list to the above
                
                if suffixlessBaseFile in suffixlessLoFiles:
                    _relevantSuffix = loFileSuffixes[suffixlessLoFiles.index(suffixlessBaseFile)]
                    loFiles.append(os.path.join(loDir, suffixlessBaseFile + "." + _relevantSuffix))
                    mergeStats[i] += 1 # Keep track of how many sequences we've merged from this liftover file
            
            # Align and merge relevant file(s) into FASTA
            for _loFile in loFiles:
                add_FASTA_obj = ZS_SeqIO.FASTA(_loFile)
                
                # Perform MAFFT --add alignment
                m = ZS_AlignIO.MAFFT(args.mafft)
                f = m.add(f, add_FASTA_obj)
        
        # Store in fasta objects list
        fastaObjs.append(f)
    
    # Opt: Print merge stats details if applicable
    if args.liftoverDirs != []:
        print("Merging details:")
        for i in range(len(mergeStats)):
            loDir = args.liftoverDirs[i]
            stat = mergeStats[i]
            print("\t{0} sequences merged from {1}".format(stat, loDir))
    
    # Set FASTA alt IDs
    for FASTA_obj in fastaObjs:
        set_alts(FASTA_obj, metadataDict)
    
    # Add dummy missing sequences
    sequenceIDs = list(metadataDict.values())
    for FASTA_obj in fastaObjs:
        add_missing_seqs(FASTA_obj, sequenceIDs)
    
    # Sort FASTA objects to have consistent internal ordering
    for FASTA_obj in fastaObjs:
        FASTA_obj.seqs.sort(key = lambda x: sequenceIDs.index(x.alt))
    
    # Write output files
    for i in range(len(files)):
        file = os.path.basename(files[i])
        outputFileName = os.path.join(args.outputDir, file)
        
        FASTA_obj = fastaObjs[i]
        FASTA_obj.write(outputFileName, withAlt=True, asAligned=True, withConsensus=False)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
