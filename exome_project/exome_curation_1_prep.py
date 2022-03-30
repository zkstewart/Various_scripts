#! python3
# exome_curation_prep.py
# Program to enable manual curation of exome sequencing
# to occur for the Oz Mammals Genomics initiative as part
# of Matthew Phillips and Andrew Baker (et. al.'s) group.

import sys, argparse, os, math
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO
from exome_liftover import ssw_parasail

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
        metadataID = description.split(" ")[1].rsplit("_", maxsplit=1)[0] # Get just the middle part sans _S### suffix
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

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is process these files to have more informative IDs,
    to ensure every file contains all species (even if they're just blank lines), and to keep
    them in a consistent ordering.
    
    Note: This should be step 1 in the Oz Mammals project!
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-m", dest="metadataFile", required=True,
                help="Specify the metadata file location")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (default == \"1_prep\")",
                default="prep")
    # Opts
    p.add_argument("-c", dest="chunkSize", type=int, required=False,
                help="Optionally, specify how many exons should be concatenated in a file",
                default=50)
    p.add_argument("--codons", dest="showCodonPositions", required=False, action="store_true",
                help="Optionally, specify this flag if you want to produce dummy sequences containing codon numbers",
                default=False)

    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Parse metadata file
    metadataDict = get_dasyurid_metadata_dict(args.metadataFile)

    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        fastaObjs.append(f)
    
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
        file = files[i]
        outputFileName = os.path.join(args.outputDir, file)
        
        FASTA_obj = fastaObjs[i]
        FASTA_obj.write(outputFileName, withAlt=True, asAligned=True, withConsensus=False)
    
