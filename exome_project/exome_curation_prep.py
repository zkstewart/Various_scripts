#! python3
# exome_curation_prep.py
# Program to enable manual curation of exome sequencing
# to occur for the Oz Mammals Genomics initiative as part
# of Matthew Phillips and Andrew Baker (et. al.'s) group.

import sys, argparse, os, math
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO

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
    # Validate numeric inputs
    if args.chunkSize < 1:
        print("chunkSize should be at least 1")
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

def get_missing_species(fastaFiles, metadataDict):
    '''
    Function is just for debugging and getting an updated metadata file
    '''
    missingIDs = set()
    for fastaFile in fastaFiles:
        FASTA_obj = ZS_SeqIO.FASTA(fastaFile)
        for FastASeq_obj in FASTA_obj:
            description = FastASeq_obj.description
            metadataID = description.split(" ")[1].rsplit("_", maxsplit=1)[0] # Get just the middle part sans _S### suffix
            if metadataID in metadataDict:
                continue
            else:
                missingIDs.add(metadataID)
    return list(missingIDs)

def set_alts(FASTA_obj):
    '''
    This function will receive a single FASTA object and perform
    Oz Mammals-specific processing to get an alt ID appropriately
    set for each FastASeq object.
    '''
    for FastASeq_obj in FASTA_obj:
        description = FastASeq_obj.description
        metadataID = description.split(" ")[1].rsplit("_", maxsplit=1)[0] # Get just the middle part sans _S### suffix
        FastASeq_obj.alt = metadataID

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

def get_chunking_points(numberToChunk, chunkSize):
    '''
    This is a general purpose function to take in a number of "things"
    that you want to chunk, and find out how to chunk them evenly.
    
    Params:
        numberToChunk -- an integer value, possibly derived from a list length as example.
        chunkSize -- an integer value for the desired number of things per chunk.
    '''
    assert isinstance(numberToChunk, int)
    assert isinstance(chunkSize, int)
    if numberToChunk <= chunkSize:
        raise Exception("Chunking only valid if chunkSize is smaller than numberToChunk")
    
    numChunks = int(numberToChunk / chunkSize)
    rawNum = numberToChunk / numChunks # This line is more relevant in the multithreading code I took this from, but it's okay to just leave it.
    numRoundedUp = round((rawNum % 1) * numChunks, 0) # By taking the decimal place and multiplying it by the num of chunks, we can figure out how many chunks need to be rounded up
    
    chunkPoints = []
    ongoingCount = 0
    for i in range(numChunks):
        if i+1 <= numRoundedUp: # ngl I don't remember why this is needed; I'm borrowing this code from something I wrote a while back
            chunkPoints.append(math.ceil(rawNum) + ongoingCount) # Round up the rawNum, and also add our ongoingCount which corresponds to the number of things already put into a chunk
            ongoingCount += math.ceil(rawNum)
        else:
            chunkPoints.append(math.floor(rawNum) + ongoingCount)
            ongoingCount += math.floor(rawNum)
        if ongoingCount >= numberToChunk: # Without this check, if we have more chunks than things to chunk, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
            break  # This doesn't actually affect program function, but for aesthetic reasons and for clarity of how this function works, I prevent this from occurring.
    
    return chunkPoints

def get_mock_genename_sequence(FASTA_obj):
    '''
    Very specific to this Oz Mammals project, this function helps
    to format the gene name dummy sequence that Matt would like
    in the FASTA file
    '''
    geneName = os.path.basename(FASTA_obj.fileOrder[0][0]).split("-mx.fa")[0]
    return geneName.ljust(len(FASTA_obj[0].gap_seq), '-')

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to transform these alignments into a format that
    is friendly for manual inspection and modification.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-m", dest="metadataFile", required=True,
                help="Specify the metadata file location")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location")
    # Opts
    p.add_argument("-c", dest="chunkSize", type=int, required=False,
                help="Optionally, specify how many exons should be concatenated in a file",
                default=50)

    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Parse metadata file
    metadataDict = get_dasyurid_metadata_dict(args.metadataFile)
    
    ## Debugging: find missing species IDs
    #missingIDs = get_missing_species(files, metadataDict)
    #pyperclip.copy("\n".join(missingIDs))
    
    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        fastaObjs.append(f)
    
    # Set FASTA alt IDs
    for FASTA_obj in fastaObjs:
        set_alts(FASTA_obj)
    
    # Add dummy missing sequences
    sequenceIDs = list(metadataDict.keys())
    for FASTA_obj in fastaObjs:
        add_missing_seqs(FASTA_obj, sequenceIDs)
    
    # Sort FASTA objects to have consistent internal ordering
    for FASTA_obj in fastaObjs:
        FASTA_obj.seqs.sort(key = lambda x: sequenceIDs.index(x.alt))
    
    # Sort FASTA object list by conservation proportion & GC content
    ## TBD
    
    # Figure out how to concatenate FASTA files in chunks of ~50 (or whatever args.chunkSize is)
    chunkPoints = get_chunking_points(len(fastaObjs), args.chunkSize)
    
    # Perform the chunking
    concatFastaObjs = []
    prevChunkStart = 0 # Init as 0 for first loop iteration
    for i in range(0, len(chunkPoints)):
        baseFASTA = fastaObjs[prevChunkStart] # We use the first FASTA in this chunk as our base and concat to it
        geneNameSeq = get_mock_genename_sequence(baseFASTA)
        
        for x in range(prevChunkStart+1, chunkPoints[i]): # +1 since we're already using the first file as our base
            concatFASTA = fastaObjs[x]
            geneNameSeq += get_mock_genename_sequence(concatFASTA)
            for seqIndex in range(len(concatFASTA)):
                baseFASTA.seqs[seqIndex].extend(concatFASTA[seqIndex].gap_seq)
        
        dummyFastASeq_obj = ZS_SeqIO.FastASeq("GeneName", alt="GeneName", gapSeq = geneNameSeq)
        baseFASTA.insert(0, dummyFastASeq_obj)
        
        concatFastaObjs.append(baseFASTA) # Store our modified FASTA with all the concatenation performed
        prevChunkStart = chunkPoints[i]

    # Write output files
    for i in range(0, len(chunkPoints)):
        outputFileName = os.path.join(args.outputDir, "exons_chunk_{0}.fa".format(i+1))
        FASTA_obj = concatFastaObjs[i]
        FASTA_obj.write(outputFileName, withAlt=True, asAligned=True, withConsensus=False)
    
