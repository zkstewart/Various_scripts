#! python3
# gff3_coordinate_fixer.py
# Script to receive a GFF3 and a released nucleotide sequence
# file and will attempt to correct any GFF3 features where
# the coordinates don't match the real feature.

import os, argparse, sys, hashlib, pickle

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from Function_packages import ZS_GFF3IO, ZS_SeqIO, ZS_AlignIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFastaFile):
        print('I am unable to locate the genome FASTA file (' + args.genomeFastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.cdsFastaFile):
        print('I am unable to locate the CDS FASTA file (' + args.cdsFastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.exonerateExe):
        print('I am unable to locate the exonerate exe file (' + args.exonerateExe + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def find_nonequivalent_features(GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, isProtein=False):
    '''
    This function receives a GFF3 and a FASTA that is meant to
    correspond to the sequences extracted from that GFF3. The goal
    is to find cases where the sequences differ, and to return
    these feature/sequence IDs as a list for further work to be done on them.
    
    It will also resolve easy to solve issues pertaining to strand and frame
    annotations. If either appear to be wrong, this function will identify it
    easily! It's up to you to change that in the GFF3_obj though.
    
    Parameters:
        GFF3_obj -- a ZS_GFF3IO.GFF3 object containing features
                    with .ID values found within the cdsFASTA_obj
        cdsFASTA_obj -- a ZS_SeqIO.FASTA object containing sequences
                        with  .id values found within the GFF3_obj
        genomeFASTA_obj -- a ZS_SeqIO.FASTA object containing contig sequences
                           corresponding to .contig values found within the GFF3_obj
        isProtein -- a boolean indicating whether the sequences in the
                     FASTA object are proteins or not (i.e., are nucleotides)
    Returns:
        nonequalIDs -- a list containing ID values corresponding to
                       sequences that don't match between the GFF3 and
                       FASTA objects at all.
        fixesDict -- a dictionary with structure like:
                     {
                         'seqID1': [strand, frame],
                         'seqID2': [strand, frame],
                         ...
                     }
    '''
    nonequalIDs = []
    fixesDict = {}
    for FastASeq_obj in cdsFASTA_obj:
        assert FastASeq_obj.id in GFF3_obj, \
            f"'{FastASeq_obj.id}' not found in GFF3; can't proceed"
        
        # Get translation of FASTA sequence if relevant
        if not isProtein:
            fastaProtSequence, _, _ = FastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
        else:
            fastaProtSequence = FastASeq_obj.seq
        
        # Get translation of GFF3 feature
        feature = GFF3_obj[FastASeq_obj.id]
        gff3_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
            genomeFASTA_obj, feature.ID, "CDS"
        )
        gff3ProtSequence, _, _ = gff3_FastASeq_obj.get_translation(
            findBestFrame=False,
            strand=1 if feature.strand is "+" else -1,
            frame=int(feature.frame) if feature.frame.isdigit() else 0
        )
        
        # Find equivalency
        isEqual = True if gff3ProtSequence.rstrip("*") == fastaProtSequence.rstrip("*") else False
        if not isEqual:
            # Check if our findBestFrame approach works
            alt_gff3ProtSequence, alt_strand, alt_frame = gff3_FastASeq_obj.get_translation(
                findBestFrame=True
            )
            alt_isEqual = True if alt_gff3ProtSequence.rstrip("*") == fastaProtSequence.rstrip("*") else False
            
            if not alt_isEqual: # just put it in the "this is fucked" basket
                nonequalIDs.append(FastASeq_obj.id)
            else:
                fixesDict[FastASeq_obj.id] = [alt_strand, alt_frame]
    return nonequalIDs, fixesDict

def fix_nonequivalent_features(nonequalIDs, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, isProtein=False):
    '''
    Intended to follow on from find_nonequivalent_features(), this function will
    instead aim to correct the coordinates of any GFF3 features that don't match
    their respective CDS gene model. Changes will be directly made to the input
    GFF3_obj.
    
    Parameters:
        nonequalIDs -- a list containing ID values corresponding to
                       sequences that don't match between the GFF3 and
                       FASTA objects at all.
        GFF3_obj -- a ZS_GFF3IO.GFF3 object containing features
                    with .ID values found within the cdsFASTA_obj
        cdsFASTA_obj -- a ZS_SeqIO.FASTA object containing sequences
                        with  .id values found within the GFF3_obj
        genomeFASTA_obj -- a ZS_SeqIO.FASTA object containing contig sequences
                           corresponding to .contig values found within the GFF3_obj
        isProtein -- a boolean indicating whether the sequences in the
                     FASTA object are proteins or not (i.e., are nucleotides)
    '''
    for problemSeqID in nonequalIDs:
        FastASeq_obj = cdsFASTA_obj[problemSeqID]
        feature = GFF3_obj[problemSeqID]
        
        # Get translation of FASTA sequence if relevant
        if not isProtein:
            fastaProtSequence, _, _ = FastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
        else:
            fastaProtSequence = FastASeq_obj.seq
        
        # 
        stop


def get_args_hash(args):
    '''
    Simple function to receive the args object (but theoretically any will do)
    and return a consistent hash (length == 20 characters) based on its contents.
    
    Parameters:
        args -- an ArgumentParser.parse_args() response object, but also just
                any object should work, too.
    Returns:
        argsHash -- a sha256 hash of the args which will always be identical when
                    the args parameters are the same.
    '''
    hashes = []
    for key, value in args.__dict__.items():
        hashes.append(hashlib.sha256(bytes(str(key) + str(value), 'utf-8')).hexdigest())
    overallHash = hashlib.sha256(bytes("".join(hashes), 'utf-8')).hexdigest()
    return overallHash[0:20]

def main():
    # User input
    usage = """%(prog)s accepts a GFF3 file and ...
    
    Note that the IDs of everything are expected to match. That means you might need to do some
    prep work on some of the FASTA files to make their IDs match the GFF3.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-g", dest="gff3File", required=True,
                help="Specify the location of the GFF3 file")
    p.add_argument("-f", dest="genomeFastaFile", required=True,
                help="Specify the location of the genome FASTA file")
    p.add_argument("-c", dest="cdsFastaFile", required=True,
                help="Specify the location of the CDS nucleotide/protein FASTA file")
    p.add_argument("-e", dest="exonerateExe", required=True,
                help="Specify the location of the exonerate executable file")
    p.add_argument("-o", dest="outputFileName", required=True,
                help="Specify the location to write the modified GFF3 file to")
    # Opts
    p.add_argument("--isProtein", dest="isProtein", required=False, action="store_true",
                help="Optionally specify that the cdsFastaFile contains protein sequences",
                default=False)
    p.add_argument("--resume", dest="resume", required=False, action="store_true",
                help="""Optionally provide this tag to save a pickle of the exonerate
                results which can be resumed from""",
                default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Figure out the hash and exonerate pickle file for this run
    argsHash = get_args_hash(args)
    exoneratePickleFile = os.path.join(
        os.path.dirname(args.outputFileName),
        f"exonerate_results_{argsHash}.pkl"
    )

    # Parse GFF3
    GFF3_obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False) # non-strict parsing
    
    # Parse FASTA files
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.genomeFastaFile)
    cdsFASTA_obj = ZS_SeqIO.FASTA(args.cdsFastaFile)
    
    # Get exonerate results depending on whether we can resume or not
    if args.resume is False or not os.path.isfile(exoneratePickleFile):
        # Get a FASTA file for the sequences we need to search against the genome
        nonequalIDs, fixesDict = find_nonequivalent_features(GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, args.isProtein)
        nonequalFASTA_obj = ZS_SeqIO.FASTA(None)
        nonequalIDs = nonequalIDs[0:5]
        for seqID in nonequalIDs:
            nonequalFASTA_obj.add(cdsFASTA_obj[seqID])
        
        tmpFileName = ZS_AlignIO._tmp_file_name_gen("tmpExonerateQuery", "fasta")
        nonequalFASTA_obj.write(tmpFileName)
        
        # Perform exonerate search
        exonerateSearcher = ZS_AlignIO.Exonerate(args.exonerateExe, tmpFileName, args.genomeFastaFile)
        exonerateSearcher.model = "protein2genome"
        features = exonerateSearcher.run_exonerate() # returns a list of GFF3 features
        
        # Save result as pickle
        with open(exoneratePickleFile, "wb") as fileOut:
            pickle.dump([features, fixesDict], fileOut)
    else:
        # Load in the exonerate pickle
        with open(exoneratePickleFile, "rb") as fileIn:
            features, fixesDict = pickle.load(fileIn)
    
    # Attempt to fix difficult scenarios
    
    
    # Handle simple fixes
    # GFF3_obj, FASTA_obj = GFF3_obj, cdsFASTA_obj # just for testing
    
    # mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(GFF3_obj[geneID])
    # cds_FastASeq_obj, cds_featureType, cds_startingFrame = GFF3_obj.retrieve_sequence_from_FASTA(genomeFASTA_obj, mrnaFeature.ID, "CDS")
    
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
