#! python3
# gff3_coordinate_fixer.py
# Script to receive a GFF3 and a released nucleotide sequence
# file and will attempt to correct any GFF3 features where
# the coordinates don't match the real feature.

import os, argparse, sys, hashlib, pickle, math

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
            strand=1 if feature.strand == "+" else -1,
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

def fix_nonequivalent_features(nonequalIDs, exonerateResultsDict, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, isProtein=False):
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
        
        # Get exonerate results if applicable
        if problemSeqID in exonerateResultsDict:
            exonerateResults = exonerateResultsDict[problemSeqID]
        else:
            continue # nothing we can do in this scenario
        
        # Find the best exonerate result in the closest proximity to the original gene annotation
        bestResults = [] # [feature, distance, lengthDifference, similarity]
        for geneFeature in exonerateResults:
            if geneFeature.contig == feature.contig:
                mrnaFeature = geneFeature.mRNA[0]
                
                # Get the distance between the genes (if any)
                if mrnaFeature.end >= feature.start and feature.end >= mrnaFeature.start:
                    distance = 0 # if they overlap, there's essentially no distance between them
                else:
                    distance = min(abs(mrnaFeature.end - feature.start), abs(feature.end - mrnaFeature.start))
                
                # Get the difference in CDS length between the genes
                mrnaFeature_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
                    genomeFASTA_obj, mrnaFeature, "CDS"
                )
                lengthDifference = abs(len(fastaProtSequence) - (len(mrnaFeature_FastASeq_obj.seq) / 3))
                
                # Finally, store the above and the similarity score
                bestResults.append([mrnaFeature, distance, lengthDifference, float(mrnaFeature.similarity)])
        bestResults.sort(key = lambda x: (-x[3], x[1], x[2]))
        bestResult = bestResults[0][0]
        stop
        # Extend the gene to its start codon, stop codon boundaries
        
        # Validate that this result has a comparable protein
        
        ## TBD...
        stop

def feature_cds_extension_maximal(mrnaFeature, mrnaFeature_FastASeq_obj, genomeFASTA_obj, MAX_EXTENSION=30):
    '''
    This function has been heavily modified from the original logic developed as part
    of the gmap_gene_find.py / exonerate_gene_find.py programs. Its goal is to take a
    GFF3.Feature object that may not have its CDS boundaries properly predicted if e.g.,
    the stop codon has not been predicted by exonerate. It will extend the model
    a length of $MAX_EXTENSION base pairs looking for an appropriate start and/or
    stop codon. The feature object will be modified if such an instance is found.
    
    Parameters:
        mrnaFeature -- a ZS_GFF3IO.Feature object corresponding to an mRNA feature
                       or something akin to that whereby CDS subfeatures are indexed.
        mrnaFeature_FastASeq_obj -- A ZS_SeqIO.FastASeq object corresponding to the
                                    mRNA feature provided.
        genomeFASTA_obj -- a ZS_SeqIO.FASTA object containing contig sequences
                           corresponding to .contig values found within the GFF3_obj
        MAX_EXTENSION -- an integer setting the maximum amount of base pairs a CDS
                         will be extended in search of a good start or stop codon.
    '''
    START_CODONS_POS = ["ATG"] #, "CTG", "TTG"]
    START_CODONS_NEG = ["CAT"] #, "CAG", "CAA"]
    STOP_CODONS_POS = ["TAA", "TAG", "TGA"]
    STOP_CODONS_NEG = ["TTA", "CTA", "TCA"]
    
    # Determine what our start and stop codons are
    currentStart = mrnaFeature_FastASeq_obj.seq[0:3].upper()
    currentStop = mrnaFeature_FastASeq_obj.seq[-3:].upper()
    
    newCoords = [None, None]
    
    # Crawl up the genome sequence looking for a way to extend the ORF to an accepted start as determined above
    if mrnaFeature.strand == '+' and currentStart not in START_CODONS_POS:
        # Crawl back looking for the first stop codon
        genomeSeq = genomeFASTA_obj[mrnaFeature.contig].seq[0:mrnaFeature.start-1] # start is 1-based so we -1 to counter that
        for i in range(len(genomeSeq)-1, -1, -3):
            codon = genomeSeq[i-2:i+1]
            if codon.upper() in STOP_CODONS_POS:
                break
        
        # Make sure we haven't gone past the start of the contig
        if genomeSeq == "": # Handles scenario where the gene model starts at the first base of the contig
            i = 0 # i would otherwise be None or the previous iteration's value
        else:
            i = i - 2 # This walks our coordinate value back to the start of the codon (Atg) since our index currently corresponds to (atG)
        
        # Crawl back up from the stop position looking for the first accepted start codon
        accepted = None
        for x in range(i+3, len(genomeSeq), 3): # +3 to look at the next, non-stop codon
            codon = genomeSeq[x:x+3]
            if codon.upper() in START_CODONS_POS:
                accepted = x + 1 # Note that this x represents the distance in from the stop codon boundary; +1 to reconvert this to 1-based
                break
        
        # Update this in our coords value
        if accepted != None:
            newCoords = [accepted, newCoords[1]]
    
    elif mrnaFeature.strand == '-' and currentStart not in START_CODONS_POS:
        # Crawl up looking for the first stop codon
        genomeSeq = genomeFASTA_obj[mrnaFeature.contig].seq[mrnaFeature.end:] # end is 1-based; we want just after it, so accepting it as-is is correct
        for i in range(0, len(genomeSeq), 3):
            codon = genomeSeq[i:i+3]
            if codon.upper() in STOP_CODONS_NEG: # after this, i will equal the distance from the currentStart to the stop codon boundary
                break # After this, i will equal the distance from the currentStart to the stop codon boundary
        
        # Make sure we haven't gone past the end of the contig
        if genomeSeq == "": # Handles scenario where the gene model starts at the last base of the contig
            i = 0 # i would otherwise be None or the previous iteration's value
        
        # Crawl back down from the stop position looking for the first current start or ATG
        accepted = None
        for x in range(i-1, -1, -3): # Here, we go from our stop codon boundary (which we went out/to the right to derive earlier) back in to the sequence/to the left
            codon = genomeSeq[x-2:x+1]
            if codon.upper() in START_CODONS_NEG:
                accepted = x + mrnaFeature.end + 1 # Note that this x represents the distance out to the stop codon boundary; +startCoord to +1 to reconvert this to 1-based
                break
        
        # Update this in our coords value
        if accepted != None:
            newCoords = [newCoords[0], accepted]
    
    # Crawl down the genome sequence looking for a way to extend the ORF to an accepted stop as determined above
    if mrnaFeature.strand == '+' and currentStop not in STOP_CODONS_POS:
        # Trim off excess from the CDS to make sure we're in frame
        endCoord = mrnaFeature.end - (len(mrnaFeature_FastASeq_obj.seq) % 3)
        
        # Perform the coord walk
        genomeSeq = genomeFASTA_obj[mrnaFeature.contig].seq[endCoord:] # endCoord is 1-based; we want just after it, so accepting it as-is is correct
        for i in range(0, len(genomeSeq), 3):
            codon = genomeSeq[i:i+3]
            if codon.upper() in STOP_CODONS_POS:
                break
        
        if i > MAX_EXTENSION:
            pass
        else:
            i = endCoord + i + 2 + 1 # +2 to go to the end of the stop codon; +1 to make it 1-based
            newCoords = [newCoords[0], i]
    elif mrnaFeature.strand == '-' and currentStop not in STOP_CODONS_POS:
        # Trim off excess from the CDS to make sure we're in frame
        endCoord = mrnaFeature.start + (len(mrnaFeature_FastASeq_obj.seq) % 3)
        
        # Perform the coord walk
        genomeSeq = genomeFASTA_obj[mrnaFeature.contig].seq[0:endCoord-1] # endCoord is 1-based so we -1 to counter that
        for i in range(len(genomeSeq)-1, -1, -3):
            codon = genomeSeq[i-2:i+1]
            if codon.upper() in STOP_CODONS_NEG:
                    break
        
        if (endCoord - i) > MAX_EXTENSION:
            pass
        else:
            i = i - 2 + 1 # -2 to go to the start of the codon; +1 to make it 1-based
            newCoords = [i, newCoords[1]]
    
    # Update any coordinates as needed
    stophere
    ## TBD
    # coord_cds_region_update(coords, startChange, stopChange, orientation)
    
    raise NotImplementedError("AAAA")
    #return coords, cds, startCodonsPos                              # We want to return the start codons since we'll use them later

def coord_cds_region_update(coords, startChange, stopChange, orientation):
    # Part 1: Cull exons that aren't coding and figure out how far we are chopping into coding exons
    origStartChange = startChange
    origStopChange = stopChange
    startExonLen = 0
    stopExonLen = 0
    microExonSize = -3      # This value is an arbitrary measure where, if a terminal exon is less than this size, we consider it 'fake' and delete it
    for i in range(2):
        while True:
            if i == 0:                      # The GFF3 is expected to be formatted such that + features are listed lowest coord position -> highest, whereas - features are listed from highest -> lowest
                exon = coords[0]        # This is done by PASA and by the exonerate GFF3 parsing system of this code, and it basically means that our first listed exon is always our starting exon
            else:
                exon = coords[-1]
            # Extract details
            rightCoord = int(exon[1])
            leftCoord = int(exon[0])
            exonLen = rightCoord - leftCoord + 1
            # Update our change values
            if i == 0:
                startChange -= exonLen          # This helps us to keep track of how much we need to reduce startChange
                if startChange > 0:             # when we begin chopping into the first exon - if the original first exon 
                    del coords[0]           # is the one we chop, we end up with reduction value == 0
                    startExonLen += exonLen
                # Handle microexons at gene terminal
                elif startChange > microExonSize:
                    del coords[0]
                    startExonLen += exonLen
                else:
                    break
            else:
                stopChange -= exonLen
                if stopChange > 0:
                    del coords[-1]
                    stopExonLen += exonLen  # We hold onto exon lengths so we can calculate how much we chop into the new start exon
                # Handle microexons at gene terminal
                elif stopChange > microExonSize:
                    del coords[-1]
                    stopExonLen += exonLen
                else:                           # by calculating stopChange - stopExonLen... if we didn't remove an exon, stopExonLen == 0
                    break
    origStartChange -= startExonLen
    origStopChange -= stopExonLen
    # Step 2: Using the chopping lengths derived in part 1, update our coordinates
    for i in range(len(coords)):
        splitCoord = coords[i]
        if i == 0:
            if orientation == '+':
                splitCoord[0] = int(splitCoord[0]) + origStartChange
            else:
                splitCoord[1] = int(splitCoord[1]) - origStartChange
        if i == len(coords) - 1:
            if orientation == '+':
                splitCoord[1] = int(splitCoord[1]) - origStopChange
            else:
                splitCoord[0] = int(splitCoord[0]) + origStopChange
        coords[i] = splitCoord
    return coords

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
        resultsDict = exonerateSearcher.run_exonerate() # returns a list of GFF3 features
        
        # Save result as pickle & clean up temp file
        with open(exoneratePickleFile, "wb") as fileOut:
            pickle.dump([resultsDict, nonequalIDs, fixesDict], fileOut)
        os.unlink(tmpFileName)
    else:
        # Load in the exonerate pickle
        with open(exoneratePickleFile, "rb") as fileIn:
            resultsDict, nonequalIDs, fixesDict = pickle.load(fileIn)
    
    # Filter resultsDict to only have relevant results
    MIN_IDENTITY, MIN_SIMILARITY = 98.0, 98.0
    resultsDict = ZS_AlignIO.Exonerate.filter_exonerate_resultsDict(
        resultsDict, num_hits=5, identity=MIN_IDENTITY, similarity=MIN_SIMILARITY
    )
    
    # Attempt to fix genes
    exonerateResultsDict = resultsDict
    fix_nonequivalent_features(nonequalIDs, resultsDict, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, args.isProtein)
    
    
    # Handle simple fixes
    ## TBD
    
    # mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(GFF3_obj[geneID])
    # cds_FastASeq_obj, cds_featureType, cds_startingFrame = GFF3_obj.retrieve_sequence_from_FASTA(genomeFASTA_obj, mrnaFeature.ID, "CDS")
    
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
