#! python3
# gff3_coordinate_fixer.py
# Script to receive a GFF3 and a released nucleotide sequence
# file and will attempt to correct any GFF3 features where
# the coordinates don't match the real feature.

import os, argparse, sys, hashlib, pickle
from copy import deepcopy

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_GFF3IO, ZS_SeqIO, ZS_AlignIO

START_CODONS_POS = ["ATG"] #, "CTG", "TTG"]
START_CODONS_NEG = ["CAT"] #, "CAG", "CAA"]
STOP_CODONS_POS = ["TAA", "TAG", "TGA"]
STOP_CODONS_NEG = ["TTA", "CTA", "TCA"]

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
    # Handle optional parameters
    if args.gmapDir != None or args.transcriptFastaFile != None:
        assert args.transcriptFastaFile != None and args.gmapDir != None, \
            "Must receive both transcriptFASTA_obj and gmapDir if you provide one of these!"
        if not os.path.isdir(args.gmapDir):
            print('I am unable to locate the GMAP binaries directory (' + args.gmapDir + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        if not os.path.isfile(os.path.join(args.gmapDir, "gmap")):
            print('I am unable to locate the gmap executable (expected to be at ' + os.path.join(args.gmapDir, "gmap") + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        if not os.path.isfile(os.path.join(args.gmapDir, "gmap_build")):
            print('I am unable to locate the gmap_build executable (expected to be at ' + os.path.join(args.gmapDir, "gmap_build") + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        if args.transcriptFastaFile != None:
            if not os.path.isfile(args.transcriptFastaFile):
                print('I am unable to locate the transcript FASTA file (' + args.transcriptFastaFile + ')')
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
        problemIDs -- a list containing ID values corresponding to
                      sequences that don't match between the GFF3 and
                      FASTA objects at all.
        fixesDict -- a dictionary with structure like:
                     {
                         'seqID1': [strand, frame],
                         'seqID2': [strand, frame],
                         ...
                     }
    '''
    problemIDs = []
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
        gff3_FastASeq_obj, _, startingFrame = GFF3_obj.retrieve_sequence_from_FASTA(
            genomeFASTA_obj, feature.ID, "CDS"
        )
        gff3ProtSequence, _, _ = gff3_FastASeq_obj.get_translation(
            findBestFrame=False,
            strand=1,
            frame=int(startingFrame) if startingFrame.isdigit() else 0
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
                problemIDs.append(FastASeq_obj.id)
            else:
                fixesDict[FastASeq_obj.id] = [alt_strand, alt_frame]
    return problemIDs, fixesDict

def fix_features_with_exonerate_and_gmap(problemIDs, exonerateResultsDict, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, isProtein=False, transcriptFASTA_obj=None, gmapRunner=None):
    '''
    Intended to follow on from find_nonequivalent_features(), this function will
    instead aim to correct the coordinates of any GFF3 features that don't match
    their respective CDS gene model. Changes will be directly made to the input
    GFF3_obj.
    
    Parameters:
        problemIDs -- a list containing ID values corresponding to
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
        gmapRunner -- OPTIONAL; a ZS_AlignIO.GMAP object set with a dummy query file
                      and the appropriate genome FASTA set as the target.
    Returns:
        unfixedIDs -- a list indicating the sequences that we failed to fix using
                      this function
        fixedIDs -- a list indicating the sequences that we successfully fixed using
                      this function
        
    '''
    LENGTH_PCT_DIFF_ALLOWANCE = 0.10
    SIMILARITY_MINUMUM_ALLOWANCE = 90.0
    
    if transcriptFASTA_obj != None or gmapRunner != None:
        assert transcriptFASTA_obj != None and gmapRunner != None, \
            "fix_features_with_exonerate_and_gmap needs both transcriptFASTA_obj AND gmapRunner; providing only one doesn't work"
    
    unfixedIDs = []
    fixedIDs = []
    for problemSeqID in problemIDs:
        cdsFastASeq_obj = cdsFASTA_obj[problemSeqID]
        feature = GFF3_obj[problemSeqID]
        
        # Get translation of FASTA sequence if relevant
        if not isProtein:
            cdsProtSequence, _, _ = cdsFastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
        else:
            cdsProtSequence = cdsFastASeq_obj.seq
        
        # Get exonerate results if applicable
        if problemSeqID in exonerateResultsDict:
            exonerateResults = exonerateResultsDict[problemSeqID]
        else:
            exonerateResults = None
        
        # Get GMAP results if applicable
        if gmapRunner != None:
            gmapRunner.query = transcriptFASTA_obj[problemSeqID]
            gmapResults = gmapRunner.run_gmap()
            gmapRunner.clean(query=True)
        else:
            gmapResults = None
        
        # Skip if we have no results, or merge them otherwise
        if exonerateResults == None and gmapResults == None:
            continue
        elif exonerateResults != None and gmapResults != None:
            searchResults = exonerateResults + gmapResults.types["gene"]
        elif exonerateResults == None:
            searchResults = gmapResults.types["gene"]
        else:
            searchResults = exonerateResults
        
        # Find the best result in the closest proximity to the original gene annotation
        bestResults = [] # [feature, distance, lengthDifference, similarity]
        for geneFeature in searchResults:
            mrnaFeature = geneFeature.mRNA[0]
            
            # Get the distance between the genes (if any)
            if mrnaFeature.end >= feature.start and feature.end >= mrnaFeature.start:
                distance = 0 # if they overlap, there's essentially no distance between them
            else:
                distance = min(abs(mrnaFeature.end - feature.start), abs(feature.end - mrnaFeature.start))
            
            # Extend the gene to its start codon, stop codon boundaries
            mrnaFeature_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
                genomeFASTA_obj, mrnaFeature, "CDS"
            )
            mrnaFeature = feature_cds_extension_maximal(mrnaFeature, mrnaFeature_FastASeq_obj, genomeFASTA_obj, MAX_EXTENSION=30)
            
            # Validate that start and stop codons were found appropriately
            mrnaFeature_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
                genomeFASTA_obj, mrnaFeature, "CDS"
            )
            startCodon = mrnaFeature_FastASeq_obj.seq[0:3].upper()
            stopCodon = mrnaFeature_FastASeq_obj.seq[-3:].upper()
            if startCodon not in START_CODONS_POS or stopCodon not in STOP_CODONS_POS:
                continue
            
            # Validate that the difference in gene length isn't crazy
            mrnaFeature_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
                genomeFASTA_obj, mrnaFeature, "CDS"
            )
            lengthDifference = abs(len(cdsProtSequence) - (len(mrnaFeature_FastASeq_obj.seq) / 3))
            if lengthDifference > (len(cdsProtSequence) * LENGTH_PCT_DIFF_ALLOWANCE):
                continue
            
            # Validate that the similarity score is good
            similarity = float(mrnaFeature.similarity) if hasattr(mrnaFeature, "similarity") \
                else (float(mrnaFeature.identity) * float(mrnaFeature.coverage)) / 100
            if similarity < SIMILARITY_MINUMUM_ALLOWANCE:
                continue
            
            # Finally, store the result
            bestResults.append([geneFeature, distance, lengthDifference, similarity])
        bestResults.sort(key = lambda x: (-x[3], x[1], x[2], -(x[0].end - x[0].start))) # last sort will prioritise UTR prediction from GMAP
        
        # If there's no good results, note it and continue
        if bestResults == []:
            unfixedIDs.append(problemSeqID)
            continue
        # If there is a good result, update this model
        else:
            del GFF3_obj[problemSeqID]
            GFF3_obj.add_feature(bestResults[0][0])
            fixedIDs.append(problemSeqID)
    return unfixedIDs, fixedIDs

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
    Returns:
        mrnaFeature -- the same input object, modified as per the actions of this function
    '''
    # Determine what our start and stop codons are
    currentStart = mrnaFeature_FastASeq_obj.seq[0:3].upper()
    currentStop = mrnaFeature_FastASeq_obj.seq[-3:].upper()
    
    leftBorder = min([cdsFeature.start for cdsFeature in mrnaFeature.CDS])
    rightBorder = max([cdsFeature.end for cdsFeature in mrnaFeature.CDS])
    
    newCoords = [None, None]
    
    # Crawl up the genome sequence looking for a way to extend the ORF to an accepted start as determined above
    if mrnaFeature.strand == '+' and currentStart not in START_CODONS_POS:
        # Crawl back looking for the first stop codon
        genomeSeq = genomeFASTA_obj[mrnaFeature.contig].seq[0:leftBorder-1] # start is 1-based so we -1 to counter that
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
            newCoords = [accepted, newCoords[1]] # [leftBorder, rightBorder]
    
    elif mrnaFeature.strand == '-' and currentStart not in START_CODONS_POS:
        # Crawl up looking for the first stop codon
        genomeSeq = genomeFASTA_obj[mrnaFeature.contig].seq[rightBorder:] # end is 1-based; we want just after it, so accepting it as-is is correct
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
                accepted = x + rightBorder + 1 # Note that this x represents the distance out to the stop codon boundary; +startCoord to +1 to reconvert this to 1-based
                break
        
        # Update this in our coords value
        if accepted != None:
            newCoords = [newCoords[0], accepted] # [leftBorder, rightBorder]
    
    # Crawl down the genome sequence looking for a way to extend the ORF to an accepted stop as determined above
    if mrnaFeature.strand == '+' and currentStop not in STOP_CODONS_POS:
        # Trim off excess from the CDS to make sure we're in frame
        endCoord = rightBorder - (len(mrnaFeature_FastASeq_obj.seq) % 3)
        
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
            newCoords = [newCoords[0], i] # [leftBorder, rightBorder]
    elif mrnaFeature.strand == '-' and currentStop not in STOP_CODONS_POS:
        # Trim off excess from the CDS to make sure we're in frame
        endCoord = leftBorder + (len(mrnaFeature_FastASeq_obj.seq) % 3)
        
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
            newCoords = [i, newCoords[1]] # [leftBorder, rightBorder]
    
    # Update coordinates
    newLeftExonBorder = newCoords[0] if newCoords[0] != None and newCoords[0] < mrnaFeature.start else mrnaFeature.start
    newRightExonBorder = newCoords[1] if newCoords[1] != None and newCoords[1] > mrnaFeature.end else mrnaFeature.end
    
    newLeftCdsBorder = newCoords[0] if newCoords[0] != None and newCoords[0] < leftBorder else leftBorder
    newRightCdsBorder = newCoords[1] if newCoords[1] != None and newCoords[1] > rightBorder else rightBorder
    
    mrnaFeature = feature_coords_update(
        mrnaFeature,
        newLeftExonBorder,
        newRightExonBorder,
        newLeftCdsBorder,
        newRightCdsBorder
    )
    
    return mrnaFeature

def feature_coords_update(feature, newLeftExonBorder, newRightExonBorder, newLeftCdsBorder, newRightCdsBorder):
    '''
    Function to receive a feature and new start / stop coordinates and update the feature
    to cull any subfeatures/modify them to account for this change. Subfeature type must
    be specified, and is typically expected to be "exon" or "CDS"
    
    Parameters:
        mrnaFeature -- a ZS_GFF3IO.Feature object of any feature containing .exon values
        newLeftExonBorder -- an integer indicating the new start / left boundary to the feature's exons
        newLeftExonBorder -- an integer indicating the new end / right boundary to the feature's exons
        newLeftCdsBorder -- an integer indicating the new start / left boundary to the feature's CDS
        newRightCdsBorder -- an integer indicating the new end / right boundary to the feature's CDS
    '''
    assert isinstance(newLeftExonBorder, int) and isinstance(newRightExonBorder, int) \
        and isinstance(newLeftCdsBorder, int) and isinstance(newRightCdsBorder, int), \
        "new border values must all be integers"
    
    assert hasattr(feature, "CDS") and hasattr(feature, "exon"), \
        "feature object lacks CDS and exon values; updating its coords is not relevant or possible"
    
    # Make sure the parameters are sensible in the context of this function
    "Someone might mistake the left/right-ness of our start/stop idea and go for a sense/antisense based approach"
    newLeftExonBorder, newRightExonBorder = min(newLeftExonBorder, newRightExonBorder), max(newLeftExonBorder, newRightExonBorder)
    newLeftCdsBorder, newRightCdsBorder = min(newLeftCdsBorder, newRightCdsBorder), max(newLeftCdsBorder, newRightCdsBorder)
    
    # Make sure feature values are sorted ascendingly for this function
    feature.exon.sort(key = lambda x: x.start)
    feature.CDS.sort(key = lambda x: x.start)
    
    # Modify the start
    ## Exons
    exonsToDrop = 0
    for exonFeature in feature.exon:
        if exonFeature.end < newLeftExonBorder: # drop this exon
            exonsToDrop += 1
        else: # modify this exon
            exonFeature.start = newLeftExonBorder
            exonFeature.coords = [newLeftExonBorder, exonFeature.end]
            feature.start = newLeftExonBorder # the feature .start value should reflect exon, not CDS
            feature.coords = [newLeftExonBorder, feature.end] # likewise, .coords reflect exons
            break
    feature.exon = feature.exon[exonsToDrop:]
    ## CDS
    cdsToDrop = 0
    for cdsFeature in feature.CDS:
        if cdsFeature.end < newLeftCdsBorder: # drop this CDS
            cdsToDrop += 1
        else: # modify this CDS
            cdsFeature.start = newLeftCdsBorder
            cdsFeature.coords = [newLeftCdsBorder, cdsFeature.end]
            break
    feature.CDS = feature.CDS[cdsToDrop:]
    
    # Modify the end
    ## Exons
    exonsToDrop = 0
    for exonFeature in feature.exon[::-1]:
        if exonFeature.start > newRightExonBorder: # drop this exon
            exonsToDrop += 1
        else: # modify this exon
            exonFeature.end = newRightExonBorder
            exonFeature.coords = [exonFeature.start, newRightExonBorder]
            feature.end = newRightExonBorder
            feature.coords = [feature.start, newRightExonBorder]
            break
    if exonsToDrop != 0:
        feature.exon = feature.exon[:-exonsToDrop]
    ## CDS
    cdsToDrop = 0
    for cdsFeature in feature.CDS[::-1]:
        if cdsFeature.start > newRightCdsBorder: # drop this CDS
            cdsToDrop += 1
        else: # modify this CDS
            cdsFeature.end = newRightCdsBorder
            cdsFeature.coords = [cdsFeature.start, newRightCdsBorder]
            break
    if cdsToDrop != 0:
        feature.exon = feature.exon[:-cdsToDrop]
    
    # Address situations where the range is contained within an intron
    if feature.exon == []:
        newExonFeature = ZS_GFF3IO.Feature()
        newExonFeature.add_attributes({
            "ID": f"{feature.ID}.exon.1", "Parent": feature.ID,
            "contig": feature.contig, "source": feature.source, "type": "exon",
            "start": newLeftExonBorder, "end": newRightExonBorder, "coords": [newLeftExonBorder, newRightExonBorder],
            "score": ".", "strand": feature.strand, "frame": "."
        })
        feature.exon.append(newExonFeature)
        feature.coords = [newLeftExonBorder, newRightExonBorder]
        
        # Update the main feature values
        "We only do this here since, entering this if statement must mean we also entered the one above"
        feature.start = newLeftExonBorder
        feature.end = newRightExonBorder
        feature.coords = [newLeftExonBorder, newRightExonBorder]
    
    # Sort feature values as they should be in general GFF3 convention again
    "+ve stranded features will be sorted appropriately if they're in ascending order"
    if feature.strand == "-":
        feature.exon.sort(key = lambda x: -x.end)
        feature.CDS.sort(key = lambda x: -x.end)
    
    return feature

def fix_genes_by_renaming(problemIDs, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, isProtein=False):
    '''
    This function receives a GFF3 and a FASTA that is meant to
    correspond to the sequences extracted from that GFF3. The goal
    is to find cases where the sequences differ, and to return
    these feature/sequence IDs as a list for further work to be done on them.
    
    It will also resolve easy to solve issues pertaining to strand and frame
    annotations. If either appear to be wrong, this function will identify it
    easily! It's up to you to change that in the GFF3_obj though.
    
    Parameters:
        problemIDs -- a list containing ID values corresponding to
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
    Returns:
        unfixedIDs -- a list indicating the sequences that we failed to fix using
                      this function
        fixedIDs -- a list indicating the sequences that we successfully fixed using
                      this function
    '''
    unfixedIDs = []
    fixedIDs = []
    
    fastaSeqDict = {}
    for problemSeqID in problemIDs:
        cdsFastASeq_obj = cdsFASTA_obj[problemSeqID]
        feature = GFF3_obj[problemSeqID]
        
        # Get translation of FASTA sequence if relevant
        if not isProtein:
            cdsProtSequence, _, _ = cdsFastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
        else:
            cdsProtSequence = cdsFastASeq_obj.seq
        
        # Get translation of GFF3 feature
        gff3_FastASeq_obj, _, startingFrame = GFF3_obj.retrieve_sequence_from_FASTA(
            genomeFASTA_obj, feature.ID, "CDS"
        )
        gff3ProtSequence, _, _ = gff3_FastASeq_obj.get_translation(
            findBestFrame=False,
            strand=1,
            frame=int(startingFrame) if startingFrame.isdigit() else 0
        )
        
        # If this translation matches, that means we must've swapped it already, and it fixed things
        if gff3ProtSequence.rstrip("*") == cdsProtSequence.rstrip("*"):
            fixedIDs.append(problemSeqID)
            continue
        
        # Get translated sequences of all features of the same type as the problem
        for typeFeature in GFF3_obj.types[feature.type]:
            if typeFeature.ID not in fastaSeqDict:
                type_FastASeq_obj, _, startingFrame = GFF3_obj.retrieve_sequence_from_FASTA(
                    genomeFASTA_obj, typeFeature.ID, "CDS"
                )
                typeProtSequence, _, _ = type_FastASeq_obj.get_translation(
                    findBestFrame=False,
                    strand=1,
                    frame=int(startingFrame) if startingFrame.isdigit() else 0
                )
                fastaSeqDict[typeFeature.ID] = typeProtSequence
        
        # Look for a match to this sequence
        matchIDs = [seqID for seqID, protSeq in fastaSeqDict.items() if protSeq.rstrip("*") == cdsProtSequence.rstrip("*")]
        
        # Filter 1: match(es) to just a single good one
        filter1IDs = []
        for matchID in matchIDs:
            # Match must not exist in CDS file, or it must be a problem ID
            if matchID in cdsFASTA_obj and not matchID in problemIDs:
                continue
            filter1IDs.append(matchID)
        
        # Filter 2: If multiple matches, select the one on the same contig
        if len(filter1IDs) > 1:
            filter2IDs = []
            for matchID in filter1IDs:
                if GFF3_obj[matchID].contig == feature.contig:
                    filter2IDs.append(matchID)
        else:
            filter2IDs = filter1IDs
        
        # Filter 3: Only allow a single unambiguous match
        if len(filter2IDs) != 1:
            unfixedIDs.append(problemSeqID)
            continue
        matchID = filter2IDs[0]
        
        # If we have a good match and it is a problem sequence...
        if matchID in problemIDs:
            "It being a problem sequence means it IS in the CDS file"
            # ... and if it's already been fixed, we won't touch it
            if matchID in fixedIDs:
                unfixedIDs.append(problemSeqID)
                continue
            
            # ... but, if it hasn't been fixed, let's just swap them
            "This might fix both models, but if it doesn't we'll get to it later"
            GFF3_obj[problemSeqID].update_id(problemSeqID, "tmp_problem_seq_id", GFF3_obj)
            GFF3_obj[matchID].update_id(matchID, "tmp_match_seq_id", GFF3_obj)
            
            GFF3_obj["tmp_problem_seq_id"].update_id("tmp_problem_seq_id", matchID, GFF3_obj)
            GFF3_obj["tmp_match_seq_id"].update_id("tmp_match_seq_id", problemSeqID, GFF3_obj)
        
        # If we have a good match and it isn't in the CDS file...
        else:
            "If we get here, we know it ISN'T in the CDS file"
            # ... just replace the model
            del GFF3_obj[problemSeqID]
            GFF3_obj[matchID].update_id(matchID, problemSeqID, GFF3_obj)
        
        # Note it as fixed and continue
        fixedIDs.append(problemSeqID)
    
    return unfixedIDs, fixedIDs

def fix_genes_by_sliding(problemIDs, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, isProtein=False, slideLength=500):
    '''
    This function attempts to reconcile issues seen in the GFF3 annotation when compared
    to extracted CDS sequences by sliding the gene models along the contig looking
    for a match. If we found a match, we'll save it in the GFF3 object directly.
    Otherwise, we'll return a list of the sequences we couldn't fix like this.
    
    Parameters:
        slideLength -- an integer limiting how far we want to slide the feature
                       looking for the fixed sequence
    Returns:
        unfixedIDs -- a list indicating the sequences that we failed to fix using
                      this function
        fixedIDs -- a list indicating the sequences that we successfully fixed using
                      this function
    '''
    unfixedIDs = []
    fixedIDs = []
    for problemSeqID in problemIDs:
        cdsFastASeq_obj = cdsFASTA_obj[problemSeqID]
        feature = GFF3_obj[problemSeqID]
        
        # Get translation of FASTA sequence if relevant
        if not isProtein:
            cdsProtSequence, _, _ = cdsFastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
        else:
            cdsProtSequence = cdsFastASeq_obj.seq
        
        # Slide the gene left
        foundFix = False
        leftSlideFeature = deepcopy(feature) # make a backup so changes aren't permanent
        for i in range(min(slideLength, len(genomeFASTA_obj[feature.contig].seq) - feature.start)):
            slide_feature(leftSlideFeature, 1, "left")
            leftSlide_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
                genomeFASTA_obj, leftSlideFeature, "CDS"
            )
            leftSlideProtSequence, _, _ = leftSlide_FastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
            
            # If we got a hit, store it and move on
            if leftSlideProtSequence.rstrip("*") == cdsProtSequence:
                foundFix = True
                slide_feature(GFF3_obj[problemSeqID], i+1, "left")
                
                # Reset the parent gene feature with the new mRNA feature (if applicable)
                if hasattr(leftSlideFeature, "Parent"): # this means we are handling a mRNA feature
                    slide_feature(GFF3_obj[rightSlideFeature.Parent], i+1, "left")
                break
        if foundFix is True:
            fixedIDs.append(problemSeqID)
            continue
        
        # Slide the gene right
        rightSlideFeature = deepcopy(feature) # make a backup so changes aren't permanent
        for i in range(min(slideLength, len(genomeFASTA_obj[feature.contig].seq) - feature.end)):
            slide_feature(rightSlideFeature, 1, "right")
            rightSlide_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
                genomeFASTA_obj, rightSlideFeature, "CDS"
            )
            rightSlideProtSequence, _, _ = rightSlide_FastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
            
            # If we got a hit, fix the original sequence and move on
            if rightSlideProtSequence.rstrip("*") == cdsProtSequence:
                foundFix = True
                slide_feature(GFF3_obj[problemSeqID], i+1, "right")
                
                # Reset the parent gene feature with the new mRNA feature (if applicable)
                if hasattr(rightSlideFeature, "Parent"):
                    slide_feature(GFF3_obj[rightSlideFeature.Parent], i+1, "right")
                break
        if foundFix is True:
            fixedIDs.append(problemSeqID)
            continue
        
        unfixedIDs.append(problemSeqID)
    return unfixedIDs, fixedIDs

def slide_feature(mrnaFeature, slideLength, direction):
    '''
    Simply put, will slide all the coordinates of a feature to the left or right
    a set distance. That means each coordinate value will be +slideLength (right)
    or -slideLength (left).
    
    This function isn't safe with respect to sliding things out of bounds. That's
    up to you to make sure of.
    
    Parameters:
        mrnaFeature -- a ZS_GFF3IO.Feature object indicating an mRNA feature or,
                       just any feature with both CDS and exon subfeatures
        slideLength -- an integer limiting how far we want to slide the feature
                       looking for the fixed sequence
        direction -- a string in the list of ["left", "right"] which will result
                     in us sliding the feature
    '''
    # assert hasattr(mrnaFeature, "CDS") and hasattr(mrnaFeature, "exon"), \
    #     "slide_feature will only accept a feature with both CDS and exon attributes"
    assert isinstance(slideLength, int) and slideLength > 0, \
        "slideLength must be a positive integer"
    assert direction.lower() in ["left", "right"], \
        "direction value must be 'left' or 'right'"
    
    def _slide_it_left(feature, slideLength):
        feature.start = feature.start - slideLength
        feature.end = feature.end - slideLength
        feature.coords = [feature.start, feature.end]
    
    def _slide_it_right(feature, slideLength):
        feature.start = feature.start + slideLength
        feature.end = feature.end + slideLength
        feature.coords = [feature.start, feature.end]
    
    if direction.lower() == "left":
        _slide_it_left(mrnaFeature, slideLength)
        try:
            for cdsFeature in mrnaFeature.CDS:
                _slide_it_left(cdsFeature, slideLength)
        except:
            pass # we mustn't have .CDS features
        try:
            for exonFeature in mrnaFeature.exon:
                _slide_it_left(exonFeature, slideLength)
        except:
            pass # we mustn't have .exon features
    else:
        _slide_it_right(mrnaFeature, slideLength)
        try:
            for cdsFeature in mrnaFeature.CDS:
                _slide_it_right(cdsFeature, slideLength)
        except:
            pass
        try:
            for exonFeature in mrnaFeature.exon:
                _slide_it_right(exonFeature, slideLength)
        except:
            pass

def feature_contig_update(feature, contig):
    '''
    Simply put, will update the contig attribute for all subfeatures.
    
    Parameters:
        mrnaFeature -- a ZS_GFF3IO.Feature object indicating an mRNA feature or,
                       just any feature with both CDS and exon subfeatures
        contig -- a string indicating what value should be set for .contig in this
                  and all children features
    '''
    feature.contig = contig
    for childFeature in feature.children:
        feature_contig_update(childFeature, contig)

def fix_genes_by_contig_checking(problemIDs, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, isProtein=False, slideLength=500):
    '''
    This function attempts to reconcile issues seen in the GFF3 annotation when compared
    to extracted CDS sequences by checking to see if the annotation is labeled on the
    wrong contig.
    
    Returns:
        unfixedIDs -- a list indicating the sequences that we failed to fix using
                      this function
        fixedIDs -- a list indicating the sequences that we successfully fixed using
                      this function
    '''
    unfixedIDs = []
    fixedIDs = []
    for problemSeqID in problemIDs:
        cdsFastASeq_obj = cdsFASTA_obj[problemSeqID]
        feature = GFF3_obj[problemSeqID]
        
        # Get translation of FASTA sequence if relevant
        if not isProtein:
            cdsProtSequence, _, _ = cdsFastASeq_obj.get_translation(
                findBestFrame=False,
                strand=1,
                frame=0
            )
        else:
            cdsProtSequence = cdsFastASeq_obj.seq
        
        # Make a copy of the feature and begin scrolling through all available contigs for a good hit
        feature = deepcopy(feature)
        foundFix = False
        for contig in genomeFASTA_obj.ids:
            if len(genomeFASTA_obj[contig].seq) < feature.end:
                continue
            
            # Get the sequence for this contig
            feature.contig = contig
            feature_FastASeq_obj, _, _ = GFF3_obj.retrieve_sequence_from_FASTA(
                genomeFASTA_obj, feature, "CDS"
            )
            featureProtSequence, _, _ = feature_FastASeq_obj.get_translation(
                findBestFrame=True,
                strand=1
            )
            
            # If we got a hit, update the original feature and move on
            if featureProtSequence.rstrip("*") == cdsProtSequence:
                foundFix = True
                feature_contig_update(GFF3_obj[problemSeqID], contig)
                
                # Reset the parent gene feature with the new mRNA feature (if applicable)
                if hasattr(feature, "Parent"): # this means we are handling a mRNA feature
                    feature_contig_update(GFF3_obj[feature.Parent], contig)
                break
        
        if foundFix is True:
            fixedIDs.append(problemSeqID)
            continue
        
        unfixedIDs.append(problemSeqID)
    return unfixedIDs, fixedIDs

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
    ARGS_TO_IGNORE = ["reportToFile", "outputFileName", "resume"] # these don't change the operation of the program
    hashes = []
    for key, value in args.__dict__.items():
        if key in ARGS_TO_IGNORE:
            continue
        hashes.append(hashlib.sha256(bytes(str(key) + str(value), 'utf-8')).hexdigest())
    overallHash = hashlib.sha256(bytes("".join(hashes), 'utf-8')).hexdigest()
    return overallHash[0:20]

def report_program_results(GFF3_obj, outputFileName, originalProblems, numFixedByRenaming,
                           numFixedByStrandChecking, numFixedBySliding,
                           numFixedByContigChecking, numFixedByReannotation,
                           problemIDs, isLastReport=False, reportToFile=False):
    '''
    Function to print out a formatted report detailing what this program has done.
    It will automatically make the decision of whether the program is meant to continue
    or not depending on how many problems remain.
    '''
    # Get report file name if necessary
    if reportToFile is True:
        reportFileName = ZS_AlignIO._tmp_file_name_gen("gff3_fixer_report", "txt")
    
    # Format report
    report = [
        f"GFF3 being fixed = {GFF3_obj.fileLocation}",
        f"Number of problem sequences identified = {len(originalProblems)}",
        f"Number fixed by changing their sequence ID = {numFixedByRenaming}",
        f"Number fixed by changing strand annotation = {numFixedByStrandChecking}",
        f"Number fixed by sliding along their contig = {numFixedBySliding}",
        f"Number fixed by moving to another contig = {numFixedByContigChecking}",
        f"Number fixed by reannotation (GMAP/exonerate) = {numFixedByReannotation}",
        f"Number of problem sequences remaining = {len(problemIDs)}",
        "#####",
        "Original problem sequence IDs = {0}".format("\n".join(originalProblems))
    ]
    if len(problemIDs) != 0:
        report.append("#####")
        report.append("Unfixed sequence IDs = {0}".format("\n".join(problemIDs)))
    
    # Handle no problem found scenario
    if len(originalProblems) == 0 and numFixedByStrandChecking == 0:
        print("No issues have been detected; program will exit now")
        if reportToFile is True:
            with open(reportFileName, "w") as fileOut:
                fileOut.write("\n".join(report))
        quit()
    
    # Handle all problems fixed scenario
    if len(problemIDs) == 0:
        GFF3_obj.write(outputFileName)
        print("\n".join(report))
        if reportToFile is True:
            with open(reportFileName, "w") as fileOut:
                fileOut.write("\n".join(report))
        
        print("Program exited successfully after finding no more problems to fix!")
        quit()
    
    # Handle isLastReport scenario
    if isLastReport is True:
        GFF3_obj.write(outputFileName)
        print("\n".join(report))
        if reportToFile is True:
            with open(reportFileName, "w") as fileOut:
                fileOut.write("\n".join(report))

def main():
    # User input
    usage = """%(prog)s accepts a GFF3 file and ...
    
    This script is intended to work with protein CDS files (giving a nucleotide will just
    have it translated to protein), with exonerate providing the approximate gene models.
    But, if you have transcripts, this script can also employ GMAP to get the alignments
    including UTRs, from which the CDS file can be used to identify the ORF from within the
    model. This is the best case scenario.
    
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
    p.add_argument("-o", dest="outputFileName", required=True,
                help="Specify the location to write the modified GFF3 file to")
    p.add_argument("-e", dest="exonerateExe", required=True,
                help="Specify the location of the exonerate executable file")
    # Opts
    p.add_argument("--gmapDir", dest="gmapDir", required=False,
                help="""Optionally, specify the location of the GMAP binary files
                if you're going to provide a transcript FASTA file""")
    p.add_argument("--transcriptFastaFile", dest="transcriptFastaFile", required=False,
                help="""Optionally, specify the location of the transcript
                nucleotide FASTA file if you want to run GMAP""")
    p.add_argument("--isProtein", dest="isProtein", required=False, action="store_true",
                help="Optionally specify that the cdsFastaFile contains protein sequences",
                default=False)
    p.add_argument("--resume", dest="resume", required=False, action="store_true",
                help="""Optionally provide this tag to save a pickle of the exonerate
                results which can be resumed from""",
                default=False)
    p.add_argument("--reportToFile", dest="reportToFile", required=False, action="store_true",
                help="""Optionally provide this tag to write reporting information to file
                rather than printing to terminal""",
                default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Figure out the hash and exonerate pickle file for this run
    argsHash = get_args_hash(args)
    exoneratePickleFile = os.path.join(
        os.path.dirname(args.outputFileName),
        f"exonerate_results_{argsHash}.pkl"
    )
    
    # Set variables to keep track of program performance
    numFixedByStrandChecking = 0
    numFixedBySliding = 0
    numFixedByContigChecking = 0
    numFixedByRenaming = 0
    numFixedByReannotation = 0
    
    # Parse GFF3
    GFF3_obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False) # non-strict parsing
    
    # Parse FASTA files
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.genomeFastaFile)
    cdsFASTA_obj = ZS_SeqIO.FASTA(args.cdsFastaFile)
    
    # Find the problem sequence IDs
    problemIDs, fixesDict = find_nonequivalent_features(GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, args.isProtein)
    originalProblems = problemIDs
    numFixedByStrandChecking = len(fixesDict)
    
    # Enact fixes to sequence strands
    for featureID, strand_frame_pair in fixesDict.items():
        strand, frame = strand_frame_pair
        GFF3_obj[featureID].strand = "+" if strand == 1 else "-"
    
    # If no problems existed, or no more exist after strand fixing, report that and exit program
    report_program_results(GFF3_obj, args.outputFileName, originalProblems, numFixedByRenaming,
                           numFixedByStrandChecking, numFixedBySliding,
                           numFixedByContigChecking, numFixedByReannotation,
                           problemIDs, isLastReport=False, reportToFile=args.reportToFile)
    
    # Try to fix genes by renaming the model
    problemIDs, fixedIDs = fix_genes_by_renaming(problemIDs, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, args.isProtein)
    numFixedByRenaming = len(fixedIDs)
    report_program_results(GFF3_obj, args.outputFileName, originalProblems, numFixedByRenaming,
                           numFixedByStrandChecking, numFixedBySliding,
                           numFixedByContigChecking, numFixedByReannotation,
                           problemIDs, isLastReport=False, reportToFile=args.reportToFile)
    
    # Try to fix genes by sliding existing models up/down their contig
    problemIDs, fixedIDs = fix_genes_by_sliding(problemIDs, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, args.isProtein)
    numFixedBySliding = len(fixedIDs)
    report_program_results(GFF3_obj, args.outputFileName, originalProblems, numFixedByRenaming,
                           numFixedByStrandChecking, numFixedBySliding,
                           numFixedByContigChecking, numFixedByReannotation,
                           problemIDs, isLastReport=False, reportToFile=args.reportToFile)
    
    # Try to fix genes by checking if their contig is misannotated
    problemIDs, fixedIDs = fix_genes_by_contig_checking(problemIDs, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, args.isProtein)
    numFixedByContigChecking = len(fixedIDs)
    report_program_results(GFF3_obj, args.outputFileName, originalProblems, numFixedByRenaming,
                           numFixedByStrandChecking, numFixedBySliding,
                           numFixedByContigChecking, numFixedByReannotation,
                           problemIDs, isLastReport=False, reportToFile=args.reportToFile)
    
    # Get exonerate results to help with fixing any remaining IDs
    if args.resume is False or not os.path.isfile(exoneratePickleFile):
        # Get a FASTA file for the sequences we need to search against the genome
        problem_FASTA_obj = ZS_SeqIO.FASTA(None)
        for seqID in problemIDs:
            problem_FASTA_obj.add(cdsFASTA_obj[seqID])
        
        tmpFileName = ZS_AlignIO._tmp_file_name_gen("tmpExonerateQuery", "fasta")
        problem_FASTA_obj.write(tmpFileName)
        
        # Perform exonerate search
        exonerateSearcher = ZS_AlignIO.Exonerate(args.exonerateExe, tmpFileName, args.genomeFastaFile)
        exonerateSearcher.model = "protein2genome"
        resultsDict = exonerateSearcher.run_exonerate() # returns a list of GFF3 features
        
        # Save result as pickle & clean up temp file
        with open(exoneratePickleFile, "wb") as fileOut:
            pickle.dump(resultsDict, fileOut)
        os.unlink(tmpFileName)
    else:
        # Load in the exonerate pickle
        with open(exoneratePickleFile, "rb") as fileIn:
            resultsDict = pickle.load(fileIn)
    
    # Ensure that the GMAP database exists and we have an object to run it with if applicable
    if args.gmapDir != None:
        gmapRunner = ZS_AlignIO.GMAP(args.gmapDir, "ATCG", args.genomeFastaFile) # ATCG as placeholder query
        gmapRunner.gmap_build()
        gmapRunner.clean(query=True)
    else:
        gmapRunner = None
    
    # Filter resultsDict to only have relevant results
    MIN_IDENTITY, MIN_SIMILARITY = 95.0, 95.0
    resultsDict = ZS_AlignIO.Exonerate.filter_exonerate_resultsDict(
        resultsDict, num_hits=5, identity=MIN_IDENTITY, similarity=MIN_SIMILARITY
    )
    
    # Attempt to fix remaining genes via exonerate/GMAP search
    transcriptFASTA_obj = ZS_SeqIO.FASTA(args.transcriptFastaFile) \
        if args.transcriptFastaFile != None else None
    problemIDs, fixedIDs = fix_features_with_exonerate_and_gmap(problemIDs, resultsDict, GFF3_obj, cdsFASTA_obj, genomeFASTA_obj, args.isProtein, transcriptFASTA_obj, gmapRunner)
    numFixedByReannotation = len(fixedIDs)
    report_program_results(GFF3_obj, args.outputFileName, originalProblems, numFixedByRenaming,
                           numFixedByStrandChecking, numFixedBySliding,
                           numFixedByContigChecking, numFixedByReannotation,
                           problemIDs, isLastReport=True, reportToFile=args.reportToFile)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
