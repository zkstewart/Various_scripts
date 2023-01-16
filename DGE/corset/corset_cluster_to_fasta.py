#! python3
# corset_cluster_to_fasta.py
# Script to parse the corset results-clusters.txt file and the
# transcriptome FASTA, producing an output file containing
# representative sequences from each cluster

import os, argparse, math
from Bio import SeqIO

def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.corsetClustersFile):
        print(f"I was unable to locate the Corset clusters file at '{args.corsetClustersFile}'")
        print("Make sure you specified the right file name and location, then try again.")
        quit()
    if not os.path.isfile(args.fastaFile):
        print(f"I was unable to locate the transcriptome FASTA file at '{args.fastaFile}'")
        print("Make sure you specified the right file name and location, then try again.")
        quit()
    if args.blastFile != None:
        if not os.path.isfile(args.blastFile):
            print(f"I was unable to locate the BLAST file at '{args.blastFile}'")
            print("Make sure you specified the right file name and location, then try again.")
            quit()
    # Ensure that the output location is sensible
    if args.outputTSV != None:
        if os.path.isfile(args.outputTSV):
            print(f"The specified output TSV '{args.outputTSV}' already exists. This script will not overwrite existing files.")
            print("Make sure to move this file or specify a different file name, then try again.")
            quit()
    if os.path.isfile(args.outputFile):
        print(f"The specified output file '{args.outputFile}' already exists. This script will not overwrite existing files.")
        print("Make sure to move this file or specify a different file name, then try again.")
        quit()

def parse_corset_clusters(corsetClustersFile):
    clustersDict = {}
    with open(corsetClustersFile , 'r') as fileIn:
        for line in fileIn:
            # Parse this line or skip if irrelevant
            sl = line.rstrip("\r\n ").split("\t")
            if len(sl) != 2:
                continue
            seqID, clustID = sl
            
            clustersDict.setdefault(clustID, [])
            clustersDict[clustID].append(seqID)
    
    return clustersDict

def parse_best_blast(blastFile, evalueCutoff):
    blastDict = {}
    with open(blastFile , 'r') as fileIn:
        for line in fileIn:
            # Parse this line or skip if irrelevant
            sl = line.rstrip("\r\n ").split("\t")
            if len(sl) != 12:
                continue
            qid, tid, identityPct, alignLength, \
                mismatch, gapopen, qstart, qend, \
                tstart, tend, evalue, bitscore = sl
            evalue = float(evalue)
            bitscore = float(bitscore)
            
            # Skip hits that don't pass E-value cutoff
            if evalue > evalueCutoff:
                continue
            
            # Store hit if it's superior to previous ones
            blastDict.setdefault(qid, [tid, bitscore])
            if blastDict[qid][1] < bitscore:
                blastDict[qid] = [tid, bitscore]
    
    return blastDict

def parse_annot_blast(annotTableFile, evalueCutoff):
    TID_COLUMN = "Target_accessions"
    EVALUE_COLUMN = "Expect_value"
    BIT_COLUMN = "Bit_score"
    
    blastDict = {}
    with open(annotTableFile , 'r') as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Parse header row
            if sl[0].startswith("#"):
                tidIndex = sl.index(TID_COLUMN)
                evalueIndex = sl.index(EVALUE_COLUMN)
                bitIndex = sl.index(BIT_COLUMN)
                continue
            # Skip if irrelevant
            if len(sl) != 12:
                continue
            
            # Parse details and skip if no hit recorded
            qid, tids, evalues, bitscores = sl[0], sl[tidIndex], sl[evalueIndex], sl[bitIndex]
            if tids == ".":
                continue
            tids = tids.replace(" ", "").replace("]", "").split("[")
            evalues = evalues.replace(" ", "").replace("]", "").split("[")
            bitscores = bitscores.replace(" ", "").replace("]", "").split("[")
            assert len(tids) == len(evalues)
            
            # Sort details to get the best hit
            hits = [[tids[i], float(evalues[i]), float(bitscores[i])] for i in range(len(tids))]
            hits.sort(key = lambda x: -x[2])
            
            tid, evalue, bitscore = hits[0]
            
            # Skip hits that don't pass E-value cutoff
            if evalue > evalueCutoff:
                continue
            
            # Store hit
            blastDict.setdefault(qid, [tid, bitscore])
            if blastDict[qid][1] < bitscore:
                blastDict[qid] = [tid, bitscore]
    
    return blastDict

def find_corset_representatives(clustersDict, seqLens, blastDict):
    '''
    This function is manually tuned to handle EvidentialGene transcripts
    and be friendly with utrorf sequences.
    '''
    # representatives = set()
    representatives = {}
    for clustID, seqIDs in clustersDict.items():
        # Handle single member clusters
        if len(seqIDs) == 1:
            # representatives.add(seqIDs[0])
            representatives[seqIDs[0]] = clustID
        # Handle multi-member clusters
        else:
            # Get evidence criteria sorted by priority
            repEvidence = [
                    [
                        seqid,
                        blastDict[seqid][1] if seqid in blastDict 
                            else blastDict[seqid + "utrorf"][1] if seqid + "utrorf" in blastDict
                            else -math.inf, # bitscore!
                        seqLens[seqid][0] if seqid in seqLens else seqLens[seqid + "utrorf"][0], # gives seqlen
                        seqLens[seqid][1] if seqid in seqLens else seqLens[seqid + "utrorf"][1] # gives numAmbiguous
                    ] 
                for seqid in seqIDs
            ]
            repEvidence.sort(key = lambda x: (-x[1], -x[2], x[3]))
            
            # Hold onto the best based on evidence sorting
            # representatives.add(repEvidence[0][0])
            representatives[repEvidence[0][0]] = clustID
    
    return representatives

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will read in a Corset results-clusters.txt file and
    the transcriptome FASTA it was run on and will output a FASTA file containing
    the best representative sequence from each cluster. Representative sequences
    will be ranked according to 1) their BLAST score (optionally) and/or their
    2) sequence length.
    """
    
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-c", dest="corsetClustersFile",
                   required=True,
                   help="Specify the location of the corset results-clusters.txt file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Specify the FASTA transcriptome file location")
    p.add_argument("-o", dest="outputFile",
                   required=True,
                   help="Specify the name of the output FASTA file")
    # Optional
    p.add_argument("-b", dest="blastFile",
                   required=False,
                   help="""Optionally, specify a BLAST file for help with
                   picking the representative sequence""")
    p.add_argument("--isAnnotTable", dest="isAnnotTable", action="store_true",
                   required=False,
                   help="""Specify this tag if the BLAST file being read in is actually
                   an annotation table; by default, it's expected to be outfmt6.""",
                   default=False)
    p.add_argument("--evalue", dest="evalue", type=float,
                   required=False,
                   help="""Optionally, if a BLAST file is being used for representative
                   picking, specify the minimum E-value required for a BLAST result
                   to be factored into the decision (default==1e-5)""",
                   default=1e-5)
    p.add_argument("--isNucleotide", dest="isNucleotide", action="store_true",
                   required=False,
                   help="""Specify this tag if the transcriptome file being read in contains
                   nucleotide sequences; by default, it is assumed to be protein.""",
                   default=False)
    p.add_argument("--writeClusterIDs", dest="writeClusterIDs", action="store_true",
                   required=False,
                   help="""Specify this tag if you'd like to write output sequences not with
                   their original transcript ID but with the Corset cluster ID""",
                   default=False)
    p.add_argument("--outputTSV", dest="outputTSV",
                   required=False,
                   help="""Optionally, specify a file name to write a TSV file pairing
                   cluster IDs to their representative sequence ID""")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse Corset clusters file
    clustersDict = parse_corset_clusters(args.corsetClustersFile)
    
    # Load in FASTA file as dict to get lengths and number of ambiguous characters
    records = SeqIO.parse(open(args.fastaFile, "r"), "fasta")
    if args.isNucleotide:
        seqLens = {r.id: [len(r), str(r.seq).upper().count("N")] for r in records}
    else:
        seqLens = {r.id: [len(r), str(r.seq).upper().count("X")] for r in records}
    
    # Load in BLAST results if relevant
    if args.blastFile != None:
        if args.isAnnotTable == True:
            blastDict = parse_annot_blast(args.blastFile, args.evalue) # handle annotation table
        else:
            blastDict = parse_best_blast(args.blastFile, args.evalue) # handle outfmt6
    else:
        blastDict = {}
    
    # Figure out our representative sequences
    representatives = find_corset_representatives(clustersDict, seqLens, blastDict)
    
    # Write representatives to file
    found = {}
    with open(args.outputFile, "w") as fileOut:
        records = SeqIO.parse(open(args.fastaFile, "r"), "fasta") # load in again for the sequences
        for r in records:
            if r.id in representatives:
                if args.writeClusterIDs is True:
                    fileOut.write(">{0} original_id={1}\n{2}\n".format(representatives[r.id], r.id, str(r.seq)))
                else:
                    fileOut.write(r.format("fasta"))
                
                # Hold onto pairing for validation and potential TSV output
                found[representatives[r.id]] = r.id
        
        # Try to rescue sequences not found if they're utrorf only sequences
        repSet = set(representatives.keys())
        diff = repSet.difference(set(found.values()))
        if len(diff) > 0:
            records = SeqIO.parse(open(args.fastaFile, "r"), "fasta") # load in again for the sequences...
            for r in records:
                if r.id.replace("utrorf", "") in diff:
                    if args.writeClusterIDs is True:
                        fileOut.write(">{0} original_id={1}\n{2}\n".format(representatives[r.id.replace("utrorf", "")], r.id, str(r.seq)))
                    else:
                        fileOut.write(r.format("fasta"))
                    
                    # Hold onto pairing for validation and potential TSV output
                    found[representatives[r.id.replace("utrorf", "")]] = r.id
    
    # Validate that it worked
    if len(found) == len(representatives):
        print(f"Found all {len(representatives)} expected sequences")
    else:
        print(f"Warning: Only found {len(found)} sequences, instead of the expected {len(representatives)}")
    
    # Write TSV if applicable
    if args.outputTSV != None:
        with open(args.outputTSV, "w") as fileOut:
            for clustID, repID in found.items():
                fileOut.write(f"{clustID}\t{repID}\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
