#! python3
# haplotype_concat.py
# Script for working with a MSA FASTA file resulting from the
# gene_haplotyping_pipeline.py script or similar, to concatenate
# each haplotype into a single sequence for further analysis
# e.g., phylogeny.

import sys, argparse, os
import numpy as np
from collections import Counter
from Levenshtein import distance

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO

def validate_args(args):
    # Validate input data location
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the input FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print('I am unable to locate the metadata file (' + args.metadataFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a new file name or move/rename the existing file.')
        quit()

def validate_haplo_seqs(FASTA_obj):
    # Ensure that FASTA sequence IDs meet the expected format
    haploDict = {}
    matching, missing = set(), set()
    for FastASeq_obj in FASTA_obj:
        # Get details out of the sequence ID
        try:
            seqID, haplotype, haploNum = FastASeq_obj.id.rsplit("_", maxsplit=2)
        except:
            raise ValueError((f"Sequence '{FastASeq_obj.id}' does not have the expected format of " + 
                             "'geneID_haplotype_#'"))
        
        if not haplotype == "haplotype":
            raise ValueError(f"Sequence '{FastASeq_obj.id}' does not end with _haplotype* as expected")
        
        # Ensure that the haplotype ID is a number
        if not haploNum.isdigit():
            raise ValueError(f"Sequence '{FastASeq_obj.id}' does not end with a digit as expected")
        haploNum = int(haploNum)
        
        # Ensure that the haplotype ID is unique, then store
        haploDict.setdefault(seqID, set())
        if haploNum in haploDict[seqID]:
            raise ValueError(f"Sequence '{FastASeq_obj.id}' has a duplicate haplotype ID")
        haploDict[seqID].add(haploNum)
    
    # Ensure that all haplotypes are present for each gene
    haploSets = list(haploDict.values())
    for s in haploSets[1:]:
        if s != haploSets[0]:
            raise ValueError(("Not all genes have the same haplotypes e.g., the first sequence " +
                             f"has haplotypes {haploSets[0]} while another has {s}"))
    
    return s

def validate_metadata(FASTA_obj, metadataDict):
    # Update metadata to match FASTA sequence IDs
    updatedMetadataDict = {}
    found = set()
    missing = set()
    for FastASeq_obj in FASTA_obj:
        # Get details out of the sequence ID
        seqID, haplotype, haploNum = FastASeq_obj.id.rsplit("_", maxsplit=2)
        
        # Check if the sequence ID is in the metadata
        if seqID in metadataDict:
            updatedMetadataDict[seqID] = metadataDict[seqID]
            found.add(seqID)
        else:
            # See if we can trim the suffix to find a match
            added = False
            origSeqID = seqID
            while seqID.count("_") > 1:
                seqID = "_".join(seqID.split("_")[:-1])
                if seqID in metadataDict:
                    updatedMetadataDict[origSeqID] = metadataDict[seqID]
                    added = True
                    found.add(seqID)
                    break
            if not added:
                missing.add(origSeqID)
    
    # Make sure each sequence ID in the FASTA file is in the metadata file
    if len(missing) > 0:
        raise ValueError(f"The following FASTA sequence IDs were not found in the metadata file: {missing}")
    
    # Alert user to any metadata values absent from the FASTA file
    if len(updatedMetadataDict) < len(metadataDict):
        missing = set(metadataDict.keys()).difference(found)
        print(f"# WARNING: The following metadata values were not found in the FASTA file: {missing}")
        print("# This might be normal if your metadata file contains more entries than the FASTA file")
    
    return updatedMetadataDict

def parse_metadata(metadataFile):
    '''
    Function to parse a metadata file and return a dictionary
    of the data. Skips any lines that start with a hash (#).
    
    Parameters:
        metadataFile -- a string containing the location of the metadata file.
    Returns:
        metadataDict -- a dictionary with the gene ID as the key and the
                        species/variety name as the value
    '''
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            # Skip comments
            if line.startswith("#"):
                continue
            
            # Parse line
            delimiter = "\t" if "\t" in line else ","
            sl = line.strip().split(delimiter)
            if len(sl) != 2:
                raise ValueError(f"Metadata file '{metadataFile}' does not have 2 columns as expected")
            
            # Format speciesID for compliance with good FASTA formatting
            seqID, speciesID = sl
            metadataDict[seqID] = ZS_SeqIO.Conversion.format_str_as_fasta_id(speciesID)
    return metadataDict

def generate_haplotype_consensi(FASTA_obj, numHaplotypes=2, isNucleotide=False):
    '''
    Generates consensus sequences for each haplotype in the FASTA object. The
    method of doing so is heuristic in nature. The first haplotype's consensus will
    be the true consensus of the alignment. The second haplotype will eliminate any
    positions that are identical to the first haplotype's consensus, and substitute
    the next most common variant residue. This process will continue for each haplotype.
    
    Parameters:
        FASTA_obj -- a ZS_SeqIO.FASTA object.
        numHaplotypes -- an integer specifying the number of haplotypes in the FASTA object.
        isNucleotide -- a boolean specifying if the sequences are nucleotide sequences.
    Returns:
        consensi -- a list of consensus sequences (as ZS_SeqIO.FastASeq objects) for
                    each haplotype.
    '''
    stepSize = 3 if isNucleotide else 1
    seqLen = len(FASTA_obj.seqs[0].gap_seq) # all sequences should be the same length
    
    # Generate consensus sequences
    consensi = [""] * numHaplotypes
    for i in range(0, seqLen, stepSize):
        # Get the most common residue at this position
        positionList = []
        for x in range(len(FASTA_obj.seqs)):
            positionList.append(FASTA_obj.seqs[x].gap_seq[i:i+stepSize])
        
        positionCount = Counter(positionList)
        mostCommon = positionCount.most_common()
        
        # Handle non-variant columns
        if len(mostCommon) == 1:
            for x in range(numHaplotypes):
                consensi[x] += mostCommon[0][0]
        # Handle variant columns
        else:
            for x in range(numHaplotypes):
                consensi[x] += mostCommon[x][0]
    
    return [ ZS_SeqIO.FastASeq(f"haplotype_{x+1}", gapSeq=consensi[x]) for x in range(numHaplotypes) ]

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will receive a FASTA file that contains multiple haplotype variants
    for a gene from the gene_haplotyping_pipeline.py script or similar, and concatenate
    them into a MSA suitable for subsequent analyses e.g., phylogenetics. The metadata file
    should contain the gene ID and the species/variety name for each haplotype to provide
    a better identifier for each sequence; if you don't want to rename sequences, you should
    still provide a metadata file with the sequence ID (sans haplotype suffix) with this
    value duplicated into the second column. Note: this script is tolerant to some metadata
    ID mismatches, but it will alert you to any discrepancies.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="fastaFile",
                   required=True,
                   help="""Specify the location of the input FASTA file
                   containing haplotype variants for a gene""")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="""Specify the location of the metadata
                   TSV or CSV file""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the file name for the output MSA FASTA")
    # Opts
    p.add_argument("--isNucleotide", dest="isNucleotide",
                   required=False,
                   action="store_true",
                   help="Provide this flag if the input file contains nucleotide sequences",
                   default=False)
    args = p.parse_args()
    validate_args(args)
    
    # Load FASTA file
    FASTA_obj = ZS_SeqIO.FASTA(args.fastaFile, isAligned=True)
    if len(FASTA_obj) < 2:
        raise ValueError(f"FASTA file '{args.fastaFile}' must contain at least two sequences")
    
    # Load metadata file
    metadataDict = parse_metadata(args.metadataFile)
    
    # Validate haplotype sequences
    haplotypeNums = validate_haplo_seqs(FASTA_obj)
    for x in range(len(haplotypeNums)):
        if x+1 not in haplotypeNums:
            raise ValueError(f"Haplotype numbers must be sequential starting from 1; yours are {haplotypeNums}")
    
    # Validate and update metadata
    updatedMetadataDict = validate_metadata(FASTA_obj, metadataDict)
    
    # Heuristically generate consensus sequences for each haplotype
    consensi = generate_haplotype_consensi(FASTA_obj, len(haplotypeNums), args.isNucleotide)
    
    # Generate all combinations of haplotype indices
    combinations = np.array(
        np.meshgrid(*[ list(range(0, len(haplotypeNums))) for x in range(len(haplotypeNums)) ])
    ).T.reshape(-1, len(haplotypeNums))
    combinations = np.array([ x for x in combinations if len(set(x)) == len(haplotypeNums) ])
    
    # Cluster sequences into haplotype groups
    consensusGroups = [ [] for x in range(len(haplotypeNums)) ]
    for seqPrefix in updatedMetadataDict.keys():
        # Grab each haplotype for this species/variety
        queries = []
        for x in range(len(haplotypeNums)):
            seqID = f"{seqPrefix}_haplotype_{x+1}"
            for FastASeq_obj in FASTA_obj:
                if FastASeq_obj.id == seqID:
                    queries.append(FastASeq_obj)
                    break
        
        # Sanity check
        "Queries should be an empty list (metadata has extra values) or length equal to haplotypeNums"
        if queries == []:
            continue
        elif len(queries) != len(haplotypeNums):
            raise ValueError(f"Unexpected error occurred while clustering sequences for '{seqPrefix}'")
        
        # Distance each query to each consensus sequences
        results = []
        for query_FastASeq in queries:
            results.append([])
            for targetFastASeq_obj in consensi:
                dist = distance(query_FastASeq.gap_seq, targetFastASeq_obj.gap_seq)
                results[-1].append(dist)
        
        # Brute force to find optimal solution
        bestCombination = [np.inf, None]
        for combination in combinations:
            # Sum the distance value for this combination
            score = 0
            for queryIndex, targetIndex in enumerate(combination):
                score += results[queryIndex][targetIndex]
            if score < bestCombination[0]:
                bestCombination = [score, combination]
        
        # Assign queries according to score optimisation in bestCombination
        for queryIndex, targetIndex in enumerate(bestCombination[1]):
            consensusGroups[targetIndex].append(queries[queryIndex])
    
    # Turn each consensus group into a FASTA object
    consensusFASTAs = []
    for group in consensusGroups:
        groupFASTA = ZS_SeqIO.FASTA(None, isAligned=True)
        for FastASeq_obj in group:
            groupFASTA.add(FastASeq_obj)
        consensusFASTAs.append(groupFASTA)
    
    # Concatenate FASTA objects
    for x in range(1, len(consensusFASTAs)):
        consensusFASTAs[0].concat(consensusFASTAs[x]) # it all gets concatenated into the first FASTA object
    concatenatedFASTA = consensusFASTAs[0]
    
    # Update IDs for each sequence
    for FastASeq_obj in concatenatedFASTA:
        seqID, haplotype, haploNum = FastASeq_obj.id.rsplit("_", maxsplit=2)
        FastASeq_obj.id = f"{updatedMetadataDict[seqID]}"
    
    # Write output
    concatenatedFASTA.write(args.outputFileName, asAligned=True)
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
