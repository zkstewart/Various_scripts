#! python3
# orthofinder_group_best_orthologs
# Script to read in an (ideally) abbreviated Orthogroups.csv file
# and identifies the best-matching ortholog for each gene in the query
# species against a single target species

# Load packages
import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_BlastIO, ZS_SeqIO

# Define functions for later use
## Validate arguments
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.orthogroupsFileName):
        print('I am unable to locate the Orthogroups file (' + args.orthogroupsFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.queryFasta):
        print('I am unable to locate the query FASTA file (' + args.queryFasta + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.targetFasta):
        print('I am unable to locate the target FASTA file (' + args.targetFasta + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def _tmp_file_name_gen(prefix, suffix):
    '''
    Hidden function for use by this program when creating temporary files.
    Params:
        prefix -- a string for a file prefix e.g., "tmp"
        suffix -- a string for a file suffix e.g., "fasta". Note that we don't
                  use a "." in this, since it's inserted between prefix and suffix
                  automatically.
    Returns:
        tmpName -- a string for a file name which does not exist in the current dir.
    '''
    ongoingCount = 1
    while True:
        if not os.path.isfile("{0}.{1}".format(prefix, suffix)):
            return "{0}.{1}".format(prefix, suffix)
        elif os.path.isfile("{0}.{1}.{2}".format(prefix, ongoingCount, suffix)):
            ongoingCount += 1
        else:
            return "{0}.{1}.{2}".format(prefix, ongoingCount, suffix)

def parse_orthogroups_tsv(tsvFileName):
    '''
    Parameters:
        tsvFileName -- the file name / location of the Orthogroups TSV file
    Returns:
        orthogroupsDict -- a dictionary with structure like:
                           {
                               "orthogroupID1": {
                                   "species1": [species1gene1, species1gene2, ...],
                                   "species2": [...],
                                   ...
                               },
                               "orthogroupID2": {...},
                               ...
                           }
                           Notable, it also contains indices with structure like:
                           {
                               "species1": [orthogroupID_containing_species1, ...],
                               "species2": [...],
                               ...
                           }
                           These latter indices indicate which orthogroups contain representatives
                           from the species.
                           
                           To discover what those indices are, use the "_headerKey" value i.e.,
                           orthogroupsDict["_headerKey"]
    '''
    orthogroupsDict = {}
    header = None
    
    with open(tsvFileName, 'r') as fileIn:
        for line in fileIn:
            sl = line.rstrip('\r\n"').split("\t")
            if header == None:
                header = sl[1:] # remove the "Orthogroup" value in the first column
                orthogroupsDict["_headerKey"] = header
            else:
                for i in range(1, len(sl)):
                    if sl[i] != "" and sl[i] != ".":
                        orthogroupsDict.setdefault(sl[0], {})
                        orthogroupsDict[sl[0]][header[i-1]] = sl[i].split(", ")
                        
                        orthogroupsDict.setdefault(header[i-1], [])
                        orthogroupsDict[header[i-1]].append(sl[0])
    return orthogroupsDict

def generate_FASTA_obj_from_orthogroup_sequences(orthogroupsDict, orthogroupIDs, speciesID, FASTA_obj):
    '''
    Generates a new FASTA object containing only the sequences found within
    the provided orthogroup IDs list.
    
    Parameters:
        orthogroupsDict -- a dictionary as created by parse_orthogroups_tsv().
        orthogroupIDs -- a list containing orthogroup IDs to retrieve from FASTA_obj.
        speciesID -- a string indicating the species column in the orthogroupsDict.
        FASTA_obj -- a ZS_SeqIO.FASTA object containing sequences to be retrieved.
    '''
    generated_FASTA_obj = ZS_SeqIO.FASTA(None)
    
    for orthogroup in orthogroupIDs:
        assert orthogroup in orthogroupsDict, \
            f"{orthogroup} not found in orthogroupsDict; generating FASTA object failed"
        assert speciesID in orthogroupsDict[orthogroup], \
            f"{speciesID} not found in {orthogroup}; generating FASTA object failed"
        
        sequenceIDs = orthogroupsDict[orthogroup][speciesID]
        for sequence in sequenceIDs:
            generated_FASTA_obj.insert(len(generated_FASTA_obj), FASTA_obj[sequence])
    
    return generated_FASTA_obj

def main():
    ### USER INPUT
    usage = """%(prog)s will parse the Orthogroups.tsv file output by Orthofinder alongside
    the FASTA file for the querySpecies and FASTA file for the targetSpecies to identify the
    best matching query:target sequence pairs within each orthogroup.
    
    Note that all groups which lack sequences from BOTH the query and target species will
    be ignored. 
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="orthogroupsFileName", required=True,
                help="Specify the orthogroup file")
    p.add_argument("-o", "-output", dest="outputFileName", required=True,
                help="Output file name")
    p.add_argument("-qs", dest="querySpecies", required=True,
                help="Column name of the query species")
    p.add_argument("-qi", dest="queryFasta", required=True,
                help="Input FASTA file for the query species")
    p.add_argument("-ts", dest="targetSpecies", required=True,
                help="Column name of the target species")
    p.add_argument("-ti", dest="targetFasta", required=True,
                help="Input FASTA file for the target species")
    # Opts
    p.add_argument("--blastAlg", dest="blastAlgorithm", required=False, choices=["blastp", "blastn"],
                help="Optionally, specify if both files are nucleotides (blastn) or proteins (blastp); default==blastp",
                default = "blastp")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse Orthofinder file
    orthogroupsDict = parse_orthogroups_tsv(args.orthogroupsFileName)
    
    # Validate that query and target species exist
    assert args.querySpecies in orthogroupsDict, \
        f"'{args.querySpecies}' could not be located within the OrthoGroups file; program ending now"
    assert args.targetSpecies in orthogroupsDict, \
        f"'{args.targetSpecies}' could not be located within the OrthoGroups file; program ending now"
    
    # Locate orthogroups that contain both query and target species
    orthogroupIDs = [
        key 
            for key, value in orthogroupsDict.items()
                if key not in orthogroupsDict["_headerKey"]
                and key != "_headerKey"
                and args.querySpecies in value
                and args.targetSpecies in value
    ]
    
    # Load in FASTA files
    queryFASTA_obj = ZS_SeqIO.FASTA(args.queryFasta)
    targetFASTA_obj = ZS_SeqIO.FASTA(args.targetFasta)
    
    # Generate FASTA objects containing the relevant sequences
    tmp_queryFASTA_obj = generate_FASTA_obj_from_orthogroup_sequences(orthogroupsDict, orthogroupIDs, args.querySpecies, queryFASTA_obj)
    tmp_targetFASTA_obj = generate_FASTA_obj_from_orthogroup_sequences(orthogroupsDict, orthogroupIDs, args.targetSpecies, targetFASTA_obj)
    
    # BLAST sequences
    blaster = ZS_BlastIO.BLAST(tmp_queryFASTA_obj, tmp_targetFASTA_obj, args.blastAlgorithm)
    blastDict, tmpResultName = blaster.get_blast_results()
    
    # Write results to file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("query_id\tbest_target_id\n")
        
        for query_FastASeq_obj in tmp_queryFASTA_obj:
            queryID = query_FastASeq_obj.id
            if queryID not in blastDict:
                fileOut.write(f"{queryID}\t.\n")
            else:
                fileOut.write(f"{queryID}\t{blastDict[queryID][0][0]}\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
