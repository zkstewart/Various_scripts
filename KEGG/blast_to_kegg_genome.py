#! python3
# blast_to_kegg_genome.py
# This script will facilitate a BLAST of the input
# gene sequences against a KEGG-annotated genome to
# map the KEGG terms from the official genome annotation
# to the one in question.

import sys, argparse, os, requests, re
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_BlastIO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.inputFastaFile):
        print('I am unable to locate the input protein FASTA file (' + args.inputFastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.targetFastaFile):
        print('I am unable to locate the target KEGG protein FASTA file (' + args.targetFastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.blastFile != None:
        if not os.path.isfile(args.blastFile):
            print('I am unable to locate the pre-computed BLAST outfmt6 file (' + args.blastFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate KEGG ID
    isValid, keggID = kegg_id_is_valid(args.keggID)
    if not isValid:
        print('The provided KEGG ID isn\'t recognised by the KEGG API (' + args.keggID + ')')
        print('Make sure you\'ve typed it in correctly and try again.')
        quit()
    else:
        args.keggID = keggID # this will be case-corrected
    # Validate numeric inputs
    if args.threads < 1:
        print("Threads argument must be a positive integer greater than 0")
    if not 0 <= args.identity <= 100:
        print("Identity argument must be a positive float in the range 0->100 (inclusive)")
    # Handle file output
    outputSuffixes = ["_keggmap.tsv", "_keggdescriptions.tsv"]
    for suffix in outputSuffixes:
        if os.path.isfile(args.outputFilePrefix + suffix):
            print('One of the output files already exists (' + args.outputFilePrefix + suffix + ')')
            print('This program will not allowing overwriting; try again with a different prefix, or move/rename the existing file.')
            quit()

def kegg_id_is_valid(keggID):
    '''
    Parameters:
        keggID -- a string with appropriate case indicating the KEGG organism ID
    Returns:
        isValid -- a Boolean indicating whether the ID was recognised or not
        keggID -- a case-corrected KEGG ID (just in case the input wasn't in the right case)
    '''
    API_URL = "https://rest.kegg.jp/list/organism"
    
    response = requests.get(API_URL)
    for line in response.text.split("\n"):
        sl = line.split("\t")
        if sl[0] == keggID.upper():
            return True, keggID.upper()
        elif sl[1] == keggID.lower():
            return True, keggID.lower()
    
    return False, None

def download_ncbi_protein_to_kegg_gene_map(keggID):
    '''
    Parameters:
        keggID -- a string with appropriate case indicating the KEGG organism ID
    Returns:
        proteinMapDict -- a dictionary with structure like:
                            {
                                ncbi_protein_ID: kegg_gene_ID,
                                ...
                            }
    '''
    api_url = "https://rest.kegg.jp/conv/ncbi-proteinid/{0}".format(keggID)
    proteinMapDict = {}
    
    response = requests.get(api_url)
    for line in response.text.split("\n"):
        if line == "": # skip empty lines
            continue
        keggID, ncbiID = line.split("\t")
        keggID = keggID.split(":")[1]
        ncbiID = ncbiID.split(":")[1]
        proteinMapDict[ncbiID] = keggID
    
    return proteinMapDict
    
def download_kegg_gene_to_kegg_pathway_map(keggID):
    '''
    Parameters:
        keggID -- a string with appropriate case indicating the KEGG organism ID
    Returns:
        keggPathwayMap -- a dictionary with structure like:
                            {
                                kegg_gene_ID: kegg_pathway_ID,
                                ...
                            }
    '''
    api_url = "https://rest.kegg.jp/link/pathway/{0}".format(keggID)
    keggPathwayMap = {}
    
    response = requests.get(api_url)
    for line in response.text.split("\n"):
        if line == "": # skip empty lines
            continue
        keggGeneID, keggPathwayID = line.split("\t")
        keggGeneID = keggGeneID.split(":")[1]
        keggPathwayID = keggPathwayID.split(":")[1]
        keggPathwayMap[keggGeneID] = keggPathwayID
    
    return keggPathwayMap

def download_species_pathway_to_KO_information(pathwayID):
    '''
    Parameters:
        pathwayID -- a string with appropriate case indicating the KEGG pathway ID
                     for the species in question (e.g., it might be "ming01200" rather
                     than "ko01200")
    Returns:
        koPathwayID -- a string with appropriate case indicating the KO ID for the
                       pathway provided as the input parameter.
        koName -- a string indicating the name for this KO pathway
        KoDescription -- a string containing the paragraph of information detailing the
                         KO pathway.
    '''
    api_url1 = "https://rest.kegg.jp/get/{0}".format(pathwayID)
    koPathwayRegex = re.compile(r"KO_PATHWAY(.+?)\n")
    nameRegex = re.compile(r"NAME(.+?)\n")
    descriptionRegex = re.compile(r"DESCRIPTION(.+?)\n")
    
    # Get the KO pathway ID from the first response
    response1 = requests.get(api_url1)
    koPathwayID = koPathwayRegex.findall(response1.text)[0].strip(" ")
    
    # Get the KO details from a second response
    api_url2 = "https://rest.kegg.jp/get/{0}".format(koPathwayID)
    response2 = requests.get(api_url2)
    koName = nameRegex.findall(response2.text)[0].strip(" ")
    try:
        koDescription = descriptionRegex.findall(response2.text)[0].strip(" ")
    except: # handles values where no description exists
        koDescription = "."
    
    return koPathwayID, koName, koDescription

if __name__ == "__main__":
    usage = """%(prog)s receives an input FASTA containing protein sequences, and
    uses BLAST to map annotations to these sequences based on the best hit against
    an official KEGG annotated protein FASTA. This is useful when using a newer genome
    that hasn't had a KEGG annotation performed for it already!
    
    The output will be two TSV files. The first has goseq formatting i.e., "geneID\tkeggID"
    or "geneID\t0" when a hit was not found. The second relates KEGG orthology IDs to their
    name and description where applicable i.e., "tkeggID\tkeggName\tkeggDescription". When
    a description does not exist, a simple "." will be shown instead.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputFastaFile", required=True,
                help="Specify the location of the input gene files FASTA")
    p.add_argument("-t", dest="targetFastaFile", required=True,
                help="""Specify the location of the target gene files FASTA that
                KEGG annotations are being retrieved from.""")
    p.add_argument("-k", dest="keggID", required=True,
                help="Specify the ID of the KEGG organism mappings are being retrieved from (e.g., T07436)")
    p.add_argument("-o", dest="outputFilePrefix", required=True,
                help="Output file prefix for KEGG mapping TSV and KEGG description TSV files")
    # Opts
    p.add_argument("--threads", dest="threads", required=False, type=int,
                help="Optionally, specify how many threads to run BLAST with (default==\"1\")",
                default="1")
    p.add_argument("--identity", dest="identity", required=False, type=float,
                help="Optionally, specify the minimum identity percent allowed for a mapping to occur (default==\"90.0\")",
                default="90.0")
    p.add_argument("--blastFile", dest="blastFile", required=False,
                help="Optionally, specify an already existing outfmt6 file containing BLAST results")
    args = p.parse_args()
    validate_args(args)
    
    # Set up our BLAST handler
    blaster = ZS_BlastIO.BLAST(args.inputFastaFile, args.targetFastaFile, "blastp")
    blaster.set_threads(args.threads)
    
    # Run BLAST
    if args.blastFile == None:
        blastDict, _ = blaster.get_blast_results() # throw away the blastFileName return since it will be None
    else:
        blastDict = blaster.parse_blast_hit_coords(args.blastFile)
    
    # Query KEGG API for mapping information
    proteinMapDict = download_ncbi_protein_to_kegg_gene_map(args.keggID)
    keggPathwayMap = download_kegg_gene_to_kegg_pathway_map(args.keggID)
    
    # Map input FASTA IDs to KEGG pathways
    inputFasta = ZS_SeqIO.FASTA(args.inputFastaFile)
    inputMapDict = {}
    descriptionDict = {}
    foundDict = {}
    for FastASeq_obj in inputFasta:
        inputID = FastASeq_obj.id
        
        # If there's no BLAST hit, no mapping can commence
        if inputID not in blastDict:
            inputMapDict[inputID] = "0"
            continue
        
        # Retrieve the best hit, and skip if it's not good enough
        bestID, identity, qstart, qend, tstart, tend, evalue = blastDict[inputID][0]
        if identity < args.identity: # if it's not a good enough hit, skip it
            inputMapDict[inputID] = "0"
            continue
        
        # If there's no protein ID mapping, we can't proceed
        bestID = bestID.rsplit(".", maxsplit=1)[0] # correct isoform IDs
        if bestID not in proteinMapDict:
            inputMapDict[inputID] = "0"
            continue
        
        # If we can't map the KEGG gene ID to its pathway, we can't proceed
        bestKeggID = proteinMapDict[bestID]
        if bestKeggID not in keggPathwayMap:
            inputMapDict[inputID] = "0"
            continue
        
        # If we can do the mapping, grab the KEGG pathway and associated information
        else:
            keggPathwayID = keggPathwayMap[bestKeggID]
            
            # Skip API calls if this pathway has already been queried
            if keggPathwayID in foundDict:
                koPathwayID, koName, koDescription = foundDict[keggPathwayID]
            # Otherwise, make the API call
            else:
                koPathwayID, koName, koDescription = download_species_pathway_to_KO_information(keggPathwayID)
                foundDict[keggPathwayID] = [koPathwayID, koName, koDescription]
            
            # Associate results to our relevant dictionaries
            inputMapDict[inputID] = koPathwayID
            descriptionDict[koPathwayID] = [koName, koDescription]
    
    # Write output KEGG pathway mapping
    with open(args.outputFilePrefix + "_keggmap.tsv", "w") as fileOut:
        for key, value in inputMapDict.items():
            fileOut.write("{0}\t{1}\n".format(key, value))
    
    # Wrie output KEGG description and name mapping
    with open(args.outputFilePrefix + "_keggdescriptions.tsv", "w") as fileOut:
        for key, value in descriptionDict.items():
            name, description = value
            fileOut.write("{0}\t{1}\t{2}\n".format(key, name, description))
        
    
    print("Program completed successfully!")
