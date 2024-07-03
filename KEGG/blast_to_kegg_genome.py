#! python3
# blast_to_kegg_genome.py
# This script will facilitate a BLAST of the input
# gene sequences against a KEGG-annotated genome to
# map the KEGG terms from the official genome annotation
# to the one in question.

import sys, argparse, os, requests, re, pickle
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_BlastIO

API_PICKLE_FILE = ".kegg_api_data.pkl" # global this for convenience
BRITE_REGEX = re.compile(r"BRITE(.+?)(POSITION|MOTIF)", re.DOTALL) # sometimes POSITION isn't present
ORTHO_REGEX = re.compile(r"ORTHOLOGY\s+(K\d{5})\s+(.+?)\n")
KO_REGEX = re.compile(r"(\d{5})\s(.+?)\n")

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
        print("--threads must be a positive integer greater than 0")
    if not 0 <= args.identity <= 100:
        print("--identity must be a positive float in the range 0->100 (inclusive)")
    if args.evalue < 0:
        print("--evalue must be an integer >= 0")
    # Handle file output
    outputSuffixes = ["_keggmap.tsv", "_keggnames.tsv"]
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

def download_kegg_gene_to_kegg_ko_map(organismID):
    '''
    Parameters:
        organismID -- a string with appropriate case indicating the KEGG organism ID
    Returns:
        keggKoMap -- a dictionary with structure like:
                            {
                                kegg_gene_ID: kegg_KO_id,
                                ...
                            }
    '''
    api_url = "https://rest.kegg.jp/link/ko/{0}".format(organismID)
    keggKoMap = {}
    
    response = requests.get(api_url)
    for line in response.text.split("\n"):
        if line == "": # skip empty lines
            continue
        keggGeneID, keggKO = line.split("\t")
        keggGeneID = keggGeneID.split(":")[1]
        keggKO = keggKO.split(":")[1]
        keggKoMap[keggGeneID] = keggKO
    
    return keggKoMap

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

def parse_ko_get_response(responseText):
    '''
    Parameters:
        responseText -- a string containing the response from a KEGG API query
    Returns:
        responseDict -- a dictionary with structure attempting to render the API
                        response into something JSON-like. It should look like:
                        {
                            "ENTRY": "K00001",
                            "SYMBOL": "adh",
                            "NAME": "alcohol dehydrogenase [EC: ...],
                            "PATHWAY": ["map00010", "Glycolysis / Gluconeogenesis"],
                            "BRITE": [
                                ["09130", "Environmental Information Processing"],
                                ["09131", "Membrane transport"],
                                [ ... ],
                                ...
                            ]
                        }
    '''
    splitter = re.compile(r"\s+")
    briteRegex = re.compile(r"\d{5}\s")
    responseDict = {}
    
    # Iterate through response text by line and parse out the data
    lastParentKey = None
    for line in responseText.rstrip("\n ").split("\n"):
        if line.startswith("///"):
            break
        
        key, value = splitter.split(line, maxsplit=1)
        value = " ".join(value.split()) # trim whitespace
        
        # Handle new parent keys
        if line[0] != " ":
            # Index parent keys by type
            if key == "ENTRY":
                responseDict[key] = value.split(" ")[0] # remove the KO at the end
            elif key == "SYMBOL":
                responseDict[key] = value # no changes needed
            elif key == "NAME":
                responseDict[key] = value.split(" [EC")[0] # remove the enzyme code
            elif key == "PATHWAY":
                responseDict[key] = [value.split(" ", maxsplit=1)]
            elif key == "MODULE":
                responseDict[key] = [value.split(" ", maxsplit=1)]
            elif key == "REACTION":
                try:
                    reactionKey, reactionValue = value.split(" ", maxsplit=1)
                except:
                    reactionKey, reactionValue = value, "."
                responseDict[key] = [[reactionKey, reactionValue]]
            elif key == "DISEASE":
                try:
                    diseaseKey, diseaseValue = value.split(" ", maxsplit=1)
                except:
                    diseaseKey, diseaseValue = value, "."
                responseDict[key] = [[diseaseKey, diseaseValue]]
            elif key == "BRITE":
                responseDict[key] = []
            elif key == "DBLINKS":
                responseDict[key] = {}
            elif key == "GENES":
                responseDict[key] = {}
            elif key == "REFERENCE":
                responseDict.setdefault(key, [])
                responseDict[key].append({})
                if value != "" and ":" in value:
                    refKey, refValue = value.split(":", maxsplit=1)
                    responseDict[key].append({ refKey : refValue })
            else:
                raise NotImplementedError(f"Unrecognized key '{key}' in KEGG API response for '{responseDict['ENTRY']}'")
            
            # Hold onto which parent we're indexing
            lastParentKey = key
        # Handle child keys
        else:
            # Handle PATHWAY entries
            if lastParentKey == "PATHWAY":
                pathwayKey, pathwayValue = value.split(" ", maxsplit=1)
                responseDict[lastParentKey].append([pathwayKey, pathwayValue])
            
            # Handle MODULE entries
            elif lastParentKey == "MODULE":
                moduleKey, moduleValue = value.split(" ", maxsplit=1)
                responseDict[lastParentKey].append([moduleKey, moduleValue])
            
            # Handle REACTION entries
            elif lastParentKey == "REACTION":
                try:
                    reactionKey, reactionValue = value.split(" ", maxsplit=1)
                    responseDict[lastParentKey].append([reactionKey, reactionValue])
                except:
                    "This can happen when a reaction doesn't have an associated name, just a code"
                    reactionKey = value
                    responseDict[lastParentKey].append([reactionKey, "."])
            
            # Handle DISEASE entries
            elif lastParentKey == "DISEASE":
                try:
                    diseaseKey, diseaseValue = value.split(" ", maxsplit=1)
                    responseDict[lastParentKey].append([diseaseKey, diseaseValue])
                except:
                    "This can happen when a disease doesn't have an associated name, just a code"
                    diseaseKey = value
                    responseDict[lastParentKey].append([diseaseKey, "."])
            
            # Handle BRITE entries
            elif lastParentKey == "BRITE":
                if not briteRegex.match(value): # skip over non-BRITE entries
                    continue
                value = value.split(" ", maxsplit=1)
                
                if value[1] == "Brite Hierarchies": # irrelevant value
                    continue
                
                responseDict[lastParentKey].append(value)
            
            # Handle DBLINKS entries
            elif lastParentKey == "DBLINKS":
                dbKey, dbValue = value.split(": ", maxsplit=1)
                responseDict[lastParentKey].setdefault(dbKey, set())
                responseDict[lastParentKey][dbKey].add(dbValue)
            
            # Handle GENES entries
            elif lastParentKey == "GENES":
                geneKey, geneValue = value.split(": ", maxsplit=1)
                
                responseDict[lastParentKey].setdefault(geneKey, [])
                responseDict[lastParentKey][geneKey].append(geneValue)
            
            # Handle REFERENCE entries
            elif lastParentKey == "REFERENCE":
                try:
                    refKey, refValue = value.split(" ", maxsplit=1)
                    responseDict[lastParentKey][-1][refKey] = refValue
                except: # handle tagger-onners like DOIs for the JOURNAL section
                    if type(responseDict[lastParentKey][-1][refKey]) == str:
                        responseDict[lastParentKey][-1][refKey] = [responseDict[lastParentKey][-1][refKey]]
                    responseDict[lastParentKey][-1][refKey].append(value)
            else:
                raise NotImplementedError(f"Unrecognized child key '{key}' in KEGG API response for '{responseDict['ENTRY']}'")
    return responseDict

def download_ko_information(koID):
    '''
    Parameters:
        koID -- a string indicating the KO to query.
    Returns:
        responseDict -- a dictionary with structure attempting to render the API
                        response into something JSON-like. Same output as given
                        by parse_ko_get_response().
    '''
    api_url = f"https://rest.kegg.jp/get/{koID}"
    
    # Get the response
    response = requests.get(api_url)
    
    # Handle empty response value
    "This probably means the KEGG database has a mistake between protein map and kegg map data"
    if response.text == "":
        return None
    
    # Parse response into dictionary
    responseDict = parse_ko_get_response(response.text)
    
    return responseDict

def save_pickle(queriedAccs):
    if queriedAccs != {}:
        with open(API_PICKLE_FILE, "wb") as pickleOut:
            pickle.dump(queriedAccs, pickleOut)

def main():
    usage = """%(prog)s receives an input FASTA containing protein sequences, and
    uses BLAST to map annotations to these sequences based on the best hit against
    an official KEGG annotated protein FASTA. You can provide pre-computed BLAST results
    throuhg --blastFile, but you still need to provide the -i and -t files.
    
    The output will be two TSV files. The first has goseq formatting i.e., "geneID\tkeggID"
    or "geneID\t0" when a hit was not found. The second relates KEGG orthology IDs to their
    name i.e., "tkeggID\tkeggName".
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputFastaFile",
                   required=True,
                   help="Specify the location of the input gene files FASTA")
    p.add_argument("-t", dest="targetFastaFile",
                   required=True,
                   help="""Specify the location of the target gene files FASTA that
                   KEGG annotations are being retrieved from.""")
    p.add_argument("-k", dest="keggID",
                   required=True,
                   help="Specify the ID of the KEGG organism mappings are being retrieved from (e.g., T07436)")
    p.add_argument("-o", dest="outputFilePrefix",
                   required=True,
                   help="Output file prefix for KEGG mapping TSV and KEGG description TSV files")
    # Opts
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="Optionally, specify how many threads to run BLAST with (default==1)",
                   default=1)
    p.add_argument("--identity", dest="identity",
                   required=False,
                   type=float,
                   help="""Optionally, specify the minimum identity percent necessary for a
                   mapping to occur (default==90.0)""",
                   default=90.0)
    p.add_argument("--evalue", dest="evalue",
                   required=False,
                   type=float,
                   help="""Optionally, specify the minimum evalue necessary for a
                   mapping to occur (default==1e-3)""",
                   default=1e-3)
    p.add_argument("--blastFile", dest="blastFile",
                   required=False,
                   help="""Optionally, specify an already existing outfmt6 file containing
                   BLAST results""")
    args = p.parse_args()
    validate_args(args)
    
    # Query KEGG API for mapping information
    proteinMapDict = download_ncbi_protein_to_kegg_gene_map(args.keggID)
    keggKoMap = download_kegg_gene_to_kegg_ko_map(args.keggID)
    
    # Obtain target sequence IDs
    targetIDs = set()
    with open(args.targetFastaFile, "r") as fileIn:
        records = SeqIO.parse(fileIn, "fasta")
        for record in records:
            targetIDs.add(record.id)
    
    # Check that the target IDs are in the proteinMapDict
    "If this isn't true, the target file is incorrect or maybe has isoform IDs"
    overlap = targetIDs.intersection(set(proteinMapDict.keys()))
    if len(overlap) == 0:
        print(f"# ERROR: No sequence IDs in '{args.targetFastaFile}' overlap with the KEGG " + 
              f"protein mapping dictionary for '{args.keggID}'")
        print("# Your target file may have isoform IDs, the KEGG ID may be incorrect, or " +
              "something else may be wrong.")
        print("# Please try to figure this out and then try again.")
        quit()
    else:
        print(f"# Note: {len(overlap)} out of {len(targetIDs)} sequence IDs in " + 
              f"'{args.targetFastaFile}' have protein mappings for '{args.keggID}'")
    
    # Set up our BLAST handler
    blaster = ZS_BlastIO.BLAST(args.inputFastaFile, args.targetFastaFile, "blastp")
    blaster.set_threads(args.threads)
    blaster.set_clean(False) # keep the BLAST file for now, might as well?
    
    # Run BLAST
    if args.blastFile == None:
        blastDict, _ = blaster.get_blast_results() # throw away the blastFileName return since it will be None
    else:
        blastResults = ZS_BlastIO.BLAST_Results(args.blastFile)
        blastResults.parse()
        blastDict = blastResults.results
    
    # Load in any API queries that may have been performed already
    if os.path.isfile(API_PICKLE_FILE):
        with open(API_PICKLE_FILE, "rb") as pickleIn:
            apiDict = pickle.load(pickleIn)
    else:
        apiDict = {}
    
    # Obtain query sequence IDs
    queryIDs = set()
    with open(args.inputFastaFile, "r") as fileIn:
        records = SeqIO.parse(fileIn, "fasta")
        for record in records:
            queryIDs.add(record.id)
    
    # Map query IDs to KEGG information
    keggmapDict = {}
    keggnamesDict = {}
    try:
        for inputID in queryIDs:            
            # If there's no BLAST hit, no mapping can commence
            if inputID not in blastDict:
                keggmapDict[inputID] = "0"
                continue
            
            # Retrieve the best hit, and skip if it's not good enough
            bestID, identity, qstart, qend, tstart, tend, evalue = blastDict[inputID][0]
            if identity < args.identity or evalue > args.evalue:
                keggmapDict[inputID] = "0"
                continue
            
            # Identify protein mapping
            if bestID not in proteinMapDict:
                keggmapDict[inputID] = "0" # can't proceed without a KEGG gene ID
                continue
            bestKeggID = proteinMapDict[bestID]
            
            # Identify KO mapping
            if bestKeggID not in keggKoMap:
                keggmapDict[inputID] = "0" # can't proceed without a KEGG KO
                continue
            bestKO = keggKoMap[bestKeggID]
            
            # Obtain data via API (if not already done so)
            if bestKO in apiDict:
                koData = apiDict[bestKO]
            else:
                koData = download_ko_information(bestKO)
                if koData == None: # this can happen if the KEGG database has a mistake
                    keggmapDict[inputID] = "0" # can't proceed without KO data
                    continue
                else:
                    apiDict[bestKO] = koData
            
            # Store keggmap and keggnames data
            keggmapDict[inputID] = set()
            
            for parentKey in ["PATHWAY", "BRITE"]:
                if parentKey in koData:
                    for mapping, name in koData[parentKey]:
                        keggmapDict[inputID].add(mapping)
                        keggnamesDict[mapping] = name
    
    # If program is ending unsuccessfuly, save any API queries now and provide debug info
    except Exception as e:
        save_pickle(apiDict)
        
        print("## DEBUG:")
        print("## Program broke on sequence ID={0}".format(inputID))
        print("## Best KO = {0}".format(bestKO))
        print(f"## Exception message is '{str(e)}'")
        
        print("Program ended after failing =(")
        quit()
    
    # Save any API queries that may have been performed after finishing successfully
    save_pickle(apiDict)
    
    # Write output KEGG pathway mapping
    with open(args.outputFilePrefix + "_keggmap.tsv", "w") as fileOut:
        for key, value in keggmapDict.items():
            fileOut.write("{0}\t{1}\n".format(key, "; ".join(value)))
    
    # Wrie output KEGG description and name mapping
    with open(args.outputFilePrefix + "_keggnames.tsv", "w") as fileOut:
        for key, value in keggnamesDict.items():
            fileOut.write("{0}\t{1}\n".format(key, value))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
