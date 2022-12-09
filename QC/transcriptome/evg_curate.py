#! python 3
# evg_curate.py
# A script/program designed to curate a transcriptome assembly
# produced by the EvidentialGene pipeline (mine) and reduce
# the size of the transcriptome by filtering different
# EVG classes

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_BlastIO, domtblout_handling

VALID_EVG_CLASSES = [
    "althi", "althi1", "althinc",
    "altmfrag", "altmid", "main",
    "mainnc", "noclass", "noclassnc",
    "parthi", "parthi1", "perfdupl",
    "perffrag", "smallorf"
]

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.transcriptomeFile):
        print(f'I am unable to locate the input transcriptome file ({args.transcriptomeFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.blastFile != None:
        if not os.path.isfile(args.blastFile):
            print(f'I am unable to locate the input BLAST outfmt6 file ({args.blastFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.hmmerFile != None:
        if not os.path.isfile(args.hmmerFile):
            print(f'I am unable to locate the input HMMER domtblout file ({args.hmmerFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate EVG class input
    for evgClass in args.dropClasses:
        if not evgClass in VALID_EVG_CLASSES:
            print(f"{evgClass} is not recognised as a member of the valid EVG classes {VALID_EVG_CLASSES}")
            print("Make sure your input is lowercase, or exists in the list above then try again.")
            quit()
    for evgClass in args.filterClasses:
        if not evgClass in VALID_EVG_CLASSES:
            print(f"{evgClass} is not recognised as a member of the valid EVG classes {VALID_EVG_CLASSES}")
            print("Make sure your input is lowercase, or exists in the list above then try again.")
            quit()
    if len(set(args.dropClasses).intersection(set(args.filterClasses))) != 0:
        print("Detected identical class(es) in the drop and filter groups.")
        print("This would mean we will just drop the class, and filtering becomes irrelevant.")
        print("But this isn't what you probably meant to happen, so I assume it's an error on your part.")
        print("Clarify what you intend by putting each class into only 1 group, and try again.")
        quit()
    # Validate logical argument combinations
    if args.dropClasses == [] and args.filterClasses == []:
        print("At least one of the -d or -f arguments must be given values.")
        print("Otherwise, there's nothing for me to even do!")
        print("Fix this and try again.")
        quit()
    if args.filterClasses != []:
        if args.blastFile == None and args.hmmerFile == None:
            print("If you specify -f filter classes, at least one of the " 
                  "--blast or --hmmer arguments must be given values.")
            print("Otherwise, there's nothing for me to filter on.")
            print("Fix this and try again.")
            quit()
    # Validate numeric arguments
    if args.blastEvalue < 0:
        print("blastEvalue must be number >= 0")
        quit()
    if args.hmmerEvalue < 0:
        print("hmmerEvalue must be number >= 0")
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print('The specified output file already exists. This program will not allowing overwriting.')
        print('Specify a different file name, or move the existing one and try again.')
        quit()

def curate_evg_transcriptome(transcriptomeRecords, blastDict, hmmerDict,
                             blastEvalue, hmmerEvalue, dropClasses,
                             filterClasses, outputFileName):
    '''
    Parameters:
        transcriptomeRecords -- a Bio.SeqIO.parse() generator object
        blastDict -- a dictionary object obtained by parsing an outfmt6 file
                     by ZS_BlastIO.BLAST_Results (taking its .results attribute).
                     It will have a structure like:
                     {
                         'queryID1': [tid, identityPct, qstart, qend, tstart, tend, evalue, bitscore],
                         'queryID2': [ ... ],
                         ...
                     }
        hmmerDict -- a dictionary object obtained by parsing a domtblout file
                     by domtblout_handling.hmmer_parse; note that we assume the
                     domain models are protein ones.
                     It will have a structure like:
                     {
                         'queryID1': [did, dstart, dend, evalue, hmmFrom, hmmTo],
                         'queryID2': [ ... ],
                         ...
                     }
        blastEvalue -- a float value indicating the threshold to enforce when seeing
                       if a sequence we want to filter has a significant BLAST hit
        hmmerEvalue -- a float value indicating the threshold to enforce when seeing
                       if a sequence we want to filter has a significant HMMER hit
        dropClasses -- a list containing valid EVG classes which we will entirely drop
        filterClasses -- a list containing valid EVG classes that we want to filter to
                         retain ONLY sequences with a significant BLAST and/or HMMER hit
        outputFileName -- a string indicating the location to write the filtered
                          transcriptome file to
    '''
    with open(outputFileName, "w") as fileOut:
        for record in transcriptomeRecords:
            # Get EVG class tags
            tags = record.description.split(" ", maxsplit=1)[1]
            tags = {
                t.strip(" ").split("=")[0]:t.strip(" ").split("=")[1]
                    for t in tags.rstrip(";").split("; ")
            }
            evgClass = tags["evgclass"].split(",")[0]
            
            # Check if this class is in our dropset
            if evgClass in dropClasses:
                continue
            
            # Skip any downstream effort if we have no filter classes
            if filterClasses == []:
                continue
            
            # Check for significant BLAST hit
            hasBlastHit = False
            if record.id in blastDict:
                _, _, _, _, _, _, evalue, _ = blastDict[record.id]
                if evalue <= blastEvalue:
                    hasBlastHit = True
            
            # If no significant BLAST hit, check for significant HMMER hit
            hasHmmerHit = False
            if not hasBlastHit:
                if record.id in blastDict:
                    _, _, _, evalue, _, _ = hmmerDict[record.id]
                    if evalue <= hmmerEvalue:
                        hasHmmerHit = True
            
            # If we found a hit, write to file
            if hasBlastHit or hasHmmerHit:
                fileOut.write(">{0}\n{1}\n".format(
                    record.description,
                    str(record.seq)
                ))

def main():
    usage = """%(prog)s receives a transcriptome created by the EvidentialGene pipeline
    and curates it to remove low-quality sequences from specific EVG classes. Ideally,
    you run this program with BLAST and HMMER results to enable fine-tuned curation
    of classes on the basis of whether they receive significant hits or not.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="transcriptomeFile",
                   required=True,
                   help="Specify the location of the input EVG transcriptome file")
    p.add_argument("-o", dest="outputFileName", required=True,
                   help="Specify where the curated transcriptome should be written")
    # Opts (but at least 1 must be provided)
    p.add_argument("-d", dest="dropClasses", nargs="+",
                   required=False,
                   help="Optionally, specify one or more EVG classes to drop",
                   default=[])
    p.add_argument("-f", dest="filterClasses", nargs="+",
                   required=False,
                   help="Optionally, specify one or more EVG classes to filter",
                   default=[])
    # Opts
    p.add_argument("--blast", dest="blastFile",
                   required=False,
                   help="Optionally, specify a BLAST outfmt6 file")
    p.add_argument("--blast-evalue", dest="blastEvalue", type=float,
                   required=False,
                   help="""Optionally, specify an E-value to use for using BLAST
                   hits (default == 1e-5)""",
                   default=1e-5)
    p.add_argument("--hmmer", dest="hmmerFile",
                   required=False,
                   help="Optionally, specify a HMMER domtblout file")
    p.add_argument("--hmmer-evalue", dest="hmmerEvalue", type=float,
                   required=False,
                   help="""Optionally, specify an E-value to use for using HMMER
                   hits (default == 1e-3)""",
                   default=1e-3)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse BLAST file (if relevant)
    if args.blastFile != None:
        blastDict = ZS_BlastIO.BLAST_Results(args.blastFile).results
    else:
        blastDict = {}
    
    # Parse HMMER file (if relevant)
    if args.blastFile != None:
        hmmerDict = domtblout_handling.hmmer_parse(args.hmmerFile)
    else:
        hmmerDict = {}
    
    # Load transcriptome records
    records = SeqIO.parse(open(args.transcriptomeFile, "r"), "fasta")
    
    # Curate EVG transcriptome
    curate_evg_transcriptome(records, blastDict, hmmerDict,
                             args.blastEvalue, args.hmmerEvalue, args.dropClasses,
                             args.filterClasses, args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
