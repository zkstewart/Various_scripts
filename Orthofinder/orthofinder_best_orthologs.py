#! python3
# orthofinder_group_containing_species
# Script to produce an abbreviated Orthogroups.csv file which contains only
# orthogroups which contain the species of interest

# Load packages
import os, argparse, sys, math, platform, re

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_BlastIO, ZS_SeqIO, ZS_AlignIO
from OrthoGroups import OrthoGroups, fastaDict_formatter

def validate_args(args):
    # Validate input locations
    if not os.path.isfile(args.orthogroupsTSV):
        print('I am unable to locate the Orthogroups file (' + args.orthogroupsTSV + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for fastaDir in args.fastaDirs:
        if not os.path.isdir(fastaDir):
            print('I am unable to locate the FASTA directory (' + fastaDir + ')')
            print('Make sure you\'ve typed the directory name or location correctly and try again.')
            quit()
    if platform.system() == "Windows":
        if not os.path.isfile(os.path.join(args.mafftDir, "mafft.bat")):
            print("{0} does not exist".format(os.path.join(args.mafftDir, "mafft.bat")))
            print("Make sure you've indicated the correct MAFFT directory and try again.")
            quit()
    else:
        if not os.path.isfile(os.path.join(args.mafftDir, "mafft")) and not os.path.isfile(os.path.join(args.mafftDir, "mafft.exe")):
            print("mafft or mafft.exe does not exist at {0}".format(args.mafftDir))
            print("Make sure you've indicated the correct MAFFT directory and try again.")
            quit()
    
    # Clean up suffixes if needed
    if any([ "." in suffix for suffix in args.suffixes ]):
        args.suffixes = [ suffix.strip(".") for suffix in args.suffixes ]
        print("FYI: I've stripped periods from the --suffixes values since they're not wanted here.")
    
    # Validate numeric inputs
    if args.threads < 1:
        print("--threads must be a positive integer; increase this value and try again.")
        quit()
    if args.evalue < 0:
        print("--evalue must be a positive float; increase this value and try again.")
        quit()
    
    # Handle file overwrites
    if os.path.exists(args.outputDirectory) and not os.path.isdir(args.outputDirectory):
        print(f"{args.outputDirectory} already exists, but isn't a directory? I need to write to a directory...")
        print("Delete/move whatever exists here, or specify a different output directory next time.")
        quit()
    elif not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Created '{args.outputDirectory}' as part of argument validation.")
    else:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")

def extract_tsv_header(orthogroupsTSV):
    with open(orthogroupsTSV, 'r') as fileIn:
        header = fileIn.readline().strip().split('\t')
    return header[1:]

def add_seqs_to_fasta(fastASeq_list):
    FASTA_obj = ZS_SeqIO.FASTA(None)
    for fastaASeq_obj in fastASeq_list:
        FASTA_obj.add(fastaASeq_obj)
    return FASTA_obj

def main():
    ### USER INPUT
    usage = """%(prog)s will parse the Orthogroups.tsv file, before or after filtration,
    and attempt to identify the most likely 1-to-1 orthologs from among orthogroups
    containing any species of interest (SOI).
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="orthogroupsTSV",
                   required=True,
                   help="Specify the orthogroup TSV file")
    p.add_argument("-d", dest="fastaDirs",
                   required=True,
                   nargs="+",
                   help="Specify the director(y/ies) to locate FASTA files within")
    p.add_argument("-m", dest="mafftDir",
                   required=True,
                   help="Specify the directory to locate the mafft.exe or mafft.bat file within")
    p.add_argument("-o", dest="outputDirectory",
                   help="Output location for MSA files")
    p.add_argument("-s", dest="soi",
                   required=True,
                   nargs="+",
                   help="Column name(s) of the species of interest (SOI)")
    
    p.add_argument("--suffixes", dest="suffixes",
                   required=True,
                   nargs="+",
                   help="""Specify the file suffixes to look for when parsing FASTA
                   directories; this prevents issues with duplicate file finding""")
    # Optionals
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="Specify the number of BLAST threads to run; default==1",
                   default=1)
    p.add_argument("--evalue", dest="evalue",
                   required=False,
                   type=float,
                   help="Specify the E-value threshold; default==1e-3",
                   default=1e-3)
    
    args = p.parse_args()
    validate_args(args)
    
    # Extract header from TSV file
    header = extract_tsv_header(args.orthogroupsTSV)
    
    # Make sure SOIs are in the header
    for soi in args.soi:
        if soi not in header:
            print(f"SOI '{soi}' not found in OrthoGroups TSV file header.")
            print('Make sure your SOI is within the file header and try again; header shown below.')
            print(header)
            quit()
    
    # Format dictionary of species and their FASTA files
    fastaDict = fastaDict_formatter(args.fastaDirs, header, args.suffixes)
    
    # Parse OrthoFinder file
    og = OrthoGroups(args.orthogroupsTSV, fastaDict)
    
    # Find best orthologs for each sequence in each SOI
    SPLIT_STRING = "__zks__"
    
    orthologsDict = {}
    for soi in args.soi:
        groupsWithSoi = list(og.species[soi].keys())
        
        # Iterate through each orthogroup
        for orthogroupID in groupsWithSoi:
            orthologsDict.setdefault(soi, {})
            
            # Create ZS_SeqIO.FASTA object for query and target sequences
            soiSequences = [
                ZS_SeqIO.FastASeq(id, str(og.sequences[soi][id].seq))
                for id in og.species[soi][orthogroupID]
            ]
            soiFASTA_obj = add_seqs_to_fasta(soiSequences)
            
            targetSequences = [
                ZS_SeqIO.FastASeq(f"{species}{SPLIT_STRING}{id}", str(og.sequences[species][id].seq))
                for species, seqIDs in og.groups[orthogroupID].items()
                for id in seqIDs
                if species != soi
            ]
            targetFASTA_obj = add_seqs_to_fasta(targetSequences)
            
            # Run BLAST
            blaster = ZS_BlastIO.BLAST(soiFASTA_obj, targetFASTA_obj, "blastn")
            blaster.set_threads(args.threads)
            blaster.set_evalue(args.evalue)
            blastDict, _ = blaster.get_blast_results()
            
            # Get the best hit from each target species
            bestOrthologs = { queryID: {} for queryID in blastDict.keys() }
            for queryID in bestOrthologs.keys():
                for species in header:
                    if species != soi:
                        for result in blastDict[queryID]:
                            resultSpecies = result[0].split(SPLIT_STRING, maxsplit=1)[0]
                            if resultSpecies not in bestOrthologs[queryID]:
                                bestOrthologs[queryID][resultSpecies] = result
            
            # Store result
            orthologsDict[soi].update(bestOrthologs)
    
    # Generate MSAs for each ortholog group
    aligner = ZS_AlignIO.MAFFT(args.mafftDir)
    aligner.set_threads(args.threads)
    
    for soi, orthologs in orthologsDict.items():
        for queryID, bestOrthologs in orthologs.items():
            # Figure out an E-value cut-off to remove worse hits
            bestResult = min([ result[6] for result in bestOrthologs.values() ]) # result[6] == E-value
            if bestResult == 0:
                bestResult = 1e-323 # can't go lower, is rounded to 0 otherwise
            
            power = math.floor(math.log10(bestResult))
            powerCutoff = power / 2
            
            # Begin making the ZS_SeqIO.FASTA object
            FASTA_obj = add_seqs_to_fasta([ZS_SeqIO.FastASeq(queryID, str(og.sequences[soi][queryID].seq))])
            
            # Loop through each target species' best ortholog
            for targetSpecies, resultsList in bestOrthologs.items():
                tID, identity, qstart, qend, tstart, tend, evalue = resultsList
                
                # Skip if this hit sucks (by comparison)
                if evalue > 10**powerCutoff:
                    continue
                
                # Get target sequence in the correct orientation
                tID = tID.split(SPLIT_STRING, maxsplit=1)[1]
                targetSequence = str(og.sequences[targetSpecies][tID].seq)
                if tstart > tend: # +ve orientation
                    targetSequence = ZS_SeqIO.FastASeq.get_reverse_complement(None, targetSequence)
                
                # Store it in FASTA object
                FASTA_obj.add(ZS_SeqIO.FastASeq(tID, targetSequence))
            
            # Align with MAFFT
            aligner.run(FASTA_obj)
            
            # Write output
            outputFileName = os.path.join(args.outputDirectory, f"{queryID}.fasta")
            if os.path.isfile(outputFileName):
                print(f"Skipping '{outputFileName}' since it already exists.")
            else:
                FASTA_obj.write(outputFileName, asAligned = True)
    
    # Generate a file for the manual curation process
    orderedSeqIDs = sorted([ y for x in orthologsDict.values() for y in x.keys() ],
                           key = lambda x: int(re.sub(r"[^\d]+", "", x)))
    orderedSOIs = [
        soi
        for seqID in orderedSeqIDs
        for soi in orthologsDict.keys()
        if seqID in orthologsDict[soi]
    ]
    assert len(orderedSeqIDs) == len(orderedSOIs), "Mismatch between ordered sequence IDs and SOIs"
    
    with open(os.path.join(args.outputDirectory, "manual_curation.tsv"), "w") as fileOut:
        fileOut.write("Sequence_ID\tSOI\trepresentative\n")
        for seqID, soi in zip(orderedSeqIDs, orderedSOIs):
            fileOut.write(f"{seqID}\t{soi}\t\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
