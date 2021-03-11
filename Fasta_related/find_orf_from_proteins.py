#! python3
# find_orf_from_proteins.py

# Load packages
import re, os, argparse
from Bio import SeqIO

# ### USER INPUT
# usage = """%(prog)s reads in a FASTA formatted transcriptome file and another FASTA file with proteins. It will search
# through the transcriptome for matching CDS regions and output these as a new FASTA file. Note that this script assumes the 
# sequence IDs are roughly equivalent between files, so it doesn't perform a deep search. Maybe I'll make that another time.
# """
# # Reqs
# p = argparse.ArgumentParser(description=usage)
# p.add_argument("-t", dest="transcriptomeFile",
#     help="Input transcriptome fasta file name")
# p.add_argument("-t", dest="proteinFile",
#     help="Input protein fasta file name")
# p.add_argument("-o", "-output", dest="outputFileName",
#     help="Output fasta file name")
# args = p.parse_args()

## HARD-CODED TESTING
class Args:
    def __init__(self, t, p, o):
        self.transcriptomeFile = t
        self.proteinFile = p
        self.outputFileName = o

args = Args(r"F:\toxins_annot\proteomes\nts\Telmatactis_transcriptome_nt.fasta", r"F:\toxins_annot\proteomes\telmatactis_secretome.fasta", "Telmatactis_CDSs_nt.fasta")

# Load the protein fasta file as a generator objects
proteinRecords = SeqIO.parse(open(args.proteinFile, 'r'), 'fasta')

# Get the transcriptome as a dictionary
transcriptomeDict = SeqIO.to_dict(SeqIO.parse(open(args.transcriptomeFile, 'r'), 'fasta'))

# Get the nucleotide (record) out of our generator (records) and grab them ORFs!
with open(args.outputFileName, "w") as fileOut:
    for record in proteinRecords:
        # Get relevant details
        protein = str(record.seq)
        seqid = record.description
        # Handle TransDecoder IDs
        if seqid.rsplit("|", maxsplit=1)[-1].startswith("CDS") or seqid.rsplit("|", maxsplit=1)[-1].startswith("Telmatactis"):
            seqid = seqid.rsplit("|", maxsplit=1)[0]
        # Find the transcript record
        transcriptomeRecord = transcriptomeDict[seqid]
        transcript = str(transcriptomeRecord.seq)
        # Find frame where translation was obtained
        found = False
        for strand, nuc in [(+1, transcriptomeRecord.seq), (-1, transcriptomeRecord.seq.reverse_complement())]:
            for frame in range(3):
                length = 3 * ((len(transcriptomeRecord)-frame) // 3)
                frameNuc = str(nuc[frame:frame+length])
                frameProt = str(nuc[frame:frame+length].translate(table=1))
                if protein[1:] in frameProt: # find the frame by a fragment, in cases where the start codon has been altered
                    found = True
                    break
            if found == True:
                break
        if found == False:
            stophere
        # Split protein/nucleotide into corresponding ORFs
        ongoingLength = 0                                       # The ongoingLength will track where we are along the unresolvedProt sequence for getting the nucleotide sequence
        splitNucleotide = []
        splitProtein = []
        frameProt = frameProt.split('*')
        for i in range(len(frameProt)):
                if len(frameProt) == 1 or i + 1 == len(frameProt):    # This means the splitProtein has no stop codons or we're looking at the last ORF which also has no stop codon
                        splitProtein.append(frameProt[i])
                        splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i])*3])      # This will grab the corresponding nucleotide region
                        ongoingLength += len(frameProt[i])*3
                else:
                        splitProtein.append(frameProt[i] + '*')
                        splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i] + '*')*3])       
                        ongoingLength += (len(frameProt[i]) + 1)*3
        # Find the index for the split segment
        for splitIndex in range(len(splitProtein)):
            if protein[1:] in splitProtein[splitIndex]: # [1:] for start codon alteration
                break
        nuclSegment = splitNucleotide[splitIndex]
        protSegment = splitProtein[splitIndex]
        # Find the start site
        startCodon = protSegment.index(protein[1:]) - 1 # -1 to account for [1:]
        nuclSegment = nuclSegment[startCodon*3:]
        # Ensure equality
        assert len(protSegment[startCodon:]) + 1 >= len(protein) # +1 incase the original protein didn't have an asterisk stop
        # Output to file
        fileOut.write(">{0}\n{1}\n".format(record.description, nuclSegment))

#### SCRIPT ALL DONE
print('Program completed successfully!')
