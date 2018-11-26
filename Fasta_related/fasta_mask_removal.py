#! python3
# NCBI_contam_fix.py
# Code to automatically fix contaminants identified by NCBI during
# genome/metagenome file upload

# Import external packages
import argparse, os
from Bio import SeqIO

# Define functions for later use
# Argument validation
def validate_args(args):
        # Ensure all arguments have been specified
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument was not specified, fix this and try again.')
                        quit()
        # Validate file locations
        if not os.path.isfile(args.fastaFile):
                print('I am unable to locate the input FASTA file (' + args.fastaFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Ensure lengthCutoff value is sensible
        if args.lengthCutoff < 0:
                print('lengthCutoff value must be an integer >= 0. Fix this argument and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print('There is already a file named ' + args.outputFileName + '. Move/delete/rename this file and try again.')
                quit()

def fasta_mask_remove_withinloop(record, maskType, minLength):
        # Setup
        import re
        # Ensure maskType is sensible
        if maskType.lower() not in ['soft', 'lowerx', 'upperx']:
                print('fasta_mask_remove_withinloop: maskType value is not valid; should be "soft" or "lowerx" or "upperx"; fix the code here.')
                quit()
        # Get our sequence as a string
        seq = str(record.seq)
        # Produce our regular expression for finding masked regions
        if maskType.lower() == 'soft':
                regex = re.compile(r'[a-z]+')
        elif maskType.lower() == 'lowerx':
                regex = re.compile(r'x+')
        elif maskType.lower() == 'upperx':
                regex = re.compile(r'X+')
        # Find masked region coordinates
        coordList = []
        for maskMatch in regex.finditer(seq):
                coord = maskMatch.span()
                coordList.append(coord)
        # Handle scenario where no masked regions exist in contig
        if coordList == []:
                return False
        # Sort our coordinate list so we can trim from end to start
        coordList.sort(reverse=True)
        # Perform trimming/splitting operation
        tmpTrimList = []
        for i in range(len(coordList)):
                coord = coordList[i]
                splitSeqs = [seq[:coord[0]], seq[coord[1]:]]    # Regex matches act as 0-based ranges, so we don't need to +1 or -1 to anything
                seq = seq[:coord[0]]                            # We need to remove the end portion of the sequence now so we don't reobtain it (seq[coord[1]:] has no end boundary)
                if i+1 != len(coordList):                       # If this isn't the last coord, we just hold onto the end bit since we're going to be trimming/splitting the front bit more
                        if len(splitSeqs[-1]) >= minLength:
                                tmpTrimList.append(splitSeqs[-1])
                else:                                           # If this is the last coord (or the only coord), we want the end and start bits, but only if they meet our minLength cutoff
                        if len(splitSeqs[-1]) >= minLength:
                                tmpTrimList.append(splitSeqs[-1])
                        if len(splitSeqs[0]) >= minLength:
                                tmpTrimList.append(splitSeqs[0])
        # Handle possible scenario where we trim a contig and no bits >= minLength long result
        if tmpTrimList == []:
                return None
        # Reverse our tmpTrimList since we added sequence bits into this list from end -> front
        tmpTrimList.reverse()
        return tmpTrimList

##### USER INPUT SECTION

usage = """%(prog)s reads in a FASTA file that has been repeat masked
and removes these segments from the file. Required arguments include the type
of masking employed (i.e., "soft" == lowercase masking, lower/upper x refer to
whether the masked characters are represented by lower or uppercase x's), and
the minimum length of contigs to be returned in your FASTA file (if there are
segments inbetween or before/after masked sections that are below this length,
they will be exluded from the output file)
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fastaFile",
                  help="Input FASTA file")
p.add_argument("-m", "-mask", dest="maskType", choices=["soft", "lowerx", "upperx"],
                  help="Specify the type of masking used in your FASTA file")
p.add_argument("-l", "-length", dest="lengthCutoff", type=int,
                  help="Specify the minimum length of contigs to return in your mask-removed FASTA file")
p.add_argument("-o", "-output", dest="outputFileName",
             help="output fasta file name containing transcript sequences")

args = p.parse_args()
validate_args(args)

# Parse fasta file for iteration
records = SeqIO.parse(open(args.fastaFile, 'r'), 'fasta')

# Loop through our fasta file and generate output
with open(args.outputFileName, 'w') as fileOut:
        for record in records:
                recordSeqs = fasta_mask_remove_withinloop(record, args.maskType, args.lengthCutoff)
                if recordSeqs != None and recordSeqs != False:
                        if len(recordSeqs) == 1:
                                seqid = '>' + record.description + '_maskTrim'
                                fileOut.write(seqid + '\n' + recordSeqs[0] + '\n')
                        else:
                                for i in range(len(recordSeqs)):
                                        seqid = '>' + record.description + '_maskSplit_fragment' + str(i+1)
                                        fileOut.write(seqid + '\n' + recordSeqs[i] + '\n')
                        print('Trimmed and/or split ' + record.description)
                # Handle sequences that, after trimming, cannot have bits >= 200bp in length
                elif recordSeqs == None:
                        print('Excluded ' + record.description + ' after mask removal reduced length to <' + str(args.lengthCutoff))
                        continue
                # Output untouched sequences
                else:
                        fileOut.write('>' + record.description + '\n' + str(record.seq) + '\n')

# All done
print('Program completed successfully!')