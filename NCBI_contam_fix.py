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
        # Validate file locations
        if not os.path.isfile(args.fastaFile):
                print('I am unable to locate the input FASTA file (' + args.fastaFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.textFile):
                print('I am unable to locate the input text contaminants file (' + args.textFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print('There is already a file named ' + args.outputFileName + '. Move/delete/rename this file and try again.')
                quit()

def ncbi_contam_parse(ncbiTextFile):
        # Set up
        import re
        categoryLabels = ('Exclude:', 'Duplicated:', 'Trim:')   # These correspond to the type of edits NCBI specified and this code will handle
        regex = re.compile('(' + '|'.join(categoryLabels) + r')\n.+?\n(.+?)\n\n', re.DOTALL)
        excludeList = []        # Duplicate entries for removal will go into this list
        trimDict = {}           # Need to use a dictionary structure since there might be multiple regions in a single contig
        # Read in file
        ncbiText = open(ncbiTextFile, 'r').read()
        # Pull out edit requirements with regex
        requirementList = regex.findall(ncbiText)
        if requirementList == []:
                print('Wasn\'t able to find the expected text in your contaminants file. Did you open this file in Excel or another program on a Windows system and make edits to it?')
                print('If so, you need to re-download the original file from NCBI and make sure to only make edits in an editor which uses \\n rather than \\r\\n for lines.')
                print('Program will exit now.')
                quit()
        # Organise requirements into lists of action types
        for entry in requirementList:
                if entry[0] == 'Exclude:':
                        # Split section into lines
                        table = entry[1].split('\n')
                        # Extract scaffold ID from each line and put into excludeList
                        for line in table:
                                sl = line.split('\t')
                                excludeList.append(sl[0])
                elif entry[0] == 'Duplicated:':
                        # Split section into lines
                        table = entry[1].split('\n')
                        # Add all but the first value in each line to our excludeList
                        for line in table:
                                sl = line.rsplit(' (', maxsplit=1)[0].split(' ') # The first rsplit removes the contig length e.g., (301 bp), the second will split sequence IDs
                                for i in range(len(sl)):
                                        if i != 0:
                                                if 'lcl|' in sl[i]:     
                                                        excludeList.append(sl[i][4:])   # This lcl| prefix seems to be attached to the sequence IDs I'm seeing in this NCBI report
                                                else:
                                                        excludeList.append(sl[i])       # In case it doesn't always happen I'll handle it specifically
                elif entry[0] == 'Trim:':
                        # Split section into lines
                        table = entry[1].split('\n')
                        # Loop through each line and associate contig IDs with the region for trimming
                        for line in table:
                                sl = line.split('\t')
                                # Split multiple coord spans into individual ones for handling
                                spans = sl[2].split(',')
                                for span in spans:
                                        coords = list(map(int, span.split('..')))
                                        if sl[0] not in trimDict:
                                                trimDict[sl[0]] = [coords]
                                        else:
                                                trimDict[sl[0]].append(coords)
        return excludeList, trimDict

def fasta_trim_split(fastaFile, trimDict):      # I'm not doing it this way anymore; however, it could be useful in the future so I'm leaving it here
        # Set up
        from Bio import SeqIO
        trimmedSeqs = {}
        # Load in fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Loop through fasta until we find something in our trimDict
        for record in records:
                if record.id in trimDict:
                        seqid = record.id
                elif record.description in trimDict:
                        seqid = record.description
                else:
                        continue
                # Get our sequence as a string
                seq = str(record.seq)
                # Sort our coordinate list so we can trim from end to start
                coordList = trimDict[seqid]
                coordList.sort(reverse=True)
                # Begin trimming/splitting process
                tmpTrimList = []
                for i in range(len(coordList)):
                        coord = coordList[i]
                        splitSeqs = [seq[:coord[0]-1], seq[coord[1]:]]  # -1 to first coord since (I assume) the NCBI ranges are 1-based; second one is fine cause of how 1-based->0-base works and stuff
                        seq = seq[:coord[0]-1]                          # We need to remove the end portion of the sequence now so we don't reobtain it (seq[coord[1]:] has no end boundary)
                        if i+1 != len(coordList):                       # If this isn't the last coord, we just hold onto the end bit since we're going to be trimming/splitting the front bit more
                                if splitSeqs[-1] != '':
                                        tmpTrimList.append(splitSeqs[-1].lower().strip('n').upper())    # This satisfies NCBI's requirement that stretches of N be removed from ends of split bits
                        else:                                           # If this is the last coord (or the only coord), we want the end and start bits, but only if they're not empty
                                if splitSeqs[-1] != '':                 # They'll be empty if we're trimming the start/end bits of the sequence and not something internal which would result in splitting
                                        tmpTrimList.append(splitSeqs[-1].lower().strip('n').upper())
                                if splitSeqs[0] != '':
                                        tmpTrimList.append(splitSeqs[0].lower().strip('n').upper())
                # Reverse our tmpTrimList since we added sequence bits into this list from end -> front, then store in our dict
                tmpTrimList.reverse()
                trimmedSeqs[record.description] = tmpTrimList
        return trimmedSeqs

def fasta_trim_split_withinloop(record, trimDict):
        # Find out if this sequence is for splitting
        if record.id in trimDict:
                seqid = record.id
        elif record.description in trimDict:
                seqid = record.description
        else:
                return False
        # Get our sequence as a string
        seq = str(record.seq)
        # Sort our coordinate list so we can trim from end to start
        coordList = trimDict[seqid]
        coordList.sort(reverse=True)
        # Begin trimming/splitting process
        tmpTrimList = []
        for i in range(len(coordList)):
                coord = coordList[i]
                splitSeqs = [seq[:coord[0]-1], seq[coord[1]:]]  # -1 to first coord since (I assume) the NCBI ranges are 1-based; second one is fine cause of how 1-based->0-base works and stuff
                seq = seq[:coord[0]-1]                          # We need to remove the end portion of the sequence now so we don't reobtain it (seq[coord[1]:] has no end boundary)
                if i+1 != len(coordList):                       # If this isn't the last coord, we just hold onto the end bit since we're going to be trimming/splitting the front bit more
                        if len(splitSeqs[-1].upper().strip('N')) >= 200:                # This satisfies NCBI's requirement that contigs be >= 200bp in length and handles empty bits ('') that can result from trimming start/end bits
                                tmpTrimList.append(splitSeqs[-1].upper().strip('N'))    # The upper.strip stuff satisfies NCBI's requirement that stretches of N be removed from ends of split bits
                else:                                           # If this is the last coord (or the only coord), we want the end and start bits, but only if they're not empty
                        if len(splitSeqs[-1].upper().strip('N')) >= 200:
                                tmpTrimList.append(splitSeqs[-1].upper().strip('N'))
                        if len(splitSeqs[0].upper().strip('N')) >= 200:
                                tmpTrimList.append(splitSeqs[0].upper().strip('N'))
        # Handle possible scenario where we trim a contig and no bits >= 200bp long result
        if tmpTrimList == []:
                return None
        # Reverse our tmpTrimList since we added sequence bits into this list from end -> front, then store in our dict
        tmpTrimList.reverse()
        return tmpTrimList

##### USER INPUT SECTION

usage = """%(prog)s reads in a genome fasta file in addition to an NCBI
contaminants text file and automatically fixes these errors in your file
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fastaFile",
                  help="genome fasta file")
p.add_argument("-t", "-text", dest="textFile",
                  help="NCBI contaminants text file")
p.add_argument("-o", "-output", dest="outputFileName",
             help="output fasta file name containing transcript sequences")

args = p.parse_args()
validate_args(args)

# Parse NCBI file to identify sequences for removal/trimming
excludeList, trimDict = ncbi_contam_parse(args.textFile)
# Trim sequences and store in dict
#trimmedSeqs = fasta_trim_split(args.fastaFile, trimDict)

# Parse fasta file for iteration
records = SeqIO.parse(open(args.fastaFile, 'r'), 'fasta')

# Loop through our fasta file and generate output
with open(args.outputFileName, 'w') as fileOut:
        for record in records:
                # Exclude with reference to excludeList
                if record.id in excludeList or record.description in excludeList:
                        print('Excluded ' + record.description)
                        continue
                # Handle trimmed sequence output
                recordSeqs = fasta_trim_split_withinloop(record, trimDict)
                if recordSeqs != False:
                        if len(recordSeqs) == 1:
                                seqid = '>' + record.description + '_ncbiTrim'
                                fileOut.write(seqid + '\n' + recordSeqs[0] + '\n')
                        else:
                                for i in range(len(recordSeqs)):
                                        seqid = '>' + record.description + '_ncbiSplit_fragment' + str(i+1)
                                        fileOut.write(seqid + '\n' + recordSeqs[i] + '\n')
                        print('Trimmed ' + record.description)
                # Output untouched sequences
                else:
                        fileOut.write(record.description + '\n' + str(record.seq) + '\n')

# All done
print('Program completed successfully!')