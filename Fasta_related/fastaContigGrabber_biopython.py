#! python3
# fastaContigGrabber
# Multifunctional fasta contig grabber program; the name is not entirely accurate
# since it not only grabs but also removes, but this name is catchy

import pyperclip, os, argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Define functions for later use
## Argument validation
def validate_args(args):
        # Check if this is running as a command-line program or in CMD window
        if args.fastaFileName == None and args.outputFileName == None:
                print('No input or output file names were provided; we will run as an iterative sequence retrieval program\n')
                return 'CMD'
        # Validate input FASTA locations
        if args.fastaFileName == None:
                print('You must specify an input FASTA file with -i field; we otherwise have no file to work with.')
                print('Fix your inputs and try again.')
                quit()
        if not os.path.isfile(args.fastaFileName):
                print('I am unable to locate the input FASTA file (' + args.fastaFileName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        fasta_or_fastq(args.fastaFileName)
        # Validate sequence ID arguments
        if args.idString == None and args.textFileName == None:
                print('You must specify at least one argument for -t OR -s fields; this will give us a list of sequence IDs to retrieve or remove.')
                print('Fix your inputs and try again.')
                quit()
        elif args.textFileName != None:
                if not os.path.isfile(args.textFileName):
                        print('I am unable to locate the input ID list file (' + args.textFileName + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Handle file overwrites
        if args.outputFileName == None:
                print('You must specify an output file name with -o field; we otherwise have no file to write to.')
                print('Fix your inputs and try again.')
                quit()
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
                quit()
        return 'command-line'

def fasta_or_fastq(fastaFile):
        # Get the first letter
        with open(fastaFile, 'r') as seqFile:
                for line in seqFile:
                        firstChar1 = line[0]
                        break
        # Check first letter to see if it conforms to fastq or fasta expected format
        if firstChar1 == '@':
                print('It looks like the input file is FASTQ format (starts with @); this program only handles FASTA, sorry. Program will exit now.')
                quit()
        elif firstChar1 != '>':
                print('I don\'t recognise the input file! It should start with "@" (FASTQ) or ">" (FASTA). Fix your inputs to continue.')
                quit()

## Fasta retrieve/remove functions
def fasta_retrieve_remove_tofile(fastaRecords, outputFileName, idList, behaviour, effort):
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['retrieve', 'remove']:
                print('fasta_retrieve_remove_tofile: Input behaviour value is not "retrieve" or "remove" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        effortPairs = {}
        # Main function
        foundList = []          # foundList will let us keep track of which entries we encountered; if we found none or just a handful of these, we'll report that later
        with open(outputFileName, 'w') as fileOut:
                for key in fastaRecords.keys():
                        record = fastaRecords[key]
                        # Find if sequence ID is in our idList
                        seqid = None
                        if record.id in idList:
                                seqid = record.id
                                foundList.append(seqid)
                        elif record.description in idList:
                                seqid = record.description
                                foundList.append(seqid)
                        # If relevant, put more effort into finding ID in our idList
                        elif effort == True:
                                idMatches = [seqid for seqid in idList if record.description.startswith(seqid)]                   # This will return all entries in the idList that partially match the current record's long name
                                cleanMatches = []
                                for match in idMatches:                                                                         # Here we begin to look through our partial idList matches
                                        if not match in fastaRecords:                                     # Firstly, we make sure this ID doesn't already have a perfect match; if it does, it belongs to that entry and we ignore it
                                                matchMatches = [seqid for seqid in fastaRecords.keys() if seqid.startswith(match)]   # Next we check to see if this ID matches more than one sequence in the FASTA; if it does, our program can't work
                                                if len(matchMatches) == 1:
                                                        cleanMatches.append(match)
                                                elif len(matchMatches) > 1:
                                                        print('fasta_retrieve_remove_tofile: sequence ID within the ID list (text file and/or string input) partially matches more than one sequence ID in your FASTA file with no exact matches.')
                                                        if len(matchMatches) < 10:
                                                                print('Debug info: ID file entry = "' + match + '", partial match list = "' + str(matchMatches) + '"; Program will exit now.')
                                                        else:
                                                                print('Debug info: ID file entry = "' + match + '", partial match list = "' + str(matchMatches[0:10]) + '... [more results hidden]"; Program will exit now.')
                                                        quit()
                                if cleanMatches != []:
                                        if len(cleanMatches) == 1:
                                                seqid = record.id                     # This is probably the best way to handle partial matches through effort; it's likely to be a shortened version of the seqid if truncated, so we won't use the long_name
                                                foundList.append(cleanMatches[0])       # We want the match ID here since that'll let us check if we found all IDs in our list later
                                                effortPairs[cleanMatches[0]] = seqid
                                        if len(cleanMatches) > 1:
                                                print('fasta_retrieve_remove_tofile: sequence IDs within the ID list (text file and/or string input) match multiple times against the same sequence ID in your FASTA file with no exact matches.')
                                                if len(cleanMatches) < 10:
                                                        print('Debug info: Sequence ID = "' + record.description + '", partial match list = "' + str(cleanMatches) + '"; Program will exit now.')
                                                else:
                                                        print('Debug info: Sequence ID = "' + record.description + '", partial match list = "' + str(cleanMatches[0:10]) + '... [more results hidden]"; Program will exit now.')
                                                quit()
                        # If we found the seqid, handle output according to behaviour
                        if seqid == None:                       # If seqid == None, we did not find it in our list
                                if behaviour == 'retrieve':     # Retrieve behaviour will only write seqids in our list
                                        continue
                                else:                           # Remove behaviour will only write seqids _not_ in our list
                                        fileOut.write('>' + record.description + '\n' + str(record.seq) + '\n')
                        else:                                   # If seqid != None, we did find it in our list
                                if behaviour == 'retrieve':     # Since we found it, we'll write it to file
                                        fileOut.write('>' + seqid + '\n' + str(record.seq) + '\n')
                                else:                           # Since we found it, we _won't_ write it to file
                                        continue
        # Report details re: foundList
        foundSet = set(foundList)
        idSet = set(idList)
        if foundSet == idSet:
                print('All values in the ID list were successfully ' + behaviour + 'd.')
        else:
                print('Not all values in the ID list were found within the FASTA file.')
                difference = max(len(foundSet), len(idSet)) - min(len(foundSet), len(idSet))
                if len(foundSet) < len(idSet):
                        print('We found ' + str(difference) + ' fewer entries than expected. These include:')
                        for entry in idSet:
                                if entry not in foundSet:
                                        print(entry)
                else:
                        print('We found ' + str(difference) + ' more entries than expected. This probably means your FASTA file has duplicate entries. If that is correct, no worries; if not, look into what might be wrong with your FASTA.')
        # Report details re: effortPairs
        if effortPairs != {}:
                print('Effort was used to find some sequences; below are listed pairings of the value in your ID list and the match found in your FASTA when effort was required')
                for key, value in effortPairs.items():
                        print(key + ' : ' + value)

## General purpose
def text_file_to_list(textFile):
        outList = []
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        outList.append(line.rstrip('\r\n'))
        return outList

## Custom fasta dictionary
def custom_to_dict(inputFileName):
        records = {}
        with open(inputFileName, "r") as fileIn:
                for line in fileIn:
                        if line.startswith(">"):
                                seqid = line.strip(">\r\n")
                                records[seqid] = SeqIO.SeqRecord(Seq(""), id=seqid)
                                prevSeqid = seqid
                        else:
                                records[prevSeqid].seq += line.strip("\r\n ")
        return records

#### USER INPUT SECTION
usage = """%(prog)s reads in FASTA files and allows various operations to extract
sequences. These include running this program via double-click or with no arguments to operate
within the CMD window allowing for ongoing retrieval of sequences, or as a command-line program
where sequence IDs are specified as arguments (-s) or as a text file listing these values (-t),
returning an output FASTA file where we 1) retrieved those sequences or 2) excluded
those sequences (dictated by -behaviour). At least one sequence ID input argument (i.e., -t
and/or -s) must be used; providing both will combine these into a single list.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fastaFileName",
               help="Input fasta file name")
p.add_argument("-t", "-textFile", dest="textFileName",
               help="Optionally specify the name of a text file listing sequence IDs")
p.add_argument("-s", "-string", dest="idString", nargs = "+",
               help="Optionally input sequence IDs as text here; entries will be separated at space characters")
p.add_argument("-e", "-effort", dest="effort", action = "store_true", default = False,
               help="Optionally put in extra effort to find ID matches; only do this if your input sequence IDs have been truncated and don't match the full sequence ID (including descriptions) in the FASTA exactly.")
p.add_argument("-b", "-behaviour", dest="behaviour", choices=['retrieve', 'remove'], default = 'retrieve',
               help="Specify program behaviour to either 'retrieve' or 'remove' provided sequences based on ID (default == retrieve)")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output fasta file name")

args = p.parse_args()
operationType = validate_args(args)

# CMD window argument input if relevant
if operationType == 'CMD':
        # Obtain data from FASTA file
        print('Enter the name of the FASTA file you wish to search. Include the file extension.')
        while True:
                try:
                        args.fastaFileName = input()
                        if os.path.isfile(args.fastaFileName) == False:
                                raise Exception
                        print('FASTA file located successfully')
                        print('')
                        break
                except KeyboardInterrupt:
                        print('KeyboardInterrupt received; program will exit now.')
                        quit()
                except:
                        print('Failed to locate the FASTA file. If you misspelt the name, try again. If the name is correct but still won\'t work, make sure the file is in the same directory as this program OR provide its full path')
                        continue

# Load the fasta file and parse its contents
try:
        records = SeqIO.to_dict(SeqIO.parse(open(args.fastaFileName, "r"), "fasta"))
except:
        records = custom_to_dict(args.fastaFileName)

# Load in text file if relevant and format our idList
textIds = set()
if args.textFileName != None:
        textIds = text_file_to_list(args.textFileName)
if args.idString == None or args.idString == '' or args.idString == []:
        args.idString = set()
idList = list(set(textIds).union(set(args.idString)))

# Strip > characters if they exist in our idList
for i in range(len(idList)):
        if idList[i].startswith('>'):
                idList[i] = idList[i].lstrip('>')       # Hopefully no one uses sequence IDs like ">>>>>>>>>1_contig"; these shouldn't be legal FASTA format IDs anyway AFAIK
# Remove blank entries if they exist in our idList
while '' in idList:
        del idList[idList.index('')]

# Extract sequences if handling command-line argument
if operationType == 'command-line':
        fasta_retrieve_remove_tofile(records, args.outputFileName, idList, args.behaviour, args.effort)
                
# Open CMD window and allow for ongoing sequence retrieval otherwise
elif operationType == 'CMD':
        rcMode = False
        print(r'Iterative sequence retrieval started; type "effort on" or "effort off" to turn on/off partial ID matching; type "rc on" or "rc off" to alternate between \r\n or \n output')
        while True:
                print('Copy a single ID or list of IDs then press ENTER')
                buttonPress = input()
                # Handle effort on/off switch
                if buttonPress.lower() == 'effort on':
                        args.effort = True
                        print('Effort mode = ON; truncated sequence IDs can be used to retrieve sequences')
                        continue
                elif buttonPress.lower() == 'effort off':
                        args.effort = False
                        print('Effort mode = OFF; exact sequence IDs must be used to retrieve sequences')
                        continue
                # Handle rc on/off switch
                if buttonPress.lower() == 'rc on':
                        rcMode = True
                        print(r'RC mode = ON; output sequences will be joined with \r\n')
                        continue
                elif buttonPress.lower() == 'rc off':
                        rcMode = False
                        print(r'RC mode = OFF; output sequences will be joined with \n')
                        continue
                # Parse the clipboard/input to format our idList
                if buttonPress == '':
                        idList = pyperclip.paste()
                else:
                        idList = buttonPress    # This program expects inputs to be present in the clipboard, but if the user typed something in then we'll assume this takes priority
                idList = idList.replace('\r', '').rstrip('\n').split('\n')
                # Strip > characters if they exist in our idList
                for i in range(len(idList)):
                        if idList[i].startswith('>'):
                                idList[i] = idList[i].lstrip('>')
                # Retrieve sequences to list
                clipboardList = []
                for entry in idList:
                        # Handle blank rows
                        if entry == '-' or entry == '':
                                print('Found a blank; skipped')
                                clipboardList.append('')
                                continue
                        # Find exact matches to names or long_names
                        elif entry in records:
                                print(entry + ' stored in clipboard.')
                                clipboardList.append(str(records[entry].seq))
                        # If relevant, put more effort into finding ID in our idList
                        elif args.effort == True:
                                partialMatches = [seqid for seqid in records.keys() if seqid.startswith(entry)]
                                if len(entry) > 100:
                                        entry = entry[0:100] + '...'    # This can help us deal with mistakes the user makes, perhaps feeding a contig sequence back into the loop which we don't want to print out
                                if len(partialMatches) == 0:
                                        print('Using effort, no sequences were found in your FASTA file with ID starting with "' + entry + '"; skipped')
                                        clipboardList.append('')
                                elif len(partialMatches) == 1:
                                        print('Found a match using effort; "' + entry + '" uniquely partially matches "' + partialMatches[0] + '"; stored in clipboard')
                                        clipboardList.append(str(records[partialMatches[0]].seq))
                                else:
                                        print('Found multiple possible matches using effort for "' + entry + '"; look through the list below to see which one you want; skipped')
                                        if len(partialMatches) < 10:
                                                print(', '.join(partialMatches))
                                        else:
                                                print(', '.join(partialMatches[0:10]) + '... [more results hidden]')
                                        clipboardList.append('')
                        # If no effort is being used, give up
                        else:
                                print('No exact matches were found in your FASTA file for "' + entry + '"; skipped')
                                clipboardList.append('')
                # Put into clipboard
                if rcMode == True:
                        pyperclip.copy('\r\n'.join(clipboardList))
                else:
                        pyperclip.copy('\n'.join(clipboardList))

# ALL DONE!
