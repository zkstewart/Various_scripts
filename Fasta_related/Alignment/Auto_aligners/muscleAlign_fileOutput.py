import os, re, pyperclip, textwrap
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO


# Global values
transIDs = []
aligns = []
tempFirstAA = []
chunkedList = []
records = []

asteriskHead = ''
alignsForPyperclip = ''

mode = 1

# Create regex and add character so regex can function
stdoutRegex = re.compile(r'(>[^>].+?)\n([RHKDESTNQCUGPAILMFWYVX\n-]+?)>', re.DOTALL)


# Definition for string chunker function
def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

# Default file location. Will need to be changed by the user if they wish to use a different folder to house their muscle .exe
muscle_exe = r"D:\Bioinformatics\Protein_analysis\muscle"
#in_file = r"C:\Users\lythl\Desktop\DGE\HMMER and BLAST+\final domain search\temporary_unaligned_1234567890.fasta"

# Prompt user to call for sequences
print('Do you want to refer to a fasta file, or something in your clipboard? Enter 1 for fasta, or 2 for clipboard')
print('If you want to change this option after this, you\'ll need to re-run this script')
while True:
        clipOrFasta = input()
        if clipOrFasta != '1' and clipOrFasta != '2':
                print('You mistyped something. Enter 1 for fasta or 2 for clipboard then press ENTER')
        else:
                break

# Depending on above, retrieve sequences
while True:
    try:
        if clipOrFasta == '1':
            print('Enter the name of the fasta file. Do not include the file extension.')
            while True:
                    try:
                            fileName = input()
                            if os.path.isfile(fileName + '.fasta') == False:
                                    raise Exception
                            print('Fasta file located successfully')
                            print('')
                            break
                    except:
                            print('Fasta file failed to load. If you misspelt the name, try again. If this script isn\'t in the same directory as this file, move the script there then re-run this script')
                            continue

            # Format the path to fasta file to feed into muscle
            unalignedFasta = os.getcwd() + '\\' + fileName + '.fasta'
            muscle_cline = MuscleCommandline(muscle_exe, input = unalignedFasta)
            stdout, stderr = muscle_cline()
            
        else:
            print('Copy the sequences you wish to align from within an Excel environment then press ENTER. This script operates in two modes: one accepts transcript ID + sequence, the other accepts sequences alone. Enter \'123\' to switch between modes.')
            buttonPress = input()
            if buttonPress == '123':        # This will tell the script to swap between modes
                if mode == 1:
                    mode = 2
                    print('Changed mode to accept raw sequences. Copy your sequence now.')
                    buttonPress = input()
                else:
                    mode = 1
                    print('Changed mode to accept IDs and sequences. Copy your sequence now.')
                    buttonPress = input()
            rawText = pyperclip.paste()
            newLineSplit = rawText.split(sep='\r\n')
            del newLineSplit[-1]                                                            # There's always a blank row at the end when splitting by new line with Excel

            # Search through each row for sequences to be aligned
            for entry in newLineSplit:
                tempSplit = entry.split(sep='\t')
                if mode == 1:
                    for i in range(0, len(tempSplit), 2):
                        records.append(SeqRecord(Seq(tempSplit[i+1]), id= '>' + tempSplit[i], description = ""))           # Creates a SeqRecord object for each sequence with just the sequence and its transcript ID
                else:
                    for i in range(0, len(tempSplit)):
                        records.append(SeqRecord(Seq(tempSplit[i]), id= '>sequence' + str(i), description = ""))           # Creates a SeqRecord object for each sequence with just the sequence and an invented transcript ID
            # Write data to RAM to feed into muscle. Note that I don't fully understand how this works, took it from the biopython manual
            handle = StringIO()
            SeqIO.write(records, handle, "fasta")
            data = handle.getvalue()

            # Feed this data into muscle
            muscle_cline = MuscleCommandline(muscle_exe)
            stdout, stderr = muscle_cline(stdin=data)

        # Use regex to capture sequences and format each sequence
        stdout = stdout + '>'                                       # Need a '>' at the end of stdout so the regex can capture each sequence accurately
        muscleOutAligns = stdoutRegex.findall(stdout)
        for entry in muscleOutAligns:
            muscleAlign = entry[1].replace('\n', '')
            transIDs.append(entry[0])
            aligns.append(muscleAlign)

        # Format the asterisk header similar to MEGA alignment
        for x in range(len(aligns[0])):                             # All hits in aligns should be of equal length, choosing the first is the safest option
            for i in range(len(aligns)):             
                tempFirstAA.append(aligns[i][x])                    # i corresponds to the sequence, x corresponds to the position in sequence  
            setTempFirstAA = set(tempFirstAA)
            if len(setTempFirstAA) == 1:                            # If they're all the same, calling set() will return just that one letter.
                asteriskHead += '*'
            else:
                asteriskHead += ' '                                 # Otherwise there will be multiple entries, and thus a blank                            
            tempFirstAA = []                                        # Clear this from memory to restart for next loop

        # Produce readable output within the cmd window
        transLength = len(max(transIDs, key=len))                   # Get the longest transcript ID so we can justify the other ones to make them look good when printed in the cmd window
        transIDs = [' '.ljust(transLength)] + transIDs              # Add a blank line to this list to justify the asterisk header
        aligns = [asteriskHead] + aligns                            # Add the asterisk head to the aligns list 

        for entry in aligns:
                chunkedList.append(list(chunkstring(entry, 90)))    # Here we chunk the sequence into a shorter length so it can fit within a relatively narrow cmd window

        # Print the function to the cmd window
        for x in range(len(chunkedList[0])):
            for i in range(len(chunkedList)):
                print(transIDs[i].ljust(transLength) + '\t' + chunkedList[i][x])       # Calling it like this ensures the transcript ID is shown on each line before the alignment, and everything stays in its correct orientation

        # Produce output to put into clipboard if closer inspection is wanted
        for i in range(1, len(aligns)):                    
            alignsForPyperclip += transIDs[i] + '\n' + aligns[i] + '\n'            # Get all the alignments into something that can be called in a single function. We +1 to the transIDs call since we added a blank line previously
        outfile = open('muscle_out.fasta', 'w')
        outfile.write(alignsForPyperclip)
        outfile.close()

        # Let user know the loop is done, and prompt for beginning again
        print('All done! If you want to do this again, go ahead!')
        print('')

        # Dump all temporary values before starting the loop again
        transIDs = []
        aligns = []
        chunkedList = []
        records = []
        asteriskHead = ''
        alignsForPyperclip = ''
    except:
        print('Something went wrong. You might have had something highlighted in the cmd window, or perhaps the section you copied or the fasta file you specified is in the wrong format. Try again.')
        # Dump all temporary values before starting the loop again
        transIDs = []
        aligns = []
        chunkedList = []
        records = []
        asteriskHead = ''
        alignsForPyperclip = ''
        print('')
