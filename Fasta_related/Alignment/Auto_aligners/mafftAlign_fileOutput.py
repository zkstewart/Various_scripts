import os, re, pyperclip, textwrap, shutil, subprocess
from Bio.Align.Applications import MafftCommandline
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

tempFileText = ''
mafft_align_parse = ''

asteriskHead = ''
alignsForPyperclip = ''

mode = 1

# Create regex and add character so regex can function
stdoutRegex = re.compile(r'(>[^>].+?)\n([RHKDESTNQCUGPAILMFWYVX\n-]+?)>', re.DOTALL)


# Definition for string chunker function
def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

# Default file location. Will need to be changed by the user if they wish to use a different folder to house their muscle .exe
mafft_exe = r'"D:\Bioinformatics\Protein_analysis\mafft-7.311-win64-signed\mafft-win"'

# Define settings to use
print('MAFFT SETTINGS SPECIFICATION')
print('')
print('Do you want to use E-INSi or L-INSi? Type \'E\' or \'L\'')
while True:
        insi_setting = input()
        if insi_setting.lower() != 'e' and insi_setting.lower() != 'l':
                print('You mistyped something. Enter \'E\' for E-INSi or \'L\' for L-INSi then press ENTER')
        else:
                break
print('')

# Prompt user to call for sequences
print('Do you want to refer to a fasta file, or something in your clipboard? Enter 1 for fasta, or 2 for clipboard')
print('If you want to change this option after this, you\'ll need to re-run this script')
while True:
        clipOrFasta = input()
        if clipOrFasta != '1' and clipOrFasta != '2':
                print('You mistyped something. Enter 1 for fasta or 2 for clipboard then press ENTER')
        else:
                break
print('')

# Depending on above, retrieve sequences
while True:
    ##### FASTA FILE #####
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
            # Create a temporary file to be aligned by MAFFT
            shutil.copy(fileName + '.fasta', mafft_exe[1:-1])
            os.rename(mafft_exe[1:-1] + '\\' + fileName + '.fasta', mafft_exe[1:-1] + '\\temp_for_mafft_align.fasta')
            
            # Feed this data into mafft through the CMD window
            if insi_setting.lower() == 'e':
                os.system(r'cd ' + mafft_exe + ' & mafft.bat --genafpair temp_for_mafft_align.fasta > temp_mafft_align_out.fasta')
            else:
                os.system(r'cd ' + mafft_exe + ' & mafft.bat --localpair temp_for_mafft_align.fasta > temp_mafft_align_out.fasta')

            # Read MAFFT output file
            mafft_align_file = open(mafft_exe[1:-1] + '\\temp_mafft_align_out.fasta', 'rU')     # Reads in the output file
            records = list(SeqIO.parse(mafft_align_file, 'fasta'))                              # Parse the sequences
            mafft_align_file.close()

            # Create new MAFFT output file
            print('What do you want to call the alignment output file?')
            outName = input()
            
            for record in records:
                mafft_align_parse += '>' + record.id + '\n' + str(record.seq) + '\n'
                transIDs.append(record.id)                                                      # These two appends are for the asterisk header ouput formatting later
                aligns.append(str(record.seq))
            mafft_align_parse = mafft_align_parse[:-1]                                          # Remove the last newline
            mafftOut = open(outName + '.fasta', 'w')
            mafftOut.write(mafft_align_parse)
            mafftOut.close()

            # Clean up (i.e. remove temporary file and move output file)
            os.remove(mafft_exe[1:-1] + r'\temp_for_mafft_align.fasta')
            os.remove(mafft_exe[1:-1] + r'\temp_mafft_align_out.fasta')

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

            # Prompt user to close this file
            print('When you\'re done with looking at this window, either click out of it (X at the top right) or press ENTER in this window')
            buttonPress=input()
            quit()


    ##### CLIPBOARD #####            
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
                        if tempSplit[i] == '-' or tempSplit[i+1] == '-':
                            continue
                        records.append(SeqRecord(Seq(tempSplit[i+1]), id= '>' + tempSplit[i], description = ""))           # Creates a SeqRecord object for each sequence with just the sequence and its transcript ID
                else:
                    for i in range(0, len(tempSplit)):
                        if tempSplit[i] == '-':
                            continue
                        records.append(SeqRecord(Seq(tempSplit[i]), id= '>sequence' + str(i+1), description = ""))           # Creates a SeqRecord object for each sequence with just the sequence and an invented transcript ID
            # Create a temporary file to be aligned by MAFFT
            for entry in records:
                tempFileText += entry.id + '\n' + str(entry.seq) + '\n'
            file = open(mafft_exe[1:-1] + '\\temp_for_mafft_align.txt', 'w')        # Strip the "" from the start and end of mafft_exe
            file.write(tempFileText)
            file.close()

            # Feed this data into mafft through the CMD window
            if insi_setting.lower() == 'e':
                mafft_cmd = os.path.join(mafft_exe, 'mafft.bat') + ' --genafpair ' + os.path.join(mafft_exe, 'temp_for_mafft_align.txt') + ' > ' + os.path.join(mafft_exe, 'temp_mafft_align_out.fasta')
                run_mafft = subprocess.Popen(mafft_cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
                mafftout, maffterr = run_mafft.communicate()
            else:
                mafft_cmd = os.path.join(mafft_exe, 'mafft.bat') + ' --localpair ' + os.path.join(mafft_exe, 'temp_for_mafft_align.txt') + ' > ' + os.path.join(mafft_exe, 'temp_mafft_align_out.fasta')
                run_mafft = subprocess.Popen(mafft_cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, shell = True)
                mafftout, maffterr = run_mafft.communicate()
                
            # Read MAFFT output file
            mafft_align_file = open(mafft_exe[1:-1] + '\\temp_mafft_align_out.fasta', 'rU')     # Reads in the output file
            records = list(SeqIO.parse(mafft_align_file, 'fasta'))                              # Parse the sequences
            mafft_align_file.close()                                                            

            # Create new MAFFT output file
            print('What do you want to call the alignment output file?')
            outName = input()
            
            for record in records:
                mafft_align_parse += '>' + record.id + '\n' + str(record.seq) + '\n'
                transIDs.append(record.id)                                                      # These two appends are for the asterisk header ouput formatting later
                aligns.append(str(record.seq))
            mafft_align_parse = mafft_align_parse[:-1]                                          # Remove the last newline
            mafftOut = open(outName + '.fasta', 'w')
            mafftOut.write(mafft_align_parse)
            mafftOut.close()

            # Clean up (i.e. remove temporary file and move output file)
            os.remove(mafft_exe[1:-1] + r'\temp_for_mafft_align.txt')
            os.remove(mafft_exe[1:-1] + r'\temp_mafft_align_out.fasta')

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
        for i in range(len(aligns)):                    
            alignsForPyperclip += transIDs[i] + '\t' + aligns[i] + '\n'            # Get all the alignments into something that can be called in a single function. We +1 to the transIDs call since we added a blank line previously
        pyperclip.copy(alignsForPyperclip)

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
        tempFileText = ''
        mafft_align_parse = ''
        
        print('')
