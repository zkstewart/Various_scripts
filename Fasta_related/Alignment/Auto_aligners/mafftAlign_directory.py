import os, re
from Bio import SeqIO

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

# Depending on above, retrieve sequences
while True:
    ##### FASTA FILE #####
        # Get file names
        filenames = []
        for file in os.listdir('.'):
            if not file.endswith('.py'):
                filenames.append(file)
        for fileName in filenames:            
            # Feed this data into mafft through the CMD window
            outName = fileName.rsplit('.',maxsplit=1)[0] + '_align.' + fileName.rsplit('.',maxsplit=1)[1]
            if insi_setting.lower() == 'e':
                os.system(r'cd ' + mafft_exe + ' & mafft.bat --genafpair ' + os.path.join(os.getcwd(), fileName) + ' > ' + outName)
            else:
                os.system(r'cd ' + mafft_exe + ' & mafft.bat --localpair ' + os.path.join(os.getcwd(), fileName) + ' > ' + outName)
