#! python3

# micro_table_analysis
# Current goals - 1) parse table file and generate dictionaries containing
# all relevant details. 2) convert dictionaries into table files amenable
# to R statistical analysis. 3) automatically perform R analysis as part of
# this program.

# Load packages
import os, argparse

### Various functions to perform operations throughout the program
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.tableFile):
                print('I am unable to locate the input table file (' + args.tableFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def micro_table_chunks(tableFileName):
        # Split table file into chunks for parsing - each chunk should be separated by a completely empty line
        chunks = []
        with open(tableFileName, 'r') as fileIn:
                currChunk = []
                for line in fileIn:
                        if line.startswith('Plate ID'):
                                if currChunk == []:
                                        # This is the first chunk, just start it up
                                        currChunk.append(line)
                                else:
                                        # Hold onto current chunk
                                        chunks.append(currChunk)
                                        currChunk = [line]
                        elif set(line) == {'\t', '\r', '\n'} or set(line) == {'\t', '\n'}:   # i.e., if this is a blank line
                                continue
                        elif currChunk != []:
                                currChunk.append(line)
        chunks.append(currChunk)        # Add the last chunk
        return chunks

def drop_hidden_chunks(chunks):
        for i in range(len(chunks)-1,-1,-1):
                line3 = chunks[i][2].split('\t')
                if line3[13] == '' or '?' in line3[13]:
                        del chunks[i]
                # Hard-coded extra handling
                elif line3[0] == 'YM181 phone Yvette Xiao':
                        del chunks[i]
        return chunks

def micro_chunk_to_dict(chunks):
        # Function setup
        chunkDict = {}
        # Loop
        for chunk in chunks:
                # Plate ID and code
                plateID, plateCode = chunk[2].split('\t')[0], chunk[2].split('\t')[1]
                assert plateID not in chunkDict
                chunkDict[plateID] = {'Code': plateCode}
                # Colony morphology and counts - before / after indices
                line1 = chunk[0].split('\t')
                for i in range(len(line1)):
                        if 'before' in line1[i].lower():
                                beforeIndex = i
                        elif 'after' in line1[i].lower():
                                afterIndex = i
                # Colony morphology and counts - group letters and indices
                line2 = chunk[1].split('\t')
                beforeLetters = {}
                afterLetters = {}
                for i in range(len(line2)):
                        if len(line2[i]) == 1 and line2[i] != '\n':     # i.e., if we're looking at a single identifying letter
                                if i < afterIndex - 1:                  # -1 since 'AFTER' is present 1 cell after the total plate count for after
                                        #beforeCounts[line2[i]] = i
                                        beforeLetters[i] = line2[i]
                                else:
                                        #afterCounts[line2[i]] = i
                                        afterLetters[i] = line2[i]
                # Colony morphology and counts - dict start and counts
                line3 = chunk[2].split('\t')
                beforeValues = {}
                afterValues = {}
                for i in range(len(line3)):
                        if i in beforeLetters:
                                currLetter = beforeLetters[i]
                                beforeValues[currLetter] = {'Count': line3[i]}
                        elif i in afterLetters:
                                currLetter = afterLetters[i]
                                afterValues[currLetter] = {'Count': line3[i]}
                        elif i == 2:
                                beforeValues['Total count'] = line3[i]
                        elif i == 8:
                                afterValues['Total count'] = line3[i]
                # Colony morphology and counts - colony morphology
                for x in range(3, len(chunk)):
                        lineX = chunk[x].split('\t')
                        if lineX[0].lower() == 'media':
                                break
                        for i in range(len(lineX)):
                                if i in beforeLetters:
                                        currLetter = beforeLetters[i]
                                        beforeValues[currLetter][lineX[2]] = lineX[i]
                                elif i in afterLetters:
                                        currLetter = afterLetters[i]
                                        afterValues[currLetter][lineX[2]] = lineX[i]
                # Add colony morphology details to chunkDict
                chunkDict[plateID]['Before'] = beforeValues
                chunkDict[plateID]['After'] = afterValues
                # Colony survey data
                line3 = chunk[2].rstrip('\t\r\n').split('\t')
                chunkDict[plateID]['Survey data'] = line3[13:]
                # Media reactions
                chunkDict[plateID]['Media'] = {}
                for x in range(len(chunk)):
                        lineX = chunk[x].split('\t')
                        if lineX[0].lower() != 'media':
                                continue
                        else:
                                for i in range(3, len(lineX)):
                                        if lineX[i] != '' and lineX[i] != '\n' and lineX[i] != '\r\n':
                                                # Extract below results
                                                nbLine = chunk[x+1].split('\t')[i]
                                                hbaLine = chunk[x+2].split('\t')[i]
                                                msaLine = chunk[x+3].split('\t')[i]
                                                mcaLine = chunk[x+4].split('\t')[i]
                                                chunkDict[plateID]['Media'][lineX[i]] = [nbLine, hbaLine, msaLine, mcaLine]
        return chunkDict

def middle_size(inputSize):
        # Set up
        from statistics import mean
        # Strip extra characters
        inputSize = inputSize.strip('<>')
        # Split range
        lower, upper = inputSize.split('-')
        lower, upper = float(lower), float(upper)
        # Find middle value
        middle = mean([lower, upper])
        return str(middle)

def fix_details(chunkDict):
        # Set up
        import copy
        for key, value in chunkDict.items():
                ## HARD CODED - CONFIRM
                if key == 'IX963':
                        continue
                ## Fix details
                # Media - drop everything except A1, B1, etc.
                mediaKeys = list(value['Media'].keys())
                for mKey in mediaKeys:
                        if not mKey.endswith('1'):
                                del value['Media'][mKey]
                # Media - fill in blanks
                mediaKeys = list(value['Media'].keys())
                for mKey in mediaKeys:
                        for i in range(len(value['Media'][mKey])):
                                if value['Media'][mKey][i] == '':
                                        value['Media'][mKey][i] = 'NG'
                # Media - normalise labeling
                mediaKeys = list(value['Media'].keys())
                for mKey in mediaKeys:
                        for i in range(len(value['Media'][mKey])):
                                if value['Media'][mKey][i].lower() == 'no growth':
                                        value['Media'][mKey][i] = 'NG'
                                elif value['Media'][mKey][i].lower() == '-':    # Confirmed
                                        value['Media'][mKey][i] = 'NG'
                                elif value['Media'][mKey][i].lower() == 'g (purple)':
                                        value['Media'][mKey][i] = 'G-LF'
                                elif value['Media'][mKey][i].lower() == 'pallet' or value['Media'][mKey][i].lower() == 'pallets':
                                        value['Media'][mKey][i] = 'pellet'     # Confirmed
                                elif value['Media'][mKey][i].lower() == 'mf & nmf':
                                        value['Media'][mKey][i] = 'MF'          ## CONFIRM IF THIS IS CORRECT
                                # Handle specific typos encountered
                                elif value['Media'][mKey][i].lower() == 'rturbid':
                                        value['Media'][mKey][i] = 'turbid'
                # Before/after - drop empty values
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        for bKey in baKeys:
                                if bKey == 'Total count':
                                        continue
                                baVal = value[group][bKey]
                                allEmpty = True
                                for k, v in baVal.items():
                                        if v != '':
                                                allEmpty = False
                                                break
                                if allEmpty == True:
                                        del value[group][bKey]
                # Before/after - clean up size label and numbers
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        for bKey in baKeys:
                                if bKey == 'Total count':
                                        continue
                                baVal = copy.deepcopy(value[group][bKey])       # Making a deepcopy lets us delete values in the original
                                for k, v in baVal.items():
                                        if 'size' in k.lower():
                                                value[group][bKey]['Size phi (mm)'] = value[group][bKey].pop(k)
                                                # Handle size ranges
                                                if '-' in v:
                                                        ## HARD CODED HANDLE ##
                                                        if v == "2'0.5-3":
                                                                v = '0.5-3'     # Confirmed
                                                        middle = middle_size(v)
                                                        value[group][bKey]['Size phi (mm)'] = middle
                # Before/after - replace blanks with 'ND'
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        for bKey in baKeys:
                                if bKey == 'Total count':
                                        continue
                                for k, v in value[group][bKey].items():
                                        if v == '':
                                                value[group][bKey][k] = 'ND'
                # Before/after - clean up count numbers
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        for bKey in baKeys:
                                if bKey == 'Total count':
                                        continue
                                #baVal = copy.deepcopy(value[group][bKey])       # Making a deepcopy lets us delete values in the original
                                #for k, v in baVal.items():
                                        #if k == 'Count':
                                value[group][bKey]['Count'] = value[group][bKey]['Count'].strip('~ ')
                                if 'tntc' in value[group][bKey]['Count'].lower():
                                        value[group][bKey]['Count'] = 'TNTC'    ## CONFIRM IF THIS IS CORRECT - HANDLES CZ124
                # Before/after - handle missing count numbers
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        for bKey in baKeys:
                                if bKey == 'Total count':
                                        continue
                                if value[group][bKey]['Count'] == '' or value[group][bKey]['Count'] == 'ND':    # It will == 'ND' here since the replace blanks function comes first currently
                                        if value[group]['Total count'] == 'TNTC':
                                                value[group][bKey]['Count'] = 'TNTC'    ## CONFIRM IF THIS IS CORRECT
                                        else:
                                                print('Don\'t know how to handle this missing count cell.')
                                                print(key)
                                                print(value)
                                                stophere
                # Survey data - replace ? with 'Other' value
                questionOthers = ['', '6', '4', '6', '6', '6', '5', '5']        # These numbers refer to the position of the 'Other' option
                for i in range(len(value['Survey data'])):
                        if value['Survey data'][i] == '?':
                                value['Survey data'][i] = questionOthers[i]     # Variable name doubles as a deep insight into humanity
                ## Validate details
                # Media - check labeling for consistency
                mediaKeys = list(value['Media'].keys())
                for mKey in mediaKeys:
                        for i in range(len(value['Media'][mKey])):
                                if i == 0:      # NB 1ml turbid/clear
                                        options = ['turbid', 'clear', 'NG', 'N/A', 'pellet']
                                        if value['Media'][mKey][i] not in options:
                                                print('Problem with values - NB 1ml')
                                                stophere
                                elif i == 1:      # HBA
                                        options = ['alpha', 'beta', 'gamma', 'spreader', 'NG']  # Spreader is a valid group according to Elise email
                                        if value['Media'][mKey][i] not in options:
                                                print('Problem with values - HBA')
                                                stophere
                                elif i == 2:      # MSA
                                        options = ['MF', 'NMF', 'NG']
                                        if value['Media'][mKey][i] not in options:
                                                print('Problem with values - MSA')
                                                stophere
                                elif i == 3:      # MCA
                                        options = ['NG', 'G-LF', 'G']
                                        if value['Media'][mKey][i] not in options:
                                                print('Problem with values - MCA')
                                                stophere
                # Before/after - make sure everything is a number
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        for bKey in baKeys:
                                if bKey == 'Total count':
                                        continue
                                if value[group][bKey]['Count'] == 'TNTC':
                                        continue
                                try:
                                        int(value[group][bKey]['Count'])
                                except:
                                        print('Can\'t turn this into an integer.')
                                        print(key)
                                        print(value)
                                        stophere
        return chunkDict

### USER INPUT
usage = """%(prog)s
"""

# Reqs
p = argparse.ArgumentParser(description=usage)

p.add_argument("-t", "-tableFile", dest="tableFile",
               help="Input microbiology table file name (should be saved as a tab-delimited text file).")
p.add_argument("-o", "--output", dest="outputFileName",
               help="Output results file name")
args = p.parse_args()

## Hard coded for testing
args.tableFile = r'E:\elise\RAW_LQB301_ProjectData25-06-2018(9673).txt'
#args.tableFile = r'E:\elise\DELETED_LQB301_ProjectData25-06-2018(9673).txt'
args.outputFileName = r'E:\elise\test_micro_parse.txt'
validate_args(args)

# Parse table file into chunks
chunks = micro_table_chunks(args.tableFile)

# Drop chunks that were hidden
chunks = drop_hidden_chunks(chunks)

# Turn chunks into dictionaries
chunkDict = micro_chunk_to_dict(chunks)

# Fix details for the dictionary
chunkDict = fix_details(chunkDict)

#### SCRIPT ALL DONE
print('Program completed successfully!')
