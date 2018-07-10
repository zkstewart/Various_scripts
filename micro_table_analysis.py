#! python3

# micro_table_analysis
# Current goals - 1) parse table file and generate dictionaries containing
# all relevant details. 2) convert dictionaries into table files amenable
# to R statistical analysis. 3) automatically perform R analysis as part of
# this program.

# Load packages
import os, argparse

# Various functions to perform operations throughout the program
## Validate arguments
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

## Parse and process raw data
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
                                        beforeLetters[i] = line2[i]
                                else:
                                        afterLetters[i] = line2[i]
                        elif '/' in line2[i]:   # i.e., if we're handling a shared group - only found on plate CZ124
                                letters = line2[i].split('/')
                                for letter in letters:
                                        if i < afterIndex - 1:                  # -1 since 'AFTER' is present 1 cell after the total plate count for after
                                                beforeLetters[i] = letters
                                        else:
                                                afterLetters[i] = letters
                # Colony morphology and counts - dict start and counts
                line3 = chunk[2].split('\t')
                beforeValues = {}
                afterValues = {}
                for i in range(len(line3)):
                        if i in beforeLetters:
                                currLetter = beforeLetters[i]
                                for letter in currLetter:                       # This and similar lines of code are specifically to handle CZ124, but future cases might also crop up
                                        beforeValues[letter] = {'Count': line3[i]}
                        elif i in afterLetters:
                                currLetter = afterLetters[i]
                                for letter in currLetter:
                                        afterValues[letter] = {'Count': line3[i]}
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
                                        for letter in currLetter:
                                                beforeValues[letter][lineX[2]] = lineX[i]
                                elif i in afterLetters:
                                        currLetter = afterLetters[i]
                                        for letter in currLetter:
                                                afterValues[letter][lineX[2]] = lineX[i]
                # Add colony morphology details to chunkDict
                chunkDict[plateID]['Before'] = beforeValues
                chunkDict[plateID]['After'] = afterValues
                # Colony survey data
                line3 = chunk[2].rstrip('\t\r\n').split('\t')
                assert len(line3[13:]) == 8     # Make sure the data was entered correctly
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
        if '-' in inputSize:
                lower, upper = inputSize.split('-')
        else:
                lower, upper = inputSize.split('&')
                lower, upper = lower.strip(' '), upper.strip(' ')
        lower, upper = float(lower), float(upper)
        # Find middle value
        middle = mean([lower, upper])
        return str(middle)

def fix_details(chunkDict):
        # Set up
        import copy
        for key, value in chunkDict.items():
                ## HARD CODED - CONFIRM
                #if key == 'CZ124' or key == 'CA152' or key == 'FJ000' or key == 'K123M':
                #        continue
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
                                                if '-' in v or '&' in v:
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
                                        if 'tntc' in value[group][bKey].lower():
                                                value[group][bKey] = 'TNTC'
                                        continue
                                value[group][bKey]['Count'] = value[group][bKey]['Count'].strip('~ ')
                                if 'tntc' in value[group][bKey]['Count'].lower():
                                        value[group][bKey]['Count'] = 'TNTC'    # Not strictly confirmed, but assumed correct - HANDLES CZ124
                # Before/after - add in missing numbers when a total is present
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        if 'Total count' in baKeys:
                                del baKeys[baKeys.index('Total count')]
                        if value[group]['Total count'].isdigit():
                                totalSum = 0
                                numGroups = 0
                                for i in range(len(baKeys)):
                                        if value[group][baKeys[i]]['Count'] == '' or value[group][baKeys[i]]['Count'] == 'ND' or value[group][baKeys[i]]['Count'] == 'TNTC':
                                                numGroups += 1
                                        else:
                                                totalSum += int(value[group][baKeys[i]]['Count'])
                                if numGroups == 0:
                                        continue
                                # Divide the total count across each morphology as evenly as possible
                                totalNum = int(value[group]['Total count']) - totalSum
                                valueList = divide_num_to_list(totalNum, numGroups)
                                # Update values
                                ongoingCount = 0
                                for i in range(len(baKeys)):
                                        if value[group][baKeys[i]]['Count'] == '' or value[group][baKeys[i]]['Count'] == 'ND' or value[group][baKeys[i]]['Count'] == 'TNTC':
                                                value[group][baKeys[i]]['Count'] = str(valueList[ongoingCount])
                                                ongoingCount += 1
                # Before/after - handle missing count numbers when total is TNTC
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        if 'Total count' in baKeys:
                                del baKeys[baKeys.index('Total count')]
                        for i in range(len(baKeys)):
                                if value[group][baKeys[i]]['Count'] == '' or value[group][baKeys[i]]['Count'] == 'ND':    # It will == 'ND' here since the replace blanks function comes first currently
                                        if value[group]['Total count'] == 'TNTC':
                                                value[group][baKeys[i]]['Count'] = 'TNTC'    # Confirmed (for BM998)
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
                # Survey data - replace absent values with 0
                for i in range(len(value['Survey data'])):
                        if value['Survey data'][i] == '':
                                value['Survey data'][i] = '0'
                # Survey data - replace incorrect values with 'Other' value
                for i in range(len(value['Survey data'])):
                        if i == 0:
                                continue
                        if int(value['Survey data'][i]) > int(questionOthers[i]):
                                value['Survey data'][i] = questionOthers[i]     # This handles MAC97, CONFIRM with Elise?
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
                                        if value[group][bKey] == 'TNTC':
                                                continue
                                        else:
                                                try:
                                                        int(value[group][bKey])
                                                except:
                                                        print('Can\'t turn this into an integer.')
                                                        print(key)
                                                        print(value)
                                                        stophere
                                else:
                                        if value[group][bKey]['Count'] == 'TNTC':
                                                continue
                                        else:
                                                try:
                                                        int(value[group][bKey]['Count'])
                                                except:
                                                        print('Can\'t turn this into an integer.')
                                                        print(key)
                                                        print(value)
                                                        stophere
                # Before/after - make sure individual numbers sum to total count
                for group in ['Before', 'After']:
                        baKeys = list(value[group].keys())
                        if 'Total count' in baKeys:
                                del baKeys[baKeys.index('Total count')]
                        totalCount = value[group]['Total count']
                        # Scenario 1: total count is TNTC
                        if totalCount == 'TNTC':
                                found = False
                                for bKey in baKeys:
                                        if value[group][bKey]['Count'] == 'TNTC':
                                                found = True
                                if found == False:
                                        print('Doesn\'t add up to TNTC.')
                                        print(key)
                                        print(value)
                                        stophere
                        # Scenario 2: total count is an integer
                        else:
                                totalCount = int(totalCount)
                                sumCount = 0
                                for bKey in baKeys:
                                        sumCount += int(value[group][bKey]['Count'])
                                if sumCount != totalCount:
                                        print('Doesn\'t add up to total count')
                                        print(key)
                                        print(value)
                                        stophere
        return chunkDict

## Data extraction and reformating
def count_categoriser(inputValue, returnType):
        #categoriesLabel = ['ND', '0', '1-10', '11-50', '50-100', '100-200', '200-400', 'TNTC']
        #categoriesRange = ['ND', 0, range(1,11), range(11,51), range(51,101), range(101,201), range(201,401), 'TNTC']
        categoriesLabel = ['ND', '0', '1-10', '11-100', '100-340', 'TNTC']
        categoriesAsNumeric = ['ND', '1', '2', '3', '4', '5']
        categoriesRange = ['ND', 0, range(1,11), range(11,101), range(101,341), 'TNTC']
        for i in range(len(categoriesRange)):
                if inputValue == categoriesRange[i]:
                        if returnType == 'label':
                                return categoriesLabel[i]
                        else:
                                return categoriesAsNumeric[i]
                elif inputValue.isdigit() and categoriesRange[i] != 'ND' and categoriesRange[i] != 'TNTC':
                        # Handle 0 comparison
                        if int(inputValue) == categoriesRange[i]:
                                if returnType == 'label':
                                        return categoriesLabel[i]
                                else:
                                        return categoriesAsNumeric[i]
                        # Handle other number comparisons
                        if categoriesRange[i] == 0:
                                continue
                        elif int(inputValue) in categoriesRange[i]:
                                if returnType == 'label':
                                        return categoriesLabel[i]
                                else:
                                        return categoriesAsNumeric[i]

def growth_categoriser(inputValue, mediaType):
        nbCatLabels = ['ND', 'N/A', 'NG', 'turbid', 'clear', 'pellet']  # Need to handle plate 9K78A, mould has N/A value
        nbCatAsNumeric = ['ND', '0', '1', '2', '3', '4']
        hbaCatLabels = ['ND', 'NG', 'alpha', 'beta', 'gamma', 'spreader']
        hbaCatAsNumeric = ['ND', '1', '2', '3', '4', '5']
        msaCatLabels = ['ND', 'NG', 'NMF', 'MF']
        msaCatAsNumeric = ['ND', '1', '2', '3']
        mcaCatLabels = ['ND', 'NG', 'G', 'G-LF']
        mcaCatAsNumeric = ['ND', '1', '2', '3']
        if mediaType == 'NB':
                for i in range(len(nbCatLabels)):
                        if inputValue == nbCatLabels[i]:
                                return nbCatAsNumeric[i]
        elif mediaType == 'HBA':
                for i in range(len(hbaCatLabels)):
                        if inputValue == hbaCatLabels[i]:
                                return hbaCatAsNumeric[i]
        elif mediaType == 'MSA':
                for i in range(len(msaCatLabels)):
                        if inputValue == msaCatLabels[i]:
                                return msaCatAsNumeric[i]
        elif mediaType == 'MCA':
                for i in range(len(mcaCatLabels)):
                        if inputValue == mcaCatLabels[i]:
                                return mcaCatAsNumeric[i]
        else:
                print('Something went wrong. Check this out (growth categoriser).')
                stophere

def growth_before_after_extract(dictValue):
        # Figure out which letters are before/after
        beforeKeys = list(dictValue['Before'].keys())
        del beforeKeys[beforeKeys.index('Total count')]
        afterKeys = list(dictValue['After'].keys())
        del afterKeys[afterKeys.index('Total count')]
        # Categorise media values by before/after
        mediaDict = {'Before': {'NB': [], 'HBA': [], 'MSA': [], 'MCA': []}, 'After': {'NB': [], 'HBA': [], 'MSA': [], 'MCA': []}}
        for key, value in dictValue['Media'].items():
                letter = key[0]
                if letter in beforeKeys:
                        mediaDict['Before']['NB'].append(growth_categoriser(value[0], 'NB'))
                        mediaDict['Before']['HBA'].append(growth_categoriser(value[1], 'HBA'))
                        mediaDict['Before']['MSA'].append(growth_categoriser(value[2], 'MSA'))
                        mediaDict['Before']['MCA'].append(growth_categoriser(value[3], 'MCA'))
                elif letter in afterKeys:
                        mediaDict['After']['NB'].append(growth_categoriser(value[0], 'NB'))
                        mediaDict['After']['HBA'].append(growth_categoriser(value[1], 'HBA'))
                        mediaDict['After']['MSA'].append(growth_categoriser(value[2], 'MSA'))
                        mediaDict['After']['MCA'].append(growth_categoriser(value[3], 'MCA'))
                else:
                        print('Something went wrong. Check this out.')
                        print(dictValue)
                        print(key)
        return mediaDict

def total_count_quiz_response(chunkDict):
        # Set up
        questionOthers = ['', 6, 4, 6, 6, 6, 5, 5] # These are the numbers of responses possible for each quiz question
        # Format a list to store values in
        beforeFormatList = []
        afterFormatList = []
        for i in range(len(questionOthers)):
                beforeFormatList.append({})
                afterFormatList.append({})
                if i == 0:
                        mainMedItemDict, varMedItemDict = medical_item_dicts()
                        for value in mainMedItemDict.values():
                                beforeFormatList[i][value] = []
                                afterFormatList[i][value] = []
                else:
                        for x in range(1, questionOthers[i] + 1):
                                beforeFormatList[i][x] = []
                                afterFormatList[i][x] = []
        # Loop through our dictionary of values and extract relevant details
        for key, value in chunkDict.items():
                for i in range(len(value['Survey data'])):
                        # Specifically handle medical items
                        if i == 0:
                                medItem = value['Survey data'][i]
                                medItem = medical_item_categoriser(medItem)
                                # Convert the count to a category
                                countCatBefore = count_categoriser(value['Before']['Total count'], 'number')
                                countCatAfter = count_categoriser(value['After']['Total count'], 'number')
                                # Add it to our list
                                if medItem != False and countCatBefore != 'ND' and countCatAfter != 'ND':
                                        beforeFormatList[i][medItem].append(countCatBefore)        # This means we should have ordered lists
                                        afterFormatList[i][medItem].append(countCatAfter)
                        # Skip any samples which lack a response for this question
                        elif value['Survey data'][i] == '0':
                                continue
                        else:
                                # Convert the count to a category
                                countCatBefore = count_categoriser(value['Before']['Total count'], 'number')
                                countCatAfter = count_categoriser(value['After']['Total count'], 'number')
                                # Add it to our respective list
                                if countCatBefore != 'ND' and countCatAfter != 'ND':      # Skip anything which lacks values
                                        beforeFormatList[i][int(value['Survey data'][i])].append(countCatBefore)        # This means we should have ordered lists
                                        afterFormatList[i][int(value['Survey data'][i])].append(countCatAfter)
        return beforeFormatList, afterFormatList

def microbe_diversity_quiz_response(chunkDict):
        # Set up
        questionOthers = ['', 6, 4, 6, 6, 6, 5, 5] # These are the numbers of responses possible for each quiz question
        # Format a list to store values in
        beforeFormatList = []
        afterFormatList = []
        for i in range(len(questionOthers)):
                beforeFormatList.append({})
                afterFormatList.append({})
                if i == 0:
                        mainMedItemDict, varMedItemDict = medical_item_dicts()
                        for value in mainMedItemDict.values():
                                beforeFormatList[i][value] = []
                                afterFormatList[i][value] = []
                else:
                        for x in range(1, questionOthers[i] + 1):
                                beforeFormatList[i][x] = []
                                afterFormatList[i][x] = []
        # Loop through our dictionary of values and extract relevant details
        for key, value in chunkDict.items():
                for i in range(len(value['Survey data'])):
                        # Specifically handle medical items
                        if i == 0:
                                medItem = value['Survey data'][i]
                                medItem = medical_item_categoriser(medItem)
                                # Compute the diversity count
                                beforeDiv, afterDiv = diversity_number(value)
                                # Add it to our list
                                if medItem != False:
                                        beforeFormatList[i][medItem].append(beforeDiv)        # This means we should have ordered lists
                                        afterFormatList[i][medItem].append(afterDiv)
                        # Skip any samples which lack a response for this question
                        elif value['Survey data'][i] == '0':
                                continue
                        else:
                                # Compute the diversity count
                                beforeDiv, afterDiv = diversity_number(value)
                                # Add it to our respective list
                                beforeFormatList[i][int(value['Survey data'][i])].append(beforeDiv)        # This means we should have ordered lists
                                afterFormatList[i][int(value['Survey data'][i])].append(afterDiv)
        return beforeFormatList, afterFormatList

def media_reaction_quiz_response(chunkDict):
        # Set up
        questionOthers = ['', 6, 4, 6, 6, 6, 5, 5] # These are the numbers of responses possible for each quiz question
        # Format a dict to store values in
        beforeFormatDict = {'NB': [], 'HBA': [], 'MSA': [], 'MCA': []}
        afterFormatDict = {'NB': [], 'HBA': [], 'MSA': [], 'MCA': []}
        for x in beforeFormatDict.keys():
                for i in range(len(questionOthers)):
                        beforeFormatDict[x].append({})
                        afterFormatDict[x].append({})
                        if i == 0:
                                mainMedItemDict, varMedItemDict = medical_item_dicts()
                                for value in mainMedItemDict.values():
                                        beforeFormatDict[x][i][value] = []
                                        afterFormatDict[x][i][value] = []
                        else:
                                for y in range(1, questionOthers[i] + 1):
                                        beforeFormatDict[x][i][y] = []
                                        afterFormatDict[x][i][y] = []
        # Loop through our dictionary of values and extract relevant details
        for key, value in chunkDict.items():
                # Pull out relevant details for before/after
                mediaDict = growth_before_after_extract(value)
                beforeKeys = list(beforeFormatDict.keys())
                for x in range(len(beforeKeys)):
                        for i in range(len(value['Survey data'])):
                                # Extract data from media dict
                                beforeMedia = mediaDict['Before'][beforeKeys[x]]
                                afterMedia = mediaDict['After'][beforeKeys[x]]
                                # Remove ND entries
                                while 'ND' in beforeMedia:
                                        del beforeMedia[beforeMedia.index('ND')]
                                while 'ND' in afterMedia:
                                        del afterMedia[afterMedia.index('ND')]
                                # Skip entries with no media results
                                #if beforeMedia == [] or afterMedia == []:
                                #        continue
                                # Specifically handle medical items
                                if i == 0:
                                        medItem = value['Survey data'][i]
                                        medItem = medical_item_categoriser(medItem)
                                        # Add it to our list
                                        if medItem != False:
                                                beforeFormatDict[beforeKeys[x]][i][medItem] += beforeMedia        # This means we should have ordered lists
                                                afterFormatDict[beforeKeys[x]][i][medItem] += afterMedia
                                # Skip any samples which lack a response for this question
                                elif value['Survey data'][i] == '0':
                                        continue
                                else:
                                        # Compute the diversity count
                                        beforeDiv, afterDiv = diversity_number(value)
                                        # Add it to our respective list
                                        beforeFormatDict[beforeKeys[x]][i][int(value['Survey data'][i])] += beforeMedia        # This means we should have ordered lists
                                        afterFormatDict[beforeKeys[x]][i][int(value['Survey data'][i])] += afterMedia
        return beforeFormatDict, afterFormatDict

def medical_item_dicts():
        # Categorise and normalise naming for most abundant items
        mainMedItemDict = {'STETHOSCOPE': 'STETHOSCOPE',
                       'PEN': 'PEN', 'MARKER PEN': 'PEN', 'BIRO': 'PEN', 'INK PEN': 'PEN', 'WRITING PEN': 'PEN',
                       'SAFETY GLASSES': 'SAFETY GLASSES',
                       'PEN LIGHT': 'PEN LIGHT',
                       'MOBILE PHONE': 'MOBILE PHONE', 'PHONE': 'MOBILE PHONE', 'MOBILE PHONE CASE': 'MOBILE PHONE',
                       'GLASSES': 'GLASSES', 'EYE GLASSES': 'GLASSES', 'PRESCRIPTION GLASSES': 'GLASSES',
                       'WORK SHOES': 'CLOTHES', 'COAT COLLAR': 'CLOTHES', 'LAB COAT': 'CLOTHES', 'WATCH': 'CLOTHES',
                       'SCISSORS': 'SCISSORS', 'TRAUMA SCISSORS': 'SCISSORS', 'TRAUMA SHEERS': 'SCISSORS'}
        # Categorise various items - unlikely we can do meaningful statistics with these items due to low count
        varMedItemDict = {'LAPTOP': 'BACKPACK STORED', 'METAL RULER': 'BACKPACK STORED',
                       'ID CARD': 'ID CARD',
                       'EPIPEN': 'EPIPEN',
                       'BLOOD PRESSURE CUFF': 'BLOOD PRESSURE CUFF',
                       'STETHOSCOPE-VET': 'STETHOSCOPE-VET',
                       'MICROPIPETTOR': 'MICROPIPETTOR',        # Drop this item - it isn't the students?
                       'UNSPECIFIED': 'UNSPECIFIED'}
        return mainMedItemDict, varMedItemDict

def medical_item_categoriser(medItem):
        mainMedItemDict, varMedItemDict = medical_item_dicts()
        if medItem in mainMedItemDict:
                return mainMedItemDict[medItem]
        elif medItem in varMedItemDict:
                return False
        else:
                print('Failed to find this item in the dictionary!')
                print(medItem)
                print('Need to fix the hard coding of this script to address this problem. Ask Zac if he\'s not the one using this and he\'s still around.')
                quit()

def diversity_number(chunkValue):
        beforeDiv = len(chunkValue['Before'].keys()) - 1        # -1 for 'Total count' value
        afterDiv = len(chunkValue['After'].keys()) - 1
        return beforeDiv, afterDiv


### Overview of data for tabulation
def tabulate_data(chunkDict):
        # Set up
        responseOccurrence = {1: {'OTHER': 0}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}, 7: {}, 8: {}}
        countOccurBefore = {'ND': 0, '0': 0, '1-10': 0, '11-100': 0, '100-340': 0, 'TNTC': 0}
        countOccurAfter = {'ND': 0, '0': 0, '1-10': 0, '11-100': 0, '100-340': 0, 'TNTC': 0}
        medItemCountBefore = {'OTHER': {'0': 0, '1-10': 0, '100-340': 0, '11-100': 0, 'TNTC': 0}}
        medItemCountAfter = {'OTHER': {'0': 0, '1-10': 0, '100-340': 0, '11-100': 0, 'TNTC': 0}}
        for key, value in chunkDict.items():
                # Extract relevant details for tabulation
                ## Quiz responses
                survey = value['Survey data']
                #### R1
                medicalItem = medical_item_categoriser(survey[0])
                if medicalItem == False:
                        responseOccurrence[1]['OTHER'] += 1
                elif medicalItem not in responseOccurrence[1]:
                        responseOccurrence[1][medicalItem] = 1
                else:
                        responseOccurrence[1][medicalItem] += 1
                #### R2
                if survey[1] not in responseOccurrence[2]:
                        responseOccurrence[2][survey[1]] = 1
                else:
                        responseOccurrence[2][survey[1]] += 1
                #### R3
                if survey[2] not in responseOccurrence[3]:
                        responseOccurrence[3][survey[2]] = 1
                else:
                        responseOccurrence[3][survey[2]] += 1
                #### R4
                if survey[3] not in responseOccurrence[4]:
                        responseOccurrence[4][survey[3]] = 1
                else:
                        responseOccurrence[4][survey[3]] += 1
                #### R5
                if survey[4] not in responseOccurrence[5]:
                        responseOccurrence[5][survey[4]] = 1
                else:
                        responseOccurrence[5][survey[4]] += 1
                #### R6
                if survey[5] not in responseOccurrence[6]:
                        responseOccurrence[6][survey[5]] = 1
                else:
                        responseOccurrence[6][survey[5]] += 1
                #### R7
                if survey[6] not in responseOccurrence[7]:
                        responseOccurrence[7][survey[6]] = 1
                else:
                        responseOccurrence[7][survey[6]] += 1
                #### R8
                if survey[7] not in responseOccurrence[8]:
                        responseOccurrence[8][survey[7]] = 1
                else:
                        responseOccurrence[8][survey[7]] += 1
                # Microbial counts
                countCatBefore = count_categoriser(value['Before']['Total count'], 'label')
                countOccurBefore[countCatBefore] += 1
                countCatAfter = count_categoriser(value['After']['Total count'], 'label')
                countOccurAfter[countCatAfter] += 1
                # Other format
                if medicalItem == False:
                        if countCatBefore not in medItemCountBefore['OTHER']:
                                medItemCountBefore['OTHER'][countCatBefore] = 1
                        else:
                                medItemCountBefore['OTHER'][countCatBefore] += 1
                        if countCatAfter not in medItemCountAfter['OTHER']:
                                medItemCountAfter['OTHER'][countCatAfter] = 1
                        else:
                                medItemCountAfter['OTHER'][countCatAfter] += 1
                else:
                        # Establish key-value pairs for medical items
                        if medicalItem not in medItemCountBefore:
                                medItemCountBefore[medicalItem] = {'0': 0, '1-10': 0, '100-340': 0, '11-100': 0, 'TNTC': 0}     # used to be {}
                        if medicalItem not in medItemCountAfter:
                                medItemCountAfter[medicalItem] = {'0': 0, '1-10': 0, '100-340': 0, '11-100': 0, 'TNTC': 0}
                        # Add counts to medical items
                        #if countCatBefore not in medItemCountBefore[medicalItem]:
                        #        medItemCountBefore[medicalItem][countCatBefore] = 1
                        #else:
                        medItemCountBefore[medicalItem][countCatBefore] += 1
                        #if countCatAfter not in medItemCountAfter[medicalItem]:
                        #        medItemCountAfter[medicalItem][countCatAfter] = 1
                        #else:
                        medItemCountAfter[medicalItem][countCatAfter] += 1
        # Format data
        ## Quiz response
        ### R1
        r1Table = dict_tabulate(responseOccurrence[1], None)
        ### R2
        r2four5sixText = {'0': 'No response', '1': 'Never', '2': 'Once a year', '3': 'Once a week', '4': 'Daily', '5': 'After each patient', '6': 'Other'}
        r2Table = dict_tabulate(responseOccurrence[2], r2four5sixText)
        ### R3
        r3Text = {'0': 'No response', '1': 'Water only', '2': 'Soap and water', '3': '70% ethanol/isopropanol', '4': 'Other'}
        r3Table = dict_tabulate(responseOccurrence[3], r3Text)
        ### R4
        r4Table = dict_tabulate(responseOccurrence[4], r2four5sixText)
        ### R5
        r5Table = dict_tabulate(responseOccurrence[5], r2four5sixText)
        ### R6
        r6Table = dict_tabulate(responseOccurrence[6], r2four5sixText)
        ### R7
        r7eightText = {'0': 'No response', '1': 'Strongly disagree', '2': 'Disagree', '3': 'Neutral', '4': 'Agree', '5': 'Strongly agree'}
        r7Table = dict_tabulate(responseOccurrence[7], r7eightText)
        ### R8
        r8Table = dict_tabulate(responseOccurrence[8], r7eightText)
        ### Combine
        quizROut = ['Variable\tn=\t%', 'Medical equipment item\t\t']
        quizROut += r1Table
        quizROut += ['Frequency of medical item disinfection\t\t']
        quizROut += r2Table
        quizROut += ['Agents used for cleaning medical equipment\t\t']
        quizROut += r3Table
        quizROut += ['Frequency of medical equipment item disinfection of components in direct contact with patients\t\t']
        quizROut += r4Table
        quizROut += ['Frequency of medical equipment item disinfection of components in indirect contact with patients\t\t']
        quizROut += r5Table
        quizROut += ['Frequency of medical equipment item disinfection of components in direct contact with myself\t\t']
        quizROut += r6Table
        quizROut += ['Infection control practice is critical to protecting me from infectious diseases\t\t']
        quizROut += r7Table
        quizROut += ['Infection control practice is critical to protecting my patients/clients from infectious diseases\t\t']
        quizROut += r8Table
        ## Microbial counts
        microOut = micro_count_tabulate(medItemCountBefore, medItemCountAfter)
        # Join results tables
        combinedOut = ['Table 1. Infection control perceptions and practices of student HCWs'] + quizROut + ['Table 2. Microbial enumeration before and after prescribed medical equipment disinfection'] + microOut
        combinedOut = '\n'.join(combinedOut)
        # Return results
        return combinedOut

def micro_count_tabulate(medItemCountBefore, medItemCountAfter):
        outList = ['Medical equipment item\tNo growth\t1-10\t11-100\t100-340\tTNTC']
        order = ['CLOTHES', 'MOBILE PHONE', 'PEN', 'STETHOSCOPE', 'SCISSORS', 'SAFETY GLASSES', 'GLASSES', 'PEN LIGHT', 'OTHER']
        for key in order:
                beforeRow = dict_tabulate(medItemCountBefore[key], None)
                afterRow = dict_tabulate(medItemCountAfter[key], None)
                # Reformat rows a bit
                outRow = [['Before'], ['After']]
                for entry in beforeRow:
                        splitEntry = entry.split('\t')
                        outRow[0].append(splitEntry[1] + ' (' + splitEntry[2] + ')')
                for entry in afterRow:
                        splitEntry = entry.split('\t')
                        outRow[1].append(splitEntry[1] + ' (' + splitEntry[2] + ')')
                outList.append(key)
                # Format the outrow a bit more
                for i in range(len(outRow)):
                        outRow[i] = '\t'.join(outRow[i])
                outList += outRow
        return outList

def dict_tabulate(inputDict, replaceDict):
        pairList = []
        # Preliminary format
        for key, value in inputDict.items():
                pairList.append([key, str(value)])
        pairList.sort()
        # Replace keys with their proper values
        if replaceDict != None:
                pairList = pair_replace(pairList, replaceDict)
        else:
                # Handle dict_tabulate calls for quiz responses
                otherEntry = None
                for entry in pairList:
                        if entry[0] == 'OTHER':
                                otherEntry = entry
                pairList.sort(key = lambda x: -int(x[1]))
                if otherEntry != None:
                        del pairList[pairList.index(otherEntry)]
                        pairList.append(otherEntry)
                # Handle dict_tabulate calls for medical item counts
                else:
                        order = ['0', '1-10', '11-100', '100-340', 'TNTC']
                        newPairList = []
                        for entry in order:
                                for pair in pairList:
                                        if pair[0] == entry:
                                                newPairList.append(pair)
                        pairList = newPairList
        # Tally count
        total = 0
        for i in range(len(pairList)):
                total += int(pairList[i][1])
        # Reformat list and compute percent proportions
        for i in range(len(pairList)):
                proportion = round(round(int(pairList[i][1]) / total, 3) * 100, 2)      # Python is REALLY weird here - has to be some sort of a bug...
                pairList[i] = '\t'.join([pairList[i][0], str(pairList[i][1]), str(proportion) + '%'])
        return pairList

#def key_replace(inputDict, replaceDict):
#        newDict = {}
#        for key, value in inputDict.items():
#                newKey = replaceDict[key]
#                newDict[newKey] = value
#        return newDict

def pair_replace(pairList, replaceDict):
        for i in range(len(pairList)):
                pairList[i][0] = replaceDict[pairList[i][0]]
        return pairList

### Test 1 & 4
def quiz_cats_to_csv(quizCatDict, prefix, suffix, header, resultDir):
        import os
        fileNames = []
        for i in range(len(quizCatDict)):
                #if i == 0:
                #        continue
                # Get the output directory
                outDirPath = os.path.join(os.getcwd(), resultDir)
                if not os.path.isdir(outDirPath):
                        os.mkdir(outDirPath)
                # Generate an output name
                name = file_name_gen(os.path.join(outDirPath, prefix + suffix), '_' + str(i+1) + '.csv')
                fileNames.append(name)
                # Produce output file
                with open(name, 'w') as fileOut:
                        fileOut.write(header)
                        for key, value in quizCatDict[i].items():
                                for val in value:
                                        fileOut.write(str(key) + ',' + str(val) + '\n')
### Test 2
def compare_total_to_csv(chunkDict,  resultDir):
        import os
        # Get the output directory
        outDirPath = os.path.join(os.getcwd(), resultDir)
        if not os.path.isdir(outDirPath):
                os.mkdir(outDirPath)
        # Generate an output name
        name = file_name_gen(os.path.join(outDirPath, 'quizcat_totals_overall'), '.csv')
        with open(name, 'w') as fileOut:
                fileOut.write('count_category_before,count_category_after\n')
                for key, value in chunkDict.items():
                        # Convert the count to a category
                        countCatBefore = count_categoriser(value['Before']['Total count'], 'number')
                        countCatAfter = count_categoriser(value['After']['Total count'], 'number')
                        # Add it to our respective list
                        if countCatBefore != 'ND' and countCatAfter != 'ND':      # Skip anything which lacks values
                                fileOut.write(countCatBefore + ',' + countCatAfter + '\n')

### Test 3
def compare_quiz_cats_to_csv(quizCatDict1, quizCatDict2, resultDir):
        import os
        fileNames = []
        for i in range(len(quizCatDict1)):
                #if i == 0:
                #        continue
                # Get the output directory
                outDirPath = os.path.join(os.getcwd(), resultDir)
                if not os.path.isdir(outDirPath):
                        os.mkdir(outDirPath)
                # Generate an output name
                name = file_name_gen(os.path.join(outDirPath, 'quizcat_totals_compare'), '_' + str(i+1) + '.csv')
                fileNames.append(name)
                # Produce output file
                with open(name, 'w') as fileOut:
                        fileOut.write('quiz_response,treatment,count_category\n')
                        for key1, value1 in quizCatDict1[i].items():
                                value2 = quizCatDict2[i][key1]
                                for x in range(len(value1)):
                                        fileOut.write(str(key1) + ',before,' + value1[x] + '\n')
                                        fileOut.write(str(key1) + ',after,' + value2[x] + '\n')

## Test 5
def compare_div_to_csv(chunkDict,  resultDir):
        import os
        # Get the output directory
        outDirPath = os.path.join(os.getcwd(), resultDir)
        if not os.path.isdir(outDirPath):
                os.mkdir(outDirPath)
        # Generate an output name
        name = file_name_gen(os.path.join(outDirPath, 'quizcat_div_overall'), '.csv')
        with open(name, 'w') as fileOut:
                fileOut.write('div_before,div_after\n')
                for key, value in chunkDict.items():
                        # Compute the diversity count
                        beforeDiv, afterDiv = diversity_number(value)
                        fileOut.write(str(beforeDiv) + ',' + str(afterDiv) + '\n')

## Test 6
def compare_quiz_cats_div_to_csv(quizCatDict1, quizCatDict2, resultDir):
        import os
        fileNames = []
        for i in range(len(quizCatDict1)):
                #if i == 0:
                #        continue
                # Get the output directory
                outDirPath = os.path.join(os.getcwd(), resultDir)
                if not os.path.isdir(outDirPath):
                        os.mkdir(outDirPath)
                # Generate an output name
                name = file_name_gen(os.path.join(outDirPath, 'quizcat_div_compare'), '_' + str(i+1) + '.csv')
                fileNames.append(name)
                # Produce output file
                with open(name, 'w') as fileOut:
                        fileOut.write('quiz_response,treatment,species_num\n')
                        for key1, value1 in quizCatDict1[i].items():
                                value2 = quizCatDict2[i][key1]
                                for x in range(len(value1)):
                                        fileOut.write(str(key1) + ',before,' + str(value1[x]) + '\n')
                                        fileOut.write(str(key1) + ',after,' + str(value2[x]) + '\n')

## Test 7
def media_quiz_cats_to_csv(quizCatDict, prefix, suffix, header, resultDir):
        import os
        mediaKeys = list(quizCatDict.keys())
        for x in mediaKeys:
                for i in range(len(quizCatDict[x])):
                        # Get the output directory
                        outDirPath = os.path.join(os.getcwd(), resultDir)
                        if not os.path.isdir(outDirPath):
                                os.mkdir(outDirPath)
                        # Generate an output name
                        name = file_name_gen(os.path.join(outDirPath, prefix + suffix), '_' + x + '_' + str(i+1) + '.csv')
                        # Produce output file
                        with open(name, 'w') as fileOut:
                                fileOut.write(header)
                                for key, value in quizCatDict[x][i].items():
                                        for val in value:
                                                fileOut.write(str(key) + ',' + str(val) + '\n')

## Test 8
def compare_media_to_csv(chunkDict,  resultDir):
        import os
        # Get the output directory
        outDirPath = os.path.join(os.getcwd(), resultDir)
        if not os.path.isdir(outDirPath):
                os.mkdir(outDirPath)
        # Generate an output name
        name = file_name_gen(os.path.join(outDirPath, 'quizcat_media_overall'), '.csv')
        with open(name, 'w') as fileOut:
                fileOut.write('media_before,media_after\n')
                for key, value in chunkDict.items():
                        # Compute the diversity count
                        beforeDiv, afterDiv = diversity_number(value)
                        fileOut.write(str(beforeDiv) + ',' + str(afterDiv) + '\n')

## Basic functions
def divide_num_to_list(totalNum, numGroups):
        baseNum = int(totalNum / numGroups)
        remainder = totalNum % numGroups
        valueList = [baseNum]*numGroups
        for i in range(remainder):
                valueList[i] += 1
        return valueList

def file_name_gen(prefix, suffix):
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + suffix):
                        return prefix + suffix
                elif os.path.isfile(prefix + suffix + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + suffix + str(ongoingCount)

def write_text_to_file(fileName, text):
        with open(fileName, 'w') as fileOut:
                fileOut.write(text)

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
args.tableFile = r'E:\elise\RELABEL_LQB301_ProjectData25-06-2018(9673).txt'
#args.tableFile = r'E:\elise\DELETED_LQB301_ProjectData25-06-2018(9673).txt'
args.outputFileName = r'E:\elise\test_micro_parse.txt'
validate_args(args)

## PARSE AND PROCESS RAW DATA

# Parse table file into chunks
chunks = micro_table_chunks(args.tableFile)

# Drop chunks that were hidden
chunks = drop_hidden_chunks(chunks)

# Turn chunks into dictionaries
chunkDict = micro_chunk_to_dict(chunks)

# Fix details for the dictionary
chunkDict = fix_details(chunkDict)

# REFORMAT DATA FOR VARIOUS ANALYSES

## Prep: count by quiz response dictionary generation
quizCatBefore, quizCatAfter = total_count_quiz_response(chunkDict)

## Prep: diversity by quiz response dictionary generation
quizCatDivBefore, quizCatDivAfter = microbe_diversity_quiz_response(chunkDict)

## Prep: media reaction by quiz response dictionary generation
quizCatMediaBefore, quizCatMediaAfter = media_reaction_quiz_response(chunkDict)

# GENERATE AN OVERVIEW OF THE DATA
tabulatedData = tabulate_data(chunkDict)
#write_text_to_file('micro_data_table.txt', tabulatedData)

# PERFORM TESTS

## Overall question: What information does comparing before and after total counts tell us?

### Test 1: growth before / growth after for quiz responses [Q: Do students who provide certain answers have more/less growth before and/or after treatment? A: For medical items, yes.]
#quiz_cats_to_csv(quizCatBefore, 'quizcat_totals_', 'before', 'quiz_response,count_category\n', 'R_scripts')
#quiz_cats_to_csv(quizCatAfter, 'quizcat_totals_', 'after', 'quiz_response,count_category\n', 'R_scripts')

### Test 2: growth before VERSUS after [Q: Does treatment work to kill microbes? A: Yes.]
#compare_total_to_csv(chunkDict, 'R_scripts')

### Test 3: growth before VERSUS after for quiz responses [Q: Is cleaning more effective/important for students who provide certain answers? - A: For medical items, yes.]
#compare_quiz_cats_to_csv(quizCatBefore, quizCatAfter, 'R_scripts')     ## This is the better option, it fees into linear model generation well

## Overall question: Is there any relationship between study variables and microbial diversity?

### Test 4: microbial diversity before / after for quiz responses [Q: Do students who provide certain answers have more/less microbial diversity before and/or after treatment? A: Not technically.]
#quiz_cats_to_csv(quizCatDivBefore, 'quizcat_div_', 'before', 'quiz_response,species_num\n', 'R_scripts')
#quiz_cats_to_csv(quizCatDivAfter, 'quizcat_div_', 'after', 'quiz_response,species_num\n', 'R_scripts')

### Test 5: microbial diversity before VERSUS after [Q: Does treatment work to reduce diversity? A: Yes.]
#compare_div_to_csv(chunkDict, 'R_scripts')

### Test 6: div before VERSUS after for quiz responses [Q: Does cleaning reduce diversity more for students who provide certain answers? - A: No.]
#compare_quiz_cats_div_to_csv(quizCatDivBefore, quizCatDivAfter, 'R_scripts')

## Overall question: Is there any relationship between study variables and media results?

### Test 7: media results before / after for quiz responses [Q: Are students who provide certain answers more likely to obtain specific media results before and/or after treatment? A: Yes.]
#media_quiz_cats_to_csv(quizCatMediaBefore, 'quizcat_media_', 'before', 'quiz_response,reaction_cat\n', 'R_scripts')
#media_quiz_cats_to_csv(quizCatMediaAfter, 'quizcat_media_', 'after', 'quiz_response,reaction_cat\n', 'R_scripts')

### NOT TESTING: media results before VERSUS after [Q: Does treatment work to reduce the occurrence of certain reactions? A: Not appropriate test.]

### Test 8: Subset - of the colonies that grow, are there substantive differences in the type of reaction identified based on quiz responses? [A: .]

#### SCRIPT ALL DONE
print('Program completed successfully!')
