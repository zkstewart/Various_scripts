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
                if key == 'CZ124' or key == 'CA152' or key == 'FJ000' or key == 'K123M':
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
                        # Skip the item for now - another function will deal with this
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
                       'UNSPECIFIED': 'UNSPECIFIED', }
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

### Test 1
def quiz_cats_to_csv(quizCatDict, suffix, resultDir):
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
                name = file_name_gen(os.path.join(outDirPath, 'quizcat_totals_' + suffix), '_' + str(i+1) + '.csv')
                fileNames.append(name)
                # Produce output file
                with open(name, 'w') as fileOut:
                        fileOut.write('quiz_response,count_category\n')
                        for key, value in quizCatDict[i].items():
                                for val in value:
                                        fileOut.write(str(key) + ',' + val + '\n')

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
                if i == 0:
                        continue
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

# PERFORM TESTS

## Overall question: What information does comparing before and after total counts tell us?

### Test 1: growth before / growth after for quiz responses [Q: Do students who provide certain answers have more/less growth before and/or after treatment? A: For medical items, yes.]
#quiz_cats_to_csv(quizCatBefore, 'before', 'R_scripts')
#quiz_cats_to_csv(quizCatAfter, 'after', 'R_scripts')

### Test 2: growth before VERSUS after [Q: Does treatment work to kill microbes? A: Yes.]
#compare_total_to_csv(chunkDict, 'R_scripts')

### Test 3: growth before VERSUS after for quiz responses [Q: Is cleaning more effective/important for students who provide certain answers? - A: No.]
#compare_quiz_cats_to_csv(quizCatBefore, quizCatAfter, 'R_scripts')     ## This is the better option, it fees into linear model generation well

### Test 4: growth before VERSUS growth after for medical items [Q: Are certain medical items more difficult to clean than others?]






# Junkyard
## Test 3
#compare_quiz_cats_to_csv_byresponse(quizCatBefore, quizCatAfter, 'R_scripts')
#def compare_quiz_cats_to_csv_byresponse(quizCatDict1, quizCatDict2, resultDir):
#        import os
#        for i in range(len(quizCatDict1)):
#                if i == 0:
#                        continue
#                # Get the output directory
#                outDirPath = os.path.join(os.getcwd(), resultDir)
#                if not os.path.isdir(outDirPath):
#                        os.mkdir(outDirPath)
#                # Produce output file per response
#                for key1, value1 in quizCatDict1[i].items():
#                        if len(value1) < 10:        # Skip values with insufficient response for proper statistics
#                                continue
#                        # Generate an output name
#                        name = file_name_gen(os.path.join(outDirPath, 'quizcat_totals_compare'), '_q' + str(i+1) + '_r' + str(key1) + '.csv')
#                        with open(name, 'w') as fileOut:
#                                fileOut.write('quiz_response,count_category_before,count_category_after\n')
#                                value2 = quizCatDict2[i][key1]
#                                for x in range(len(value1)):
#                                        fileOut.write(str(key1) + ',' + value1[x] + ',' + value2[x] + '\n')

#### SCRIPT ALL DONE
print('Program completed successfully!')
