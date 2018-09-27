#! python3

# blast_handling_master_code.py

# Combination of multiple functions as a go-to script for handling
# and manipulating BLAST-tab files.

# Load packages for main
import os, argparse, time

# Define functions
def blast_hitcount(blastFile, evalue):
        # Set up
        from itertools import groupby
        grouper = lambda x: x.split('\t')[0]
        hitCount = 0
        totalCount = 0
        # Ensure that evalue is sensible
        if evalue < 0:
                print('Float value should not be less than 0 since it represents E-value (which is strictly +ve). Fix your input and try again.')
                quit()
        # Main function
        with open(blastFile) as fileIn:
                for key, value in groupby(fileIn, grouper):
                        bestHit = []
                        for entry in value:     # We do this process in case the file isn't sorted; this is the case for MMseqs2
                                sl = entry.rstrip('\r\n').split('\t')
                                if bestHit == []:
                                        bestHit = sl
                                else:
                                        # Check E-value
                                        if float(sl[10]) < float(bestHit[10]):
                                                bestHit = sl
                        # Check against cutoff
                        if float(bestHit[10]) <= evalue:
                                hitCount += 1
                        totalCount += 1
        return ['Total sequences with hits: ' + str(totalCount) + '. Number of sequences with significant hits: ' + str(hitCount) + '.'] # Return as a list so we can handle it like other functions

def blast_hitids(blastFile, evalue, whichToReturn):
        # Set up
        from itertools import groupby
        grouper = lambda x: x.split('\t')[0]
        hitList = []    # Scary function name
        # Ensure that evalue is sensible
        if evalue < 0:
                print('Float value should not be less than 0 since it represents E-value (which is strictly +ve). Fix your input and try again.')
                quit()
        # Ensure whichToReturn value is sensible
        if not whichToReturn.lower() in ['query', 'target', 'both']:
                print('hitids: String input with this function must be "query" or "target" or "both". Fix your input and try again.')
                quit()
        # Main function
        with open(blastFile) as fileIn:
                for key, value in groupby(fileIn, grouper):
                        bestHit = []
                        for entry in value:     # We do this process in case the file isn't sorted; this is the case for MMseqs2
                                sl = entry.rstrip('\r\n').split('\t')
                                if bestHit == []:
                                        bestHit = sl
                                else:
                                        # Check E-value
                                        if float(sl[10]) < float(bestHit[10]):
                                                bestHit = sl
                                        # If E-value is equivalent, check bit score
                                        elif float(sl[10]) == float(bestHit[10]):
                                                if float(sl[11]) > float(bestHit[11]):
                                                        bestHit = sl
                                                # If bit score is equivalent, check alignment length
                                                elif float(sl[11]) == float(bestHit[11]):
                                                        if int(sl[3]) > int(bestHit[3]):
                                                                bestHit = sl
                                                '''This is arbitrary, if E-value and bit score are identical, then the possible scenarios are
                                                that one hit is longer with lower identity, and the other hit is shorter with higher identity.
                                                I'm inclined to think that the first hit is better than the second if E-value/bit score are equivalent'''
                        # Process line to format it for output
                        if float(bestHit[10]) <= evalue:
                                if whichToReturn.lower() == 'query':
                                        hitList.append(sl[0])
                                elif whichToReturn.lower() == 'target':
                                        hitList.append(sl[1])
                                else:
                                        hitList += [sl[0], sl[1]]
        return hitList

def blast_fastahitretrieveremove(blastFile, fastaFile, evalue, behaviour, prefix):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Ensure that evalue is sensible
        if evalue < 0:
                print('Float value should not be less than 0 since it represents E-value (which is strictly +ve). Fix your input and try again.')
                quit()
        # Ensure that behaviour is sensible
        if not behaviour.lower() in ['retrieve', 'remove']:
                print('fastahitretrieveremove: String input with this function must be "retrieve" or "remove". Fix your input and try again.')
                quit()
        # Call hitids function to do the main work
        hitList = blast_hitids(blastFile, evalue, 'both')
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        for record in records:
                # Extract relevant details regardless of fasta or fastq
                if seqType == 'fasta':
                        seq = str(record.seq)
                        qual = ''
                else:
                        seq, qual = fastq_format_extract(record)
                # Main function action
                if behaviour.lower() == 'remove' and (record.description in hitList or record.id in hitList):
                        continue
                elif behaviour.lower() == 'retrieve' and (record.description not in hitList and record.id not in hitList):
                        continue
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + record.description + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + record.description + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

# Define general purpose functions
def fasta_or_fastq(fastaFile, prefix):
        changed = False
        # Get the first letter
        with open(fastaFile, 'r') as seqFile:
                for line in seqFile:
                        firstChar1 = line[0]
                        break
        # Check first letter to see if it conforms to fastq or fasta expected format
        if firstChar1 == '@':
                seqType = 'fastq'
                # Check the file to see if Biopython is likely to accept it
                fastaFile, changed = fastq_qual_fix(fastaFile, prefix)
        elif firstChar1 == '>':
                seqType = 'fasta'
        else:
                print('I don\'t recognise the input FASTA/Q file! It should start with "@" (fastq) or ">" (fasta). Fix your inputs to continue.')
                quit()
        # Return value
        return seqType, fastaFile, changed

def fastq_qual_fix(fastaFile, prefix):
        from Bio import SeqIO
        # Check for errors
        records = SeqIO.parse(fastaFile, 'fastq')
        try:
                for record in records:
                        break
                records.close   # Honestly not sure if this is necessary
                return [fastaFile, False]
        except ValueError:
                # Create a temporary file with modified quality lines
                tmpName = file_name_gen(str(prefix), '.fastq')
                with open(fastaFile, 'r') as fileIn, open(tmpName, 'w') as fileOut:
                        ongoingCount = 1
                        for line in fileIn:
                                if ongoingCount == 3:
                                        if not line.startswith('+'):
                                                print('Something is wrong with your fastq formatting.')
                                                print('Line number ' + str(ongoingCount) + ' (1-based) should be a comment line, but it doesn\'t start with \'+\'')
                                                print('Fix this file somehow and try again.')
                                                quit()
                                        fileOut.write('+\n')
                                else:
                                        fileOut.write(line.rstrip('\r\n') + '\n')
                                ongoingCount += 1
                                if ongoingCount == 5:
                                        ongoingCount = 1       # Reset our count to correspond to the new fastq entry
        return [tmpName, True]

def fastq_format_extract(fastqRecord):
        fqLines = fastqRecord.format('fastq').split('\n')
        fqSeq = fqLines[1]
        fqQual = fqLines[3]
        return fqSeq, fqQual

def file_name_gen(prefix, suffix):
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + '1' + suffix):
                        return prefix + '1' + suffix
                elif os.path.isfile(prefix + str(ongoingCount) + suffix):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + suffix

# Define functions for output
def output_list(outList, outputFileName):
        with open(outputFileName, 'w') as fileOut:
                fileOut.write('\n'.join(outList))

def output_list_of_lists(outList, outputPrefix, outputSuffix):
        for i in range(len(outList)):
                with open(outputPrefix + '_' + str(i+1) + '.' + outputSuffix, 'w') as fileOut:
                        fileOut.write('\n'.join(outList[i]))

# Define function for validating arguments
def validate_args(args, stringFunctions, numberFunctions, fastaFunctions, functionList):
        # Provide detailed help if specified
        if args.detailedHelp:
                import textwrap
                ## No input
                ## Number input
                hitcount = '''
                The _hitcount_ function accepts a number input. This number should
                be a float value (e.g., 0.1 or 1e-1) and refers to the E-value cutoff
                used for counting significant hits. The output is a string describing
                the total number of sequences with hits and number of sequences with
                significant hits.
                '''
                ## String input
                ## String & number input
                hitids = '''
                The _hitids_ function accepts a string and number input. String input
                should be either "query" or "target" or "both" to specify which sequence IDs to
                return. Number input should be a float value (e.g., 0.1 or 1e-1) and
                refers to the E-value cutoff used for counting significant hits. 
                The output is a text file listing sequence IDs which are matched
                significantly using the E-value cutoff.
                '''
                fastahitretrieveremove = '''
                The _fastahitretrieveremove_ function accepts a string and number input;
                it additionally requires FASTA input. String input should be either "retrieve"
                or "remove" to specify whether you want the output FASTA file to only contain
                significant hits ("retrieve") or to exclude sequences with significant hits
                ("remove"). Number input should be a float value (e.g., 0.1 or 1e-1) and
                refers to the E-value cutoff used for counting significant hits. The output
                is a FASTA file containing sequences according to specified criteria.
                Note that both query and target columns of the BLAST-tab file
                will be checked for IDs
                '''
                ## Format help message
                printList = str(functionList).replace("'", "")
                printList = eval(printList)
                for entry in printList:
                        entry = textwrap.dedent(entry)
                        entry = entry.strip('\n').replace('\n', ' ')
                        for line in textwrap.wrap(entry, width=50):
                                print(line)
                        print('')
                quit()
        # Validate file input
        if not os.path.isfile(args.blastFileName):
                print('I am unable to locate the input BLAST file (' + args.blastFileName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if args.function in fastaFunctions:
                if args.fastaFileName != None:
                        if not os.path.isfile(args.fastaFileName):
                                print('I am unable to locate the input fasta file (' + args.fastaFileName + ')')
                                print('Make sure you\'ve typed the file name or location correctly and try again.')
                                quit()
                else:
                        print('You need to specify a FASTA argument when running function \'' + args.function + '\'. Try again.')
                        quit()
        # Handle output file name & possibility that we are producing both list and fasta output
        outPrefix = args.outputFileName.rsplit('.', maxsplit=1)
        if len(outPrefix) == 1: # This probably means the user specified a prefix only; in this case we can get the suffix from the input file
                if args.fastaFileName != None:          # For the time being, if the user is performing a function that uses FASTA input we assume its output is also FASTA; this might change at some point
                        outSuffix = args.fastaFileName.rsplit('.', maxsplit=1)[-1]
                        args.outputFileName += '.' + outSuffix
                        outPrefix.append(outSuffix)
                else:
                        outSuffix = args.blastFileName.rsplit('.', maxsplit=1)[-1]
                        args.outputFileName += '.' + outSuffix
                        outPrefix.append(outSuffix)
        listOutName = outPrefix[0] + '_list.' + outPrefix[1]    # Technically this could be annoying at times if we aren't producing _list/_fasta output
        fastaOutName = outPrefix[0] + '_fasta.' + outPrefix[1]  # but it'd be more annoying to hard code each function to check if they have single or multiple outputs
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        if os.path.isfile(listOutName):
                print(listOutName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        if os.path.isfile(fastaOutName):
                print(fastaOutName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        # Validate string inputs if relevant
        if args.function in stringFunctions:
                if args.string == None:
                        print('You need to specify a string argument when running function \'' + args.function + '\'. Try again.')
                        quit()
        # Validate number inputs if relevant
        if args.function in numberFunctions:
                if args.number == None:
                        print('You need to specify a number argument when running function \'' + args.function + '\'. Try again.')
                        quit()
                # Float-based functions
                if args.function == 'hitcount' or args.function == 'hitids' or args.function == 'fastahitretrieveremove':
                        try:
                                args.number = float(args.number)
                        except:
                                print('The specified number is not accepted as a float. Check your input to make sure your input is correct (e.g., 5.00 or 1e-10) and try again.')
                                quit()
                # Integer-based functions
                if args.function == 'single2multi' or args.function == 'cullbelow' or args.function == 'cullabove' or args.function == 'chunk':
                        try:
                                args.number = int(args.number)
                        except:
                                print('The specified number is not accepted as an integer. Check your input to make sure this is a plain number (e.g., 5 or 100, not 5.00 or 1e-10) and try again.')
                                quit()
        return listOutName, fastaOutName, args.outputFileName

'''
To add a new function into this program, you need to 1) add a new description  
to the validate_args function, 2) handle it specifically if it is an integer 
or float, 3) add the actual function above, 4) add it into the function list,
and 5) enact the function below.
'''

# Function list - update as new ones are added
stringFunctions = ['hitids', 'fastahitretrieveremove']
numberFunctions = ['hitcount', 'hitids', 'fastahitretrieveremove']
basicFunctions = []
fastaFunctions = ['fastahitretrieveremove']
functionList = []
for fList in [stringFunctions + numberFunctions + basicFunctions + fastaFunctions]:     # Need more elaborate process since there can be duplicates in string and number lists
        for func in fList:
                if func not in functionList:
                        functionList.append(func)

# Hold onto program 'start' time for the purpose of temporary file generation
startTime = time.time()
changed = False         # Default this as false; if we do create a temporary file this will become True

##### USER INPUT SECTION
usage = """%(prog)s handles BLAST-tab files, producing output according to the
selected function. In addition to an input BLAST-tab file and output file name,
some functions require a string and/or number and/or FASTA fileinput as well.
Call program with -H tag for a detailed description of each function.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="blastFileName",
               help="Input BLAST-tab file")
p.add_argument("-if", "-inputFasta", dest="fastaFileName",
               help="Input FASTA file (if relevant)")
p.add_argument("-f", dest="function", choices=functionList,
               help="Function to run")
p.add_argument("-s", "-string", dest="string", type=str,
               help="String to use for various functions (if relevant)")
p.add_argument("-n", "-number", dest="number",
               help="Number to use for various functions (if relevant)")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name")
p.add_argument("-H", "-HELP", dest="detailedHelp", action='store_true',
             help="Provide detailed help for each function")

args = p.parse_args()
listOutName, fastaOutName, args.outputFileName = validate_args(args, stringFunctions, numberFunctions, fastaFunctions, functionList)

# Enact functions
outList = []    # Blank lists so we can determine what needs to be output
outFasta = []
## String functions
## String functions - FAST(A/Q) input
if args.function == 'fastahitretrieveremove':
        outFasta, args.fastaFileName, changed = blast_fastahitretrieveremove(args.blastFileName, args.fastaFileName, args.number, args.string, startTime)       # startTime is used for temporary file generation if the FASTQ file is faulty
## Number functions
if args.function == 'hitcount':
        outList = blast_hitcount(args.blastFileName, args.number)
## Number functions - FAST(A/Q) input
## Number & string functions
if args.function == 'hitids':
        outList = blast_hitids(args.blastFileName, args.number, args.string) 

# Output results
if outList != [] and outList != None:
        ## Special handling for specific functions
        # None yet
        if outFasta != [] and outFasta != None:
                output_list(outList, listOutName)
        else:
                output_list(outList, args.outputFileName)
if outFasta != [] and outFasta != None:
        ## Special handling for specific functions
        if args.function == 'chunk':    ## KEEP AS EXAMPLE - No functions need similar behaviour yet, but they might in the future
                fastaPrefix, fastaSuffix = args.outputFileName.rsplit('.', maxsplit=1)[0], args.outputFileName.rsplit('.', maxsplit=1)[1]
                output_list_of_lists(outFasta, fastaPrefix, fastaSuffix)
        elif outList != [] and outList != None:
                output_list(outFasta, fastaOutName)
        else:
                output_list(outFasta, args.outputFileName)

# Let the user know about null results
if (outList == [] or outList == None) and (outFasta == [] or outFasta == None):
        print('Looks like there is no output from this function. No output files will be generated.')

# Remove tmp files if relevant
if changed == True:
        print('Note that the output file is a bit different than the original.')
        print('In this case, the input fastq file had its description lines replaced with just the \'+\' character.') # For the time being, this print statement is correct, so we don't need specific reference to function name
        os.remove(args.fastaFileName)

print('Program completed successfully!')