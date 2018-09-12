#! python3

# fasta_handling_master_code.py

# Combination of multiple functions as a go-to script for handling
# and manipulating FASTA files.

# Load packages for main
import os, argparse, time

# Define functions
def fasta_ids(fastaFile):
        # Set up
        from Bio import SeqIO
        outList = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outList.append(record.id)
        return outList

def fasta_descriptions(fastaFile):
        # Set up
        from Bio import SeqIO
        outList = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outList.append(record.description)
        return outList

def fasta_lengths(fastaFile):
        # Set up
        from Bio import SeqIO
        outList = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outList.append(str(len(record)))
        return outList

def fasta_count(fastaFile):
        # Set up
        from Bio import SeqIO
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        ongoingCount = 0
        for record in records:
                ongoingCount += 1
        ongoingCount = [str(ongoingCount)]
        return ongoingCount

def fasta_multi2single(fastaFile):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outFasta.append('>' + record.description + '\n' + str(record.seq))
        return outFasta

def fasta_single2multi(fastaFile, multilineLength):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                sequence = str(record.seq)
                sequence = '\n'.join([sequence[i:i+multilineLength] for i in range(0, len(sequence), multilineLength)])
                outFasta.append('>' + record.description + '\n' + sequence)
        return outFasta

def fasta_cullbelow(fastaFile, length):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                sequence = str(record.seq)
                if len(sequence) < int(length):
                        continue
                # Output
                outFasta.append('>' + record.description + '\n' + sequence)
        return outFasta

def fasta_cullabove(fastaFile, length):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                sequence = str(record.seq)
                if len(sequence) > int(length):
                        continue
                # Output
                outFasta.append('>' + record.description + '\n' + sequence)
        return outFasta

## Fasta and fastq compatible functions
def fasta_rename(fastaFile, seqidprefix, prefix):
        # Set up
        from Bio import SeqIO
        outList = []
        outFasta = []
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        ongoingCount = 1
        for record in records:
                # Extract relevant details regardless of fasta or fastq
                if seqType == 'fasta':
                        seq = str(record.seq)
                        qual = ''
                else:
                        seq, qual = fastq_format_extract(record)
                # Main function action
                oldseqid = record.description
                newseqid = seqidprefix + str(ongoingCount)
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + newseqid + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + newseqid + '\n' + seq + '\n+\n' + qual) #fq
                outList.append(oldseqid + '\t' + newseqid)
                ongoingCount += 1
        return outList, outFasta, fastaFile, changed

def fasta_removestringfseqid(fastaFile, removeString, prefix):
        # Set up
        from Bio import SeqIO
        outFasta = []
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
                seqid = record.description
                seqid = seqid.replace(removeString, '')
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + seqid + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + seqid + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

def fasta_splitseqidatstring(fastaFile, splitString, prefix):
        # Set up
        from Bio import SeqIO
        outFasta = []
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
                seqid = record.description
                seqid = seqid.split(splitString)[0]
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + seqid + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + seqid + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

def fasta_trim(fastaFile, trimString, prefix):
        # Set up
        import re
        from Bio import SeqIO
        outFasta = []
        # Check the string to ensure it is correctly formatted
        stringCheck = re.compile(r'^[Ss](\d{1,10})[Ee](\d{1,10})$')
        result = stringCheck.match(trimString)
        if result == None:
                print('Format of the input trimming string is incorrect.')
                print('Call this program with -H or -HELP for assistance.')
                quit()
        # Parse the input string
        startTrim = int(result.group(1))
        endTrim = int(result.group(2))
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
                if endTrim != 0:
                        seq = seq[startTrim:-endTrim]
                        qual = qual[startTrim:-endTrim]
                else:
                        seq = seq[startTrim:]
                        qual = qual[startTrim:]
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + record.description + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + record.description + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

def fasta_retrieveseqwstring(fastaFile, retrieveString, prefix):
        # Set up
        from Bio import SeqIO
        outFasta = []
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
                if retrieveString not in seq:
                        continue
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + record.description + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + record.description + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

def fasta_retrieveseqidwstring(fastaFile, retrieveString, prefix):
        # Set up
        from Bio import SeqIO
        outFasta = []
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
                if retrieveString not in record.description:
                        continue
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + record.description + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + record.description + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

def fasta_removeseqwstring(fastaFile, removeString, prefix):
        # Set up
        from Bio import SeqIO
        outFasta = []
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
                if removeString in seq:
                        continue
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + record.description + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + record.description + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

def fasta_removeseqidwstring(fastaFile, removeString, prefix):
        # Set up
        from Bio import SeqIO
        outFasta = []
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
                if removeString in record.description:
                        continue
                # Output
                if seqType == 'fasta':
                        outFasta.append('>' + record.description + '\n' + seq)                  #fa
                else:
                        outFasta.append('@' + record.description + '\n' + seq + '\n+\n' + qual) #fq
        return outFasta, fastaFile, changed

def fasta_chunk(fastaFile, threads, prefix):
        # Set up
        import math
        from Bio import SeqIO
        outFasta = []
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Count number of sequences in file
        with open(fastaFile, 'r') as fileIn:
                if seqType == 'fasta':
                        numSeqs = 0
                        for line in fileIn:
                                if line.startswith('>'):
                                        numSeqs += 1
                else:
                        ongoingCount = 0
                        for line in fileIn:
                                ongoingCount += 1
                        numSeqs = ongoingCount / 4      # Fastq doesn't have multi-line, so it's easy to derive number of sequences this way
        # Find out where we are chunking the file
        rawNum = numSeqs / threads                              # In cases where threads > numSeqs, rawNum will be less than 1. numRoundedUp will equal the number of threads, and so we'll end up rounding these to 1. Yay!
        numRoundedUp = round((rawNum % 1) * threads, 0)         # By taking the decimal place and multiplying it by the num of threads, we can figure out how many threads need to be rounded up to process every sequence
        chunkPoints = []
        ongoingCount = 0
        for i in range(threads):
                if i+1 <= numRoundedUp:                 # i.e., if two threads are being rounded up, we'll round up the first two loops of this
                        chunkPoints.append(math.ceil(rawNum) + ongoingCount)    # Round up the rawNum, and also add our ongoingCount which corresponds to the number of sequences already put into a chunk
                        ongoingCount += math.ceil(rawNum)
                else:
                        chunkPoints.append(math.floor(rawNum) + ongoingCount)
                        ongoingCount += math.floor(rawNum)
                if ongoingCount >= numSeqs:             # Without this check, if we have more threads than sequences, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                        break                           # This doesn't actually affect program function, but for aesthetic reasons and for clarity of how this function works, I prevent this from occurring.
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        ongoingCount = 0    # This will keep track of what sequence number we are on
        for i in range(threads):
                # Flow control
                if ongoingCount == numSeqs: # This lets us stop making new files if we have more threads than we do sequences
                        break
                # Format sequences to list of lists
                outFasta.append([])
                for record in records:      # We'll run out of records if we get to a point where ongoingCount == numSeqs
                        # Extract relevant details regardless of fasta or fastq
                        if seqType == 'fasta':
                                seq = str(record.seq)
                                qual = ''
                        else:
                                seq, qual = fastq_format_extract(record)
                        # Output
                        if seqType == 'fasta':
                                outFasta[-1].append('>' + record.description + '\n' + seq)                  #fa
                        else:
                                outFasta[-1].append('@' + record.description + '\n' + seq + '\n+\n' + qual) #fq
                        # Main function action
                        ongoingCount += 1
                        if ongoingCount in chunkPoints:
                                break
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
                print('I don\'t recognise the input file! It should start with "@" (fastq) or ">" (fasta). Fix your inputs to continue.')
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
def validate_args(args, stringFunctions, numberFunctions, functionList):
        # Provide detailed help if specified
        if args.detailedHelp:
                import textwrap
                ## No input
                ids = '''
                The _ids_ function requires no special input. This function will
                produce an output text file listing sequence IDs (as parsed by Biopython
                .ids method)
                '''
                descriptions = '''
                The _descriptions_ function requires no special input. This function will
                produce an output text file listing sequence descriptions (as parsed by 
                Biopython .description method)
                '''
                lengths = '''
                The _lengths_ function requires no special input. This function 
                will produce an output text file listing the lengths of each sequence.
                '''
                count = '''
                The _count_ function requires no special input. This function 
                will produce an output text file with a single line depicting the number
                of sequences in the fasta file.
                '''
                ## Number input
                multi2single = '''
                The _multi2single_ function requires no special input. The output is 
                a singleline formatted fasta file.
                '''
                single2multi = '''
                The _single2multi_ function accepts a number input. This number refers
                to the number of sequence characters displayed per line. The output
                is a multiline formatted fasta file.
                '''
                cullbelow = '''
                The _cullbelow_ function accepts a number input. This number refers
                to the minimum length of sequence that will be present in the output
                fasta file.
                '''
                cullabove = '''
                The _cullabove_ function accepts a number input. This number refers
                to the maximum length of sequence that will be present in the output
                fasta file.
                '''
                chunk = '''
                The _chunk_ function accepts a number input. This number refers
                to the number of files to divide the original into as evenly as possible.
                Behaviour example: if you specify 10 chunks from a file with 8
                sequences, only 8 output files will be created.
                '''
                ## String input
                rename = '''
                The _rename_ function accepts a string input. This acts as prefix for 
                sequence IDs. For example, providing 'seq' as string input will result
                in an output fasta file containing 'seq1', 'seq2', etc., sequence IDs.
                '''
                removeseqwstring = '''
                The _removeseqwstring_ function accepts a string input. Any sequence
                which contains the specified string (case sensitive) will not be present
                in the output fasta file.
                '''
                removeseqidwstring = '''
                The _removeseqidwstring_ function accepts a string input. Any sequence ID
                which contains the specified string (case sensitive) will not be present
                in the output fasta file.
                '''
                retrieveseqwstring = '''
                The _retrieveseqwstring_ function accepts a string input. Only sequence IDs
                which contains the specified string (case sensitive) will be present
                in the output fasta file.
                '''
                retrieveseqidwstring = '''
                The _retrieveseqidwstring_ function accepts a string input. Only sequences
                which contains the specified string (case sensitive) will be present
                in the output fasta file.
                '''
                removestringfseqid = '''
                The _removestringfseqid_ function accepts a string input. This function will
                remove the specified string (case sensitive) from any sequence IDs
                and produce a new output fasta file.
                '''
                splitseqidatstring = '''
                The _splitseqidatstring_ function accepts a string input. This function will
                edit sequence IDs to remove any text that comes after the specified string
                (case sensitive) and produce a new output fasta file.
                '''
                trim = '''
                The _trim_ function accepts a string input in format s{digit}e{digit}. This
                function will remove {digit} from the start (s) and from the end (e)
                and produce a new output fasta file. Specify {digit} as 0 for no trimming
                to start or end - if {digit} is 1, the first or last character will be trimmed.
                '''
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
        if not os.path.isfile(args.fastaFileName):
                print('I am unable to locate the input fasta file (' + args.fastaFileName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle output file name & possibility that we are producing both list and fasta output
        outPrefix = args.outputFileName.rsplit('.', maxsplit=1)
        if len(outPrefix) == 1: # This probably means the user specified a prefix only; in this case we can get the suffix from the input file
                outSuffix = args.fastaFileName.rsplit('.', maxsplit=1)[-1]
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
                #None yet
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
stringFunctions = ['rename', 'removeseqwstring', 'removeseqidwstring', 'retrieveseqwstring', 'retrieveseqidwstring', 'removestringfseqid', 'splitseqidatstring', 'trim']
numberFunctions = ['single2multi', 'cullbelow', 'cullabove', 'chunk']
basicFunctions = ['ids', 'descriptions', 'lengths', 'count', 'multi2single']
functionList = stringFunctions + numberFunctions + basicFunctions

# Hold onto program 'start' time for the purpose of temporary file generation
startTime = time.time()
changed = False         # Default this as false; if we do create a temporary file this will become True

##### USER INPUT SECTION
usage = """%(prog)s handles fasta files, producing output according to the
selected function. For most functions, an input and output are all that is
required. A string or number is required for other functions. String input
is required for 'rename'. Number input is required for 'multi2single',
'single2multi'.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fastaFileName",
               help="Input fasta file")
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
listOutName, fastaOutName, args.outputFileName = validate_args(args, stringFunctions, numberFunctions, functionList)

# Enact functions
outList = []    # Blank lists so we can determine what needs to be output
outFasta = []
## String functions - FAST(A/Q) compatible - startTime is used for temporary file generation if the FASTQ file is faulty
if args.function == 'rename':
        outList, outFasta, args.fastaFileName, changed = fasta_rename(args.fastaFileName, args.string, startTime)
if args.function == 'removestringfseqid':
        outFasta, args.fastaFileName, changed  = fasta_removestringfseqid(args.fastaFileName, args.string, startTime)
if args.function == 'splitseqidatstring':
        outFasta, args.fastaFileName, changed  = fasta_splitseqidatstring(args.fastaFileName, args.string, startTime)
if args.function == 'retrieveseqwstring':
        outFasta, args.fastaFileName, changed = fasta_retrieveseqwstring(args.fastaFileName, args.string, startTime)
if args.function == 'retrieveseqidwstring':
        outFasta, args.fastaFileName, changed = fasta_retrieveseqidwstring(args.fastaFileName, args.string, startTime)
if args.function == 'removeseqwstring':
        outFasta, args.fastaFileName, changed = fasta_removeseqwstring(args.fastaFileName, args.string, startTime)
if args.function == 'removeseqidwstring':
        outFasta, args.fastaFileName, changed = fasta_removeseqidwstring(args.fastaFileName, args.string, startTime)
if args.function == 'trim':
        outFasta, args.fastaFileName, changed = fasta_trim(args.fastaFileName, args.string, startTime)
## Number functions
if args.function == 'single2multi':
        outFasta = fasta_single2multi(args.fastaFileName, args.number)
if args.function == 'cullbelow':
        outFasta = fasta_cullbelow(args.fastaFileName, args.number)
if args.function == 'cullabove':
        outFasta = fasta_cullabove(args.fastaFileName, args.number)
## Number functions - FAST(A/Q) compatible
if args.function == 'chunk':
        outFasta, args.fastaFileName, changed = fasta_chunk(args.fastaFileName, args.number, startTime)
## Basic functions
if args.function == 'ids':
        outList = fasta_ids(args.fastaFileName)
if args.function == 'descriptions':
        outList = fasta_descriptions(args.fastaFileName)  
if args.function == 'lengths':
        outList = fasta_lengths(args.fastaFileName)  
if args.function == 'count':
        outList = fasta_count(args.fastaFileName)  
if args.function == 'multi2single':
        outFasta = fasta_multi2single(args.fastaFileName)  

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
        if args.function == 'chunk':
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
