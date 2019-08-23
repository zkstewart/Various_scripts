#! python3

# fasta_handling_master_code.py

# Combination of multiple functions as a go-to script for handling
# and manipulating FASTA files.

# Load packages for main
import os, argparse, time

# Define functions
## Fasta ONLY functions
def fasta_multi2single(fastaFile, outputFileName):
        # Set up
        from Bio import SeqIO
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        fastaOut.write('>' + record.description + '\n' + str(record.seq) + '\n')

def fasta_single2multi(fastaFile, multilineLength, outputFileName):
        # Set up
        from Bio import SeqIO
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        sequence = str(record.seq)
                        sequence = '\n'.join([sequence[i:i+multilineLength] for i in range(0, len(sequence), multilineLength)])
                        fastaOut.write('>' + record.description + '\n' + sequence + '\n')

## Fastq ONLY functions
def fasta_q_to_a(fastqFile, outputFileName):
        # Set up
        import os
        ongoingCount = 1
        cleanExit = False
        # Perform function
        with open(fastqFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
                for line in fileIn:
                        if ongoingCount == 1:
                                if not line.startswith('@'):
                                        print('Something is wrong with your fastq formatting.')
                                        print('Line number ' + str(ongoingCount) + ' (1-based) should be an ID line, but it doesn\'t start with \'@\'')
                                        print('Fix this file somehow and try again.')
                                        cleanExit = True
                                        break
                                fileOut.write('>' + line[1:])
                        elif ongoingCount == 2:
                                fileOut.write(line)
                        elif ongoingCount == 3:
                                if not line.startswith('+'):
                                        print('Something is wrong with your fastq formatting.')
                                        print('Line number ' + str(ongoingCount) + ' (1-based) should be a comment line, but it doesn\'t start with \'+\'')
                                        print('Fix this file somehow and try again.')
                                        cleanExit = True
                                        break
                        ongoingCount += 1
                        if ongoingCount == 5:
                                ongoingCount = 1       # Reset our count to correspond to the new fastq entry
        if cleanExit == True:
                os.unlink(outputFileName)

## Fasta and fastq compatible functions
def fasta_ids(fastaFile, prefix, outputFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        listOut.write(record.id + '\n')
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_descriptions(fastaFile, prefix, outputFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        listOut.write(record.description + '\n')
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_lengths(fastaFile, prefix, outputFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        listOut.write(str(len(record)) + '\n')
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_count(fastaFile, prefix, outputFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        ongoingCount = 0
        for record in records:
                ongoingCount += 1
        with open(outputFileName, 'w') as listOut:
                listOut.write(str(ongoingCount))
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_cullbelow(fastaFile, length, prefix, outputFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        # Extract relevant details regardless of fasta or fastq
                        if seqType == 'fasta':
                                seq = str(record.seq)
                                qual = ''
                        else:
                                seq, qual = fastq_format_extract(record)
                        if len(seq) < int(length):
                                continue
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_cullabove(fastaFile, length, prefix, outputFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        # Extract relevant details regardless of fasta or fastq
                        if seqType == 'fasta':
                                seq = str(record.seq)
                                qual = ''
                        else:
                                seq, qual = fastq_format_extract(record)
                        if len(seq) > int(length):
                                continue
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_rename(fastaFile, stringInput, prefix, outputFileName, listFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        ongoingCount = 1
        with open(outputFileName, 'w') as fastaOut, open(listFileName, 'w') as listOut:
                for record in records:
                        # Extract relevant details regardless of fasta or fastq
                        if seqType == 'fasta':
                                seq = str(record.seq)
                                qual = ''
                        else:
                                seq, qual = fastq_format_extract(record)
                        # Main function action
                        oldseqid = record.description
                        if '{}' in stringInput:
                                formatCount = stringInput.count('{}')
                                subList = [str(ongoingCount)]*formatCount
                                newseqid = stringInput.format(*subList)
                        else:
                                newseqid = stringInput + str(ongoingCount)
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + newseqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + newseqid + '\n' + seq + '\n+\n' + qual + '\n') #fq
                        listOut.write(oldseqid + '\t' + newseqid + '\n')
                        ongoingCount += 1
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_listrename(fastaFile, listFileName, prefix, outputFileName):
        # Set up
        import os
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Parse list file
        listDictL = {} # L = left
        listDictR = {} # R = right
        with open(listFileName, 'r') as listIn:
                for line in listIn:
                        sl = line.rstrip('\r\n').split('\t')
                        # Validate .list format
                        try:
                                assert len(sl) == 2
                        except:
                                print('.list file is not correctly formatted; should contain two columns separated by a tab.')
                                quit()
                        # Store data in both directions
                        listDictL[sl[0]] = sl[1]
                        listDictR[sl[1]] = sl[0]
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        # Extract relevant details regardless of fasta or fastq
                        if seqType == 'fasta':
                                seq = str(record.seq)
                                qual = ''
                        else:
                                seq, qual = fastq_format_extract(record)
                        # Main function action
                        oldseqid = record.description
                        if oldseqid in listDictL:
                                newseqid = listDictL[oldseqid]
                        elif oldseqid in listDictR:
                                newseqid = listDictR[oldseqid]
                        else:
                                print('"' + oldseqid + '" not found in .list file. These files should match up.')
                                quit()
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + newseqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + newseqid + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_removestringfseqid(fastaFile, removeString, prefix, outputFileName):
        # Set up
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
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
                                fastaOut.write('>' + seqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + seqid + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_splitseqidatstring(fastaFile, splitString, prefix, outputFileName):
        # Set up
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
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
                                fastaOut.write('>' + seqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + seqid + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_trim(fastaFile, trimString, prefix, outputFileName):
        # Set up
        import re
        from Bio import SeqIO
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
        with open(outputFileName, 'w') as fastaOut:
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
                                fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_retrieveseqwstring(fastaFile, retrieveString, prefix, outputFileName):
        # Set up
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
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
                                fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_retrieveseqidwstring(fastaFile, retrieveString, prefix, outputFileName):
        # Set up
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
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
                                fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_removeseqwstring(fastaFile, removeString, prefix, outputFileName):
        # Set up
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
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
                                fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_removeseqidwstring(fastaFile, removeString, prefix, outputFileName):
        # Set up
        from Bio import SeqIO
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
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
                                fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

def fasta_chunk(fastaFile, threads, prefix, outputFileName):
        # Set up
        import math
        from Bio import SeqIO
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
                loopFileName = outputFileName.rsplit('.', maxsplit=1)[0] + '_chunk' + str(i+1) + '.' + outputFileName.rsplit('.', maxsplit=1)[1]
                with open(loopFileName, 'w') as fastaOut:
                        for record in records:      # We'll run out of records if we get to a point where ongoingCount == numSeqs
                                # Extract relevant details regardless of fasta or fastq
                                if seqType == 'fasta':
                                        seq = str(record.seq)
                                        qual = ''
                                else:
                                        seq, qual = fastq_format_extract(record)
                                # Output
                                if seqType == 'fasta':
                                        fastaOut.write('>' + record.description + '\n' + seq + '\n')                  #fa
                                else:
                                        fastaOut.write('@' + record.description + '\n' + seq + '\n+\n' + qual + '\n') #fq
                                # Main function action
                                ongoingCount += 1
                                if ongoingCount in chunkPoints:
                                        break
        # Clean temp file if relevant
        if changed == True:
                os.remove(fastaFile)

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
                # Alert user to behaviour
                print('Your input FASTQ file has qual lines (+ lines) which don\'t match the ID lines (@ lines).')
                print('This program will produce a temporary file (' + str(prefix) + ') with this issue fixed, assuming that qual lines which follow ID lines are related to each other.')
                print('If this isn\'t true, you need to fix your FASTQ file manually.')
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
                q_to_a = '''The _q_to_a_ function requires no special input. This function
                will produce an output fasta file when a fastq input is provided.
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
                The _rename_ function accepts a string input. This can behave in two ways.
                If you specify {} anywhere within the string, this will be substituted with
                an iterating number from 1-> inf. Example: string input of '{}seq' will 
                become '1seq', '2seq', etc. Not specifying {} will make this string act
                as prefix and the iterating number will be added to the end of it e.g.,
                string input of 'seq' will become 'seq1, 'seq2', etc.
                '''
                listrename = '''
                The _listrename_ function accepts a string input. This string should correspond
                to a file which has .list format (as produced by _rename_) and will swap
                sequence IDs within the provided file.
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
        if args.outputFileName == None:
                print('-o argument must be provided, fix your inputs and try again.')
                quit()
        outPrefix = args.outputFileName.rsplit('.', maxsplit=1)
        if len(outPrefix) == 1: # This probably means the user specified a prefix only; in this case we can get the suffix from the input file
                outSuffix = args.fastaFileName.rsplit('.', maxsplit=1)[-1]
                args.outputFileName += '.' + outSuffix
                outPrefix.append(outSuffix)
        listOutName = outPrefix[0] + '.list'
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        if os.path.isfile(listOutName):
                print(listOutName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        # Validate string inputs if relevant
        if args.function in stringFunctions:
                if args.string == None:
                        print('You need to specify a string argument when running function \'' + args.function + '\'. Try again.')
                        quit()
                # String-based functions
                if args.function == 'listrename':
                        if not os.path.isfile(args.string):
                                print('The specified string does not point to a file. Make sure you have typed this correctly or provided the full path and try again.')
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
        return listOutName, args.outputFileName

'''
To add a new function into this program, you need to 1) add a new description  
to the validate_args function, 2) handle it specifically if it is a special string, integer 
or float-based function, 3) add the actual function above, 4) add it into the function list,
and 5) enact the function below.
'''

# Function list - update as new ones are added
stringFunctions = ['rename', 'listrename', 'removeseqwstring', 'removeseqidwstring', 'retrieveseqwstring', 'retrieveseqidwstring', 'removestringfseqid', 'splitseqidatstring', 'trim']
numberFunctions = ['single2multi', 'cullbelow', 'cullabove', 'chunk']
basicFunctions = ['ids', 'descriptions', 'lengths', 'count', 'multi2single', 'q_to_a']
functionList = stringFunctions + numberFunctions + basicFunctions

# Hold onto program 'start' time for the purpose of temporary file generation
startTime = time.time()
changed = False         # Default this as false; if we do create a temporary file this will become True

##### USER INPUT SECTION
usage = """%(prog)s handles FASTA/Q files, producing output according to the
selected function. For most functions, an input and output are all that is
required; some require a string and/or number input as well. Call program
with -H tag for a detailed description of each function.
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
listOutName, args.outputFileName = validate_args(args, stringFunctions, numberFunctions, functionList)

# Enact functions
'''Note: startTime is used for temporary file generation if the FASTQ file is faulty'''
## String functions
if args.function == 'rename':
        fasta_rename(args.fastaFileName, args.string, startTime, args.outputFileName, listOutName)
if args.function == 'listrename':
        fasta_listrename(args.fastaFileName, args.string, startTime, args.outputFileName)     
if args.function == 'removestringfseqid':
        fasta_removestringfseqid(args.fastaFileName, args.string, startTime, args.outputFileName)
if args.function == 'splitseqidatstring':
        fasta_splitseqidatstring(args.fastaFileName, args.string, startTime, args.outputFileName)
if args.function == 'retrieveseqwstring':
        fasta_retrieveseqwstring(args.fastaFileName, args.string, startTime, args.outputFileName)
if args.function == 'retrieveseqidwstring':
        fasta_retrieveseqidwstring(args.fastaFileName, args.string, startTime, args.outputFileName)
if args.function == 'removeseqwstring':
        fasta_removeseqwstring(args.fastaFileName, args.string, startTime, args.outputFileName)
if args.function == 'removeseqidwstring':
        fasta_removeseqidwstring(args.fastaFileName, args.string, startTime, args.outputFileName)
if args.function == 'trim':
        fasta_trim(args.fastaFileName, args.string, startTime, args.outputFileName)
## Number functions
if args.function == 'single2multi':
        fasta_single2multi(args.fastaFileName, args.number, args.outputFileName)
if args.function == 'cullbelow':
        fasta_cullbelow(args.fastaFileName, args.number, startTime, args.outputFileName)
if args.function == 'cullabove':
        fasta_cullabove(args.fastaFileName, args.number, startTime, args.outputFileName)
## Number functions - FAST(A/Q) compatible
if args.function == 'chunk':
        fasta_chunk(args.fastaFileName, args.number, startTime, args.outputFileName)
## Basic functions
if args.function == 'ids':
        fasta_ids(args.fastaFileName, startTime, args.outputFileName)
if args.function == 'descriptions':
        fasta_descriptions(args.fastaFileName, startTime, args.outputFileName)
if args.function == 'lengths':
        fasta_lengths(args.fastaFileName, startTime, args.outputFileName)
if args.function == 'count':
        fasta_count(args.fastaFileName, startTime, args.outputFileName)
if args.function == 'multi2single':
        fasta_multi2single(args.fastaFileName, args.outputFileName)
if args.function == 'q_to_a':
        fasta_q_to_a(args.fastaFileName, args.outputFileName)
print('Program completed successfully!')
