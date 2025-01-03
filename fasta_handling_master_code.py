#! python3

# fasta_handling_master_code.py

# Combination of multiple functions as a go-to script for handling
# and manipulating FASTA files.

# Load packages for main
import os, argparse, re, math, textwrap
from Bio import SeqIO

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Define general purpose functions
def fasta_or_fastq(fastaFile):
        # Get the first letter
        with open(fastaFile, 'r') as seqFile:
                for line in seqFile:
                        firstChar1 = line[0]
                        break
        # Check first letter to see if it conforms to fastq or fasta expected format
        if firstChar1 == '@':
                seqType = 'fastq'
        elif firstChar1 == '>':
                seqType = 'fasta'
        else:
                print('I don\'t recognise the input file! It should start with "@" (fastq) or ">" (fasta). Fix your inputs to continue.')
                quit()
        # Return value
        return seqType

def AltFastqGeneralIterator(handle):
        '''
        Note: I have taken this code from Biopython's functions
        (https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py)
        and turned off the requirement for second_title to be just '+' or to be
        identical to the title_line value. This can be "monkey patched" into
        SeqIO by this call "SeqIO.QualityIO.FastqGeneralIterator = AltFastqGeneralIterator"
        if SeqIO was originally imported by "from Bio import SeqIO"
        '''
        # We need to call handle.readline() at least four times per record,
        # so we'll save a property look up each time:
        handle_readline = handle.readline

        line = handle_readline()
        if not line:
                return  # Premature end of file, or just empty?
        if isinstance(line[0], int):
                raise ValueError("Is this handle in binary mode not text mode?")

        while line:
                if line[0] != "@":
                        raise ValueError(
                                "Records in Fastq files should start with '@' character")
                title_line = line[1:].rstrip()
                # Will now be at least one line of quality data - in most FASTQ files
                # just one line! We therefore use string concatenation (if needed)
                # rather using than the "".join(...) trick just in case it is multiline:
                seq_string = handle_readline().rstrip()
                # There may now be more sequence lines, or the "+" quality marker line:
                while True:
                        line = handle_readline()
                        if not line:
                                raise ValueError("End of file without quality information.")
                        if line[0] == "+":
                                # The title here is optional, but if present must match!
                                ## My change is below
                                #second_title = line[1:].rstrip()
                                #if second_title and second_title != title_line: 
                                        #raise ValueError("Sequence and quality captions differ.")
                                break
                        seq_string += line.rstrip()  # removes trailing newlines
                # This is going to slow things down a little, but assuming
                # this isn't allowed we should try and catch it here:
                if " " in seq_string or "\t" in seq_string:
                        raise ValueError("Whitespace is not allowed in the sequence.")
                seq_len = len(seq_string)

                # Will now be at least one line of quality data...
                quality_string = handle_readline().rstrip()
                # There may now be more quality data, or another sequence, or EOF
                while True:
                        line = handle_readline()
                        if not line:
                                break  # end of file
                        if line[0] == "@":
                                # This COULD be the start of a new sequence. However, it MAY just
                                # be a line of quality data which starts with a "@" character.  We
                                # should be able to check this by looking at the sequence length
                                # and the amount of quality data found so far.
                                if len(quality_string) >= seq_len:
                                        # We expect it to be equal if this is the start of a new record.
                                        # If the quality data is longer, we'll raise an error below.
                                        break
                                # Continue - its just some (more) quality data.
                        quality_string += line.rstrip()

                if seq_len != len(quality_string):
                        raise ValueError("Lengths of sequence and quality values differs "
                                         " for %s (%i and %i)."
                                         % (title_line, seq_len, len(quality_string)))

                # Return the record and then continue...
                yield (title_line, seq_string, quality_string)

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

# Define quick FASTA/Q parser
class FastxParser:
    def __init__(self, file):
        self.file = file
        self.seqType = fasta_or_fastq(file)
    
    def __iter__(self):
        if self.seqType == 'fasta':
            with open(self.file, 'r') as fastaFile:
                for title, seq in SimpleFastaParser(fastaFile):
                    yield title, seq
        elif self.seqType == 'fastq':
            with open(self.file, 'r') as fastqFile:
                for title, seq, qual in FastqGeneralIterator(fastqFile):
                    yield title, seq, qual

# Define functions
## Fasta ONLY functions
def multi2single(fastaFile, outputFileName):
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        fastaOut.write('>' + record.description + '\n' + str(record.seq) + '\n')

def single2multi(fastaFile, multilineLength, outputFileName):
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        sequence = str(record.seq)
                        sequence = '\n'.join([sequence[i:i+multilineLength] for i in range(0, len(sequence), multilineLength)])
                        fastaOut.write('>' + record.description + '\n' + sequence + '\n')

def striphyphens(fastaFile, outputFileName):
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        fastaOut.write('>' + record.description + '\n' + str(record.seq).replace("-", "") + '\n')

def reversecomplement(fastaFile, outputFileName):
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        fastaOut.write('>' + record.description + '\n' + str(record.seq.reverse_complement()) + '\n')

def reversecomplement2multi(fastaFile, multilineLength, outputFileName):
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                for record in records:
                        sequence = str(record.seq.reverse_complement())
                        sequence = '\n'.join([sequence[i:i+multilineLength] for i in range(0, len(sequence), multilineLength)])
                        fastaOut.write('>' + record.description + '\n' + sequence + '\n')

def explodeintocontigs(fastaFile, outputDirectory):
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Create output directory
        if not os.path.exists(outputDirectory):
                print(f"explodeintocontigs is creating output directory '{outputDirectory}'...")
                os.makedirs(outputDirectory, exist_ok=True)
        # Perform function
        for record in records:
                outputFileName = os.path.join(outputDirectory, record.id + ".fasta")
                if os.path.exists(outputFileName):
                        raise FileExistsError('File ' + outputFileName + ' already exists. Please remove it and try again.')
                with open(outputFileName, 'w') as fastaOut:
                        fastaOut.write('>' + record.description + '\n' + str(record.seq) + '\n')

## Fastq ONLY functions
def q_to_a(fastqFile, outputFileName):
        # Set up
        ongoingCount = 0
        cleanExit = False
        # Perform function
        with open(fastqFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
                for line in fileIn:
                        if ongoingCount % 4 == 0:
                                if not line.startswith('@'):
                                        print('Something is wrong with your fastq formatting.')
                                        print('Line number ' + str(ongoingCount + 1) + ' (1-based) should be an ID line, but it doesn\'t start with \'@\'')
                                        print('Fix this file somehow and try again.')
                                        cleanExit = True
                                        break
                                fileOut.write('>' + line[1:])
                        elif ongoingCount % 4 == 1:
                                fileOut.write(line)
                        elif ongoingCount % 4 == 2:
                                if not line.startswith('+'):
                                        print('Something is wrong with your fastq formatting.')
                                        print('Line number ' + str(ongoingCount + 1) + ' (1-based) should be a comment line, but it doesn\'t start with \'+\'')
                                        print('Fix this file somehow and try again.')
                                        cleanExit = True
                                        break
                        elif ongoingCount % 4 == 3:
                                fileOut.write(line)
                        ongoingCount += 1
        if cleanExit == True:
                os.unlink(outputFileName)

## Fasta and fastq compatible functions
def crispressoalleles(fastaFile):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        seqs = []
        for record in records:
                seqs.append(str(record.seq).replace("-", ""))
        print(",".join(seqs))
        # Immediately exit to prevent contamination of stdout
        quit()

def crispressonames(fastaFile):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        names = []
        for record in records:
                names.append(record.id)
        print(",".join(names))
        # Immediately exit to prevent contamination of stdout
        quit()

def ids(fastaFile, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        listOut.write(record.id + '\n')

def descriptions(fastaFile, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        listOut.write(record.description + '\n')

def lengths(fastaFile, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        listOut.write(str(len(record)) + '\n')

def lengths_tsv(fastaFile, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        listOut.write(record.id + "\t" + str(len(record)) + '\n')

def gc(fastaFile, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        with open(outputFileName, 'w') as listOut:
                for record in records:
                        gcCount = [1 if base in ["C", "G"] else 0 for base in str(record.seq).upper()]
                        gcPct = sum(gcCount) / len(gcCount)
                        listOut.write(str(gcPct) + '\n')

def count(fastaFile, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        ongoingCount = 0
        for record in records:
                ongoingCount += 1
        with open(outputFileName, 'w') as listOut:
                listOut.write(str(ongoingCount))

def cullbelow(fastaFile, length, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                # Use performant parser
                parser = FastxParser(fastaFile)
                for seqTuple in parser:
                        # Check length
                        if len(seqTuple[1]) < int(length):
                                continue
                        # Output if length is met
                        if parser.seqType == 'fasta':
                                fastaOut.write('>' + seqTuple[0] + '\n' + seqTuple[1] + '\n')
                        else:
                                fastaOut.write('@' + seqTuple[0] + '\n' + seqTuple[1] + '\n+\n' + seqTuple[2] + '\n')

def cullabove(fastaFile, length, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Perform function
        with open(outputFileName, 'w') as fastaOut:
                # Use performant parser
                parser = FastxParser(fastaFile)
                for seqTuple in parser:
                        # Check length
                        if len(seqTuple[1]) > int(length):
                                continue
                        # Output if length is met
                        if parser.seqType == 'fasta':
                                fastaOut.write('>' + seqTuple[0] + '\n' + seqTuple[1] + '\n')
                        else:
                                fastaOut.write('@' + seqTuple[0] + '\n' + seqTuple[1] + '\n+\n' + seqTuple[2] + '\n')

def rename(fastaFile, stringInput, outputFileName, listFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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

def appendrename(fastaFile, stringInput, outputFileName, listFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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
                        oldseqid = record.id
                        newseqid = f"{oldseqid}{stringInput}"
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + newseqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + newseqid + '\n' + seq + '\n+\n' + qual + '\n') #fq
                        listOut.write(oldseqid + '\t' + newseqid + '\n')
                        ongoingCount += 1

def listrename(fastaFile, listFileName, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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
                        oldseqdesc = record.description
                        oldseqid = record.id
                        if oldseqid in listDictL:
                                newseqid = listDictL[oldseqid]
                        elif oldseqid in listDictR:
                                newseqid = listDictR[oldseqid]
                        elif oldseqdesc in listDictL:
                                newseqid = listDictL[oldseqdesc]
                        elif oldseqdesc in listDictR:
                                newseqid = listDictR[oldseqdesc]
                        else:
                                print('"' + oldseqid + '" not found in .list file. These files should match up.')
                                quit()
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + newseqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + newseqid + '\n' + seq + '\n+\n' + qual + '\n') #fq

def removestringfseqid(fastaFile, removeString, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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

def splitseqidatstring_start(fastaFile, splitString, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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
                        if splitString in seqid:
                                seqid = seqid.split(splitString, maxsplit=1)[0]
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + seqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + seqid + '\n' + seq + '\n+\n' + qual + '\n') #fq

def splitseqidatstring_end(fastaFile, splitString, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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
                        if splitString in seqid:
                                seqid = seqid.split(splitString, maxsplit=1)[1]
                        # Output
                        if seqType == 'fasta':
                                fastaOut.write('>' + seqid + '\n' + seq + '\n')                  #fa
                        else:
                                fastaOut.write('@' + seqid + '\n' + seq + '\n+\n' + qual + '\n') #fq

def trim(fastaFile, trimString, outputFileName):
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
        seqType = fasta_or_fastq(fastaFile)
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

def retrieveseqwstring(fastaFile, retrieveString, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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

def retrieveseqidwstring(fastaFile, retrieveString, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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

def removeseqwstring(fastaFile, removeString, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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

def removeseqidwstring(fastaFile, removeString, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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

def chunk(fastaFile, threads, outputFileName):
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
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

def twofastaseqidcompare(fastaFile1, fastaFile2, outputFileName):
        # Check for file type
        seqType1 = fasta_or_fastq(fastaFile1)
        seqType2 = fasta_or_fastq(fastaFile2)
        # Load fast(a/q) file
        records1 = SeqIO.parse(open(fastaFile1, 'r'), seqType1)
        records2 = SeqIO.parse(open(fastaFile2, 'r'), seqType2)
        # Obtain sequence IDs from files
        r1IDs = {}
        r2IDs = {}
        for r1 in records1:
                r1IDs.setdefault(r1.description)
        for r2 in records2:
                r2IDs.setdefault(r2.description)
        # Perform function
        r1only, r1absent, r2only, r2absent = [], [], [], []
        for k1 in r1IDs.keys():
                if k1 not in r2IDs:
                        r1only.append(k1)
                        r2absent.append(k1)
        for k2 in r2IDs.keys():
                if k2 not in r1IDs:
                        r2only.append(k2)
                        r1absent.append(k2)
        # Write output files
        with open(outputFileName + '_file1only.txt', 'w') as fileOut1:
                fileOut1.write('\n'.join(r1only))
        with open(outputFileName + '_file1absent.txt', 'w') as fileOut2:
                fileOut2.write('\n'.join(r1absent))
        with open(outputFileName + '_file2only.txt', 'w') as fileOut3:
                fileOut3.write('\n'.join(r2only))
        with open(outputFileName + '_file2absent.txt', 'w') as fileOut4:
                fileOut4.write('\n'.join(r2absent))

def twofastaseqidcompare_orthofinder(fastaFile1, fastaFile2, outputFileName):
        def sanitise_seqid(seqid):
                seqid = seqid.replace(',', '_')
                return seqid
        # Check for file type
        seqType1 = fasta_or_fastq(fastaFile1)
        seqType2 = fasta_or_fastq(fastaFile2)
        # Load fast(a/q) file
        records1 = SeqIO.parse(open(fastaFile1, 'r'), seqType1)
        records2 = SeqIO.parse(open(fastaFile2, 'r'), seqType2)
        # Obtain sequence IDs from files
        r1IDs = {}
        r2IDs = {}
        for r1 in records1:
                r1IDs.setdefault(sanitise_seqid(r1.description))
        for r2 in records2:
                r2IDs.setdefault(sanitise_seqid(r2.description))
        # Perform function
        r1only, r1absent, r2only, r2absent = [], [], [], []
        for k1 in r1IDs.keys():
                if k1 not in r2IDs:
                        r1only.append(k1)
                        r2absent.append(k1)
        for k2 in r2IDs.keys():
                if k2 not in r1IDs:
                        r2only.append(k2)
                        r1absent.append(k2)
        # Write output files
        with open(outputFileName + '_file1only.txt', 'w') as fileOut1:
                fileOut1.write('\n'.join(r1only))
        with open(outputFileName + '_file1absent.txt', 'w') as fileOut2:
                fileOut2.write('\n'.join(r1absent))
        with open(outputFileName + '_file2only.txt', 'w') as fileOut3:
                fileOut3.write('\n'.join(r2only))
        with open(outputFileName + '_file2absent.txt', 'w') as fileOut4:
                fileOut4.write('\n'.join(r2absent))

def mergefasta(fastaFile1, fastaFile2, outputFileName):
        # Check for file type
        seqType1 = fasta_or_fastq(fastaFile1)
        seqType2 = fasta_or_fastq(fastaFile2)
        # Load fast(a/q) file
        records1 = SeqIO.parse(open(fastaFile1, 'r'), seqType1)
        records2 = SeqIO.parse(open(fastaFile2, 'r'), seqType2)
        # Obtain sequence IDs from files
        r1IDs = {}
        r2IDs = {}
        for r1 in records1:
                r1IDs.setdefault(r1.description)
        for r2 in records2:
                r2IDs.setdefault(r2.description)
        # Identify unique and shared sequences from the two files
        '''We only need to check records2 since we can put record1
        completely into the merged output file
        '''
        r2only = {}
        for k2 in r2IDs.keys():
                if k2 not in r1IDs:
                        r2only.setdefault(k2)
        # Reload fast(a/q) file
        records1 = SeqIO.parse(open(fastaFile1, 'r'), seqType1)
        records2 = SeqIO.parse(open(fastaFile2, 'r'), seqType2)
        # Write output file
        with open(outputFileName, 'w') as fileOut:
                for record in records1:
                        fileOut.write(">{0}\n{1}\n".format(record.description, str(record.seq)))
                for record in records2:
                        if record.description in r2only:
                                fileOut.write(">{0}\n{1}\n".format(record.description, str(record.seq)))

def echoindex(fastaFile, index):
        # Validate index value
        if index < 1:
                raise ValueError('Index value must be 1 or greater.')
        # Check for file type
        seqType = fasta_or_fastq(fastaFile)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        ongoingCount = 0
        for record in records:
                ongoingCount += 1
                
                # Extract relevant details regardless of fasta or fastq
                if seqType == 'fasta':
                        seq = str(record.seq)
                else:
                        seq, _ = fastq_format_extract(record)
                
                # Print the sequence if it matches the index
                if ongoingCount == index:
                        print(seq)
                        # Immediately exit to prevent contamination of stdout
                        quit()
        
        # If we didn't find the index, print an error message
        raise ValueError(f"Input file '{fastaFile}' does not contain {index} sequences.")

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
                lengths_tsv = '''
                The _lengths_tsv_ function requires no special input. This function 
                will produce an output TSV pairing ID to the length of each sequence.
                '''
                count = '''
                The _count_ function requires no special input. This function 
                will produce an output text file with a single line depicting the number
                of sequences in the fasta file.
                '''
                q_to_a = '''The _q_to_a_ function requires no special input. This function
                will produce an output fasta file when a fastq input is provided.
                '''
                striphyphens = '''
                The _striphyphens_ function requires no special input. This function
                will produce an output fasta file sans any sequence hyphen characters.
                '''
                gc = '''
                The _gc_ function requires no special input. This function will produce
                and output text file listing the GC percentage of each sequence.
                '''
                explodeintocontigs = '''
                The _explodeintocontigs_ function requires no special input. This function
                will produce one output fasta file for each contig, with the output file name
                being the sequence ID. The output value for this function will specify a directory.
                '''
                crispressoalleles = '''
                The _crispressoalleles_ functions requires no special input. This function
                will take a FASTA/Q file and render into stdout a comma-delimited list of
                allele sequences suitable for input to CRISPResso2.
                '''
                crispressonames = '''
                The _crispressonames_ functions requires no special input. This function
                will take a FASTA/Q file and render into stdout a comma-delimited list of
                sequence IDs suitable for input to CRISPResso2.
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
                reversecomplement = '''
                The _reversecomplement_ function requires no special input. The output is
                a singleline formatted fasta file in which every sequence has been reverse
                complemented.
                '''
                reversecomplement2multi = '''
                The _reversecomplement2multi_ function accepts a number input. This number refers
                to the number of sequence characters displayed per line. The output is
                a multiline formatted fasta file in which every sequence has been reverse
                complemented.
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
                echoindex = '''
                The _echoindex_ function accepts a number input. This number refers to the
                index of the sequence to be echoed to stdout. The index should be given
                as a positive integer and will be 1-indexed. The result will be printed
                into stdout.
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
                appendrename = '''
                The _appendrename_ function accepts a string input. It will simply append
                the provided string to the existing sequence ID. For example, an existing
                ID of '>scaffold100' will become '>scaffold100{string}' without any
                space inbetween.
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
                splitseqidatstring_start = '''
                The _splitseqidatstring_start_ function accepts a string input. This function will
                edit sequence IDs to remove any text that comes BEFORE and including the specified
                string (case sensitive) and produce a new output fasta file.
                '''
                splitseqidatstring_end = '''
                The _splitseqidatstring_end_ function accepts a string input. This function will
                edit sequence IDs to remove any text that comes AFTER and including the specified
                string (case sensitive) and produce a new output fasta file.
                '''
                ### Special strings
                listrename = '''
                The _listrename_ function accepts a string input. This string should correspond
                to a file which has .list format (as produced by _rename_) and will swap
                sequence IDs within the provided file.
                '''
                twofastaseqidcompare = '''
                The _twofastaseqidcompare_ function accepts a string
                input. This string should correspond to a second fast(a/q) file whose sequence
                IDs are to be compared to the input file to find IDs present and absent in
                each file and produce text files listing these occurrences.
                '''
                twofastaseqidcompare_orthofinder = '''
                The _twofastaseqidcompare_orthofinder_ function accepts a string
                input. This string should correspond to a second fast(a/q) file whose sequence
                IDs are to be compared to the input file to find IDs present and absent in
                each file and produce text files listing these occurrences. This differs from
                the above by correcting changes induced by orthofinder to sequence IDs
                '''
                mergefasta = '''
                The _mergefasta_ function accepts a string input. This string should correspond
                to a second fast(a/q) file. These two files will be merged together non-redundantly
                on the basis of sequence IDs to produce a single output FASTA.
                '''
                ## Combined string & number input
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
        if args.outputFileName == None and args.function not in ["echoindex", "crispressoalleles", "crispressonames"]:
                print('-o argument must be provided, fix your inputs and try again.')
                quit()
        
        exclusions = ["explodeintocontigs", "echoindex", "crispressoalleles", "crispressonames"] # for these functions we don't want to alter the output file name
        listOutName = None
        if not args.function in exclusions:
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
                        # Special string-based functions
                        if args.function in ['listrename', 'twofastaseqidcompare', 'twofastaseqidcompare_orthofinder', 'mergefasta']:
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
                if args.function == 'single2multi' or args.function == 'cullbelow' or args.function == 'cullabove' or args.function == 'chunk' or args.function == 'reversecomplement2multi' or args.function == 'echoindex':
                        try:
                                args.number = int(args.number)
                        except:
                                print('The specified number is not accepted as an integer. Check your input to make sure this is a plain number (e.g., 5 or 100, not 5.00 or 1e-10) and try again.')
                                quit()
        return listOutName, args.outputFileName

def main():
        '''
        To add a new function into this program, you need to 1) add a new description  
        to the validate_args function, 2) handle it specifically if it is a special string, integer 
        or float-based function, 3) add the actual function above, 4) add it into the function list,
        and 5) enact the function below.
        '''
        # Fix Bio.SeqIO fastq parsing function
        SeqIO.QualityIO.FastqGeneralIterator = AltFastqGeneralIterator # This helps in cases where qual IDs differ from title IDs, preventing a program-terminating error
        
        # Function list - update as new ones are added
        stringFunctions = ['rename', 'listrename', 'appendrename', 'removeseqwstring', 'removeseqidwstring', 'retrieveseqwstring', 'retrieveseqidwstring', 'removestringfseqid', 'splitseqidatstring_start', 'splitseqidatstring_end', 'trim', 'twofastaseqidcompare', 'twofastaseqidcompare_orthofinder', 'mergefasta']
        numberFunctions = ['single2multi', 'cullbelow', 'cullabove', 'chunk', 'reversecomplement2multi', 'echoindex']
        basicFunctions = ['ids', 'descriptions', 'lengths', 'lengths_tsv', 'count', 'multi2single', 'q_to_a', 'reversecomplement', 'striphyphens', 'gc', 'explodeintocontigs', 'crispressoalleles', 'crispressonames']
        functionList = stringFunctions + numberFunctions + basicFunctions
        
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
                       help="Output file name (in some cases this should just be a prefix)")
        p.add_argument("-H", "-HELP", dest="detailedHelp", action='store_true',
                     help="Provide detailed help for each function")
        
        args = p.parse_args()
        listOutName, args.outputFileName = validate_args(args, stringFunctions, numberFunctions, functionList)
        
        # Enact functions
        ## String functions
        if args.function == 'rename':
                rename(args.fastaFileName, args.string, args.outputFileName, listOutName)
        if args.function == 'listrename':
                listrename(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'appendrename':
                appendrename(args.fastaFileName, args.string, args.outputFileName, listOutName)
        if args.function == 'removestringfseqid':
                removestringfseqid(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'splitseqidatstring_start':
                splitseqidatstring_start(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'splitseqidatstring_end':
                splitseqidatstring_end(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'retrieveseqwstring':
                retrieveseqwstring(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'retrieveseqidwstring':
                retrieveseqidwstring(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'removeseqwstring':
                removeseqwstring(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'removeseqidwstring':
                removeseqidwstring(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'trim':
                trim(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'twofastaseqidcompare':
                twofastaseqidcompare(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'twofastaseqidcompare_orthofinder':
                twofastaseqidcompare_orthofinder(args.fastaFileName, args.string, args.outputFileName)
        if args.function == 'mergefasta':
                mergefasta(args.fastaFileName, args.string, args.outputFileName)
        ## Number functions
        if args.function == 'single2multi':
                single2multi(args.fastaFileName, args.number, args.outputFileName)
        if args.function == 'cullbelow':
                cullbelow(args.fastaFileName, args.number, args.outputFileName)
        if args.function == 'cullabove':
                cullabove(args.fastaFileName, args.number, args.outputFileName)
        if args.function == 'reversecomplement2multi':
                reversecomplement2multi(args.fastaFileName, args.number, args.outputFileName)
        ## Number functions - FAST(A/Q) compatible
        if args.function == 'chunk':
                chunk(args.fastaFileName, args.number, args.outputFileName)
        if args.function == 'echoindex':
                echoindex(args.fastaFileName, args.number)
        ## Basic functions
        if args.function == 'ids':
                ids(args.fastaFileName, args.outputFileName)
        if args.function == 'descriptions':
                descriptions(args.fastaFileName, args.outputFileName)
        if args.function == 'lengths':
                lengths(args.fastaFileName, args.outputFileName)
        if args.function == 'lengths_tsv':
                lengths_tsv(args.fastaFileName, args.outputFileName)
        if args.function == 'count':
                count(args.fastaFileName, args.outputFileName)
        if args.function == 'multi2single':
                multi2single(args.fastaFileName, args.outputFileName)
        if args.function == 'q_to_a':
                q_to_a(args.fastaFileName, args.outputFileName)
        if args.function == 'reversecomplement':
                reversecomplement(args.fastaFileName, args.outputFileName)
        if args.function == 'striphyphens':
                striphyphens(args.fastaFileName, args.outputFileName)
        if args.function == 'gc':
                gc(args.fastaFileName, args.outputFileName)
        if args.function == 'explodeintocontigs':
                explodeintocontigs(args.fastaFileName, args.outputFileName)
        if args.function == 'crispressoalleles':
                crispressoalleles(args.fastaFileName)
        if args.function == 'crispressonames':
                crispressonames(args.fastaFileName)
        print('Program completed successfully!')

if __name__ == '__main__':
        main()
