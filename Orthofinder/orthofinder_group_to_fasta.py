#! python3
# orthofinder_group_to_fasta
# Script to extract individual .fasta files which contain the sequences
# which make up an orthogroup.

import os, argparse
from pyfaidx import Fasta

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Setup
        import platform
        # Ensure all arguments are specified
        for key, value in vars(args).items():
                # Special handling for MAFFT-related arguments
                if args.align:
                        if value == None or value == [] and key != 'mafftdir':  # mafftdir is allowed to be None
                                print(key + ' argument was not specified; fix your input and try again.')
                                quit()
                else:
                        if key in ['mafftdir', 'threads']:
                                continue
                        elif value == None or value == []:
                                print(key + ' argument was not specified; fix your input and try again.')
                                quit()
        # Validate that MAFFT is locatable if relevant
        if args.align:
                if platform.system() == 'Windows':
                        program_execution_check(os.path.join(args.mafftdir, 'mafft.bat -h'))
                else:
                        program_execution_check(os.path.join(args.mafftdir, 'mafft -h'))
        # Validate Orthogroup file location
        if not os.path.isfile(args.orthogroups):
                print('I am unable to locate the Orthogroups file (' + args.orthogroup + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate the FASTA file input location depending on type of input
        if len(args.inputLocation) == 1:
                # Check that specified value is a path
                if not os.path.isdir(args.inputLocation[0]):
                        print('One value was provided for -i, which means you should have provided a directory containing FASTA files.')
                        print('The provided value "' + args.inputLocation[0] + '" is not a directory; either it does not exist or it is a file; fix your input and try again.')
                        quit()
        else:
                # Check that the specified values are files
                for file in args.inputLocation:
                        if not os.path.isfile(file):
                                print('Multiple values were provided for -i, which means you should have provided the location of individual FASTA files.')
                                print('The provided value "' + file + '" is not a file; either it does not exist or it is a directory; fix your input and try again.')
                                quit()
        # Ensure that the output location is sensible
        if os.path.isfile(args.outputLocation):
                print('The specified output location "' + args.outputLocation + '" is a file. You should be specifying a directory (that may or may not exist); fix your input and try again.')
                quit()
        elif not os.path.isdir(args.outputLocation):
                pathSplit = os.path.split(args.outputLocation)
                if not os.path.isdir(pathSplit[0]):
                        print('The specified output location "' + args.outputLocation + '" is not a directory.')
                        print('The location one directory up i.e., "' + pathSplit[0] + '" is also not a directory.')
                        print('This program won\'t create new directories that aren\'t contained within an existing one.')
                        print('Either specify a new location or create one of the above two mentioned locations and try again.')
                        quit()

def program_execution_check(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '' and not 'cannot open -h' in cmderr.decode("utf-8").lower():     # Need this extra check with MAFFT since it can't run without a file input and will pop up an error to our -h; this still tells us it exists
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

## OrthoFinder-related
def parse_orthogroups_csv(orthogroupsFile):
        # Setup
        header = None
        orthoDict = {}
        # Main function
        with open(orthogroupsFile , 'r') as orthoFile:
                for line in orthoFile:
                        if header == None:
                                header = line.rstrip('\r\n').strip('"').split('\t')[1:]                 # [1:] gets rid of the blank space at the start of the file
                                # Establish subdicts for later indexing
                                for i in range(len(header)):
                                        orthoDict[header[i].strip('"')] = {}
                        else:
                                sl = line.rstrip('\r\n').strip('"').replace('""', '"').split('\t')      # For some reason sequence IDs with " in them have it replaced with ""
                                for i in range(1, len(sl)):
                                        sl[i] = sl[i].strip('"').split(', ')    # OrthoFinder separates sequence IDs like so; this works out since OrthoFinder removes commas from the IDs (which we need to handle during the parse)
                                        if sl[i] == ['']:
                                                sl[i] = []
                                        # Index species by orthogroup
                                        if i == 1:
                                                orthoDict[sl[0]] = {}
                                        orthoDict[sl[0]][header[i-1]] = sl[i]   # i-1 since sl has the orthogroup ID as its first value
                                # Index orthogroup by species
                                for i in range(len(header)):
                                        orthoDict[header[i]][sl[0]] = sl[i+1]   # i+1 since the first value in sl is the orthogroup ID
        return orthoDict, header

## Pyfaidx-related functions
def pyfaidx_long_indexing_orthofinder(pyfaidxRecords):
        # Set up
        longIndex = {}
        # Main function
        for record in pyfaidxRecords:
                orthoName = record.long_name.replace('(', '_').replace(')', '_').replace(',', '_')      # OrthoFinder replaces problem characters with underscores - bit annoying but we can handle it
                longIndex[orthoName] = record.name
        return longIndex

## MAFFT alignment-related
def mafft_align_filelist(mafftdir, fastaFileList, outputDir, outputFileNameList, threads, algorithm):   # fastaFileList and outputFileNameList are expected to be paired
        # Set up
        import os, threading, math
        from Bio.Align.Applications import MafftCommandline
        # Ensure that algorithm value is sensible
        if algorithm != None:   # If this is None, we'll just use default MAFFT
                if algorithm.lower() not in ['genafpair', 'localpair', 'globalpair']:
                        print('mafft_align: algorithm option must be an option in the below list. Fix this parameter and try again.')
                        print(['genafpair', 'localpair', 'globalpair'])
                        quit()
        # Ensure that the paired lists are equivalent
        if len(fastaFileList) != len(outputFileNameList):
                print('mafft_align_filelist: Length of the input FASTA list and output name list are not equal')
                print('fastaFileList length == ' + str(len(fastaFileList)) + '; outputFileNameList length == ' + str(len(outputFileNameList)))
                print('Fix the code for this program.')
                quit()
        # Ensure that mafftdir is usable
        if mafftdir == None:
                mafftdir = ''   # It's okay for mafftdir to not be specified since, in this case, we assume it's available from the PATH
        # Define functions integral to this one
        def run_mafft(mafftdir, outputDir, fastaFileList, outputFileNameList, startNum, endNum, algorithm):
                # Set up
                import platform
                for i in range(startNum, endNum):
                        # Run MAFFT
                        if platform.system() == 'Windows':
                                mafft_cline = MafftCommandline(os.path.join(mafftdir, 'mafft.bat'), input=fastaFileList[i])
                        else:
                                mafft_cline = MafftCommandline(os.path.join(mafftdir, 'mafft'), input=fastaFileList[i])
                        if algorithm != None:
                                if algorithm.lower() == 'genafpair':
                                        mafft_cline.genafpair = True
                                elif algorithm.lower() == 'localpair':
                                        mafft_cline.localpair = True
                                elif algorithm.lower() == 'globalpair':
                                        mafft_cline.globalpair = True
                        stdout, stderr = mafft_cline()
                        if stdout == '':
                                raise Exception('MAFFT error text below' + str(stderr))
                        # Process MAFFT output
                        stdout = stdout.split('\n')
                        while stdout[-1] == '\n' or stdout[-1] == '' or stdout[-1] == 'Terminate batch job (Y/N)?\n':   # Remove junk, sometimes MAFFT will have the 'Terminate ...' line
                                del stdout[-1]
                        stdout = '\n'.join(stdout)
                        # Create output alignment files
                        with open(outputFileNameList[i], 'w') as fileOut:
                                fileOut.write(stdout)

        # Set up threading requirements                         # This threading system is derived from chunk_fasta in (what is currently called) domfind.py
        fileListSize = len(fastaFileList)
        rawNum = fileListSize / threads                         # In cases where threads > fileListSize, rawNum will be less than 1. numRoundedUp will equal the number of threads, and so we'll end up rounding these to 1. Yay!
        numRoundedUp = round((rawNum % 1) * threads, 0)         # By taking the decimal place and multiplying it by the num of threads, we can figure out how many threads need to be rounded up to process every cluster
        chunkPoints = []
        ongoingCount = 0
        for i in range(int(threads)):
                if i+1 <= numRoundedUp:                         # i.e., if two threads are being rounded up, we'll round up the first two loops of this
                        chunkPoints.append([ongoingCount, math.ceil(rawNum) + ongoingCount])    # Round up the rawNum, and also add our ongoingCount which corresponds to the number of clusters already put into a chunk
                        ongoingCount += math.ceil(rawNum)                                       # Unlike chunk_fasta, we're storing a paired value of ongoingCount and the chunk point
                else:                                                                           # Our mafft function iterates over a range, so we go up to and not including the last value; this system is compliant with that style of sorting
                        chunkPoints.append([ongoingCount, math.floor(rawNum) + ongoingCount])   # Also note that group_dict is indexed starting from 0, so if group_dict len == 10, we want to iterate over range(0,10) since the last actual index is 9
                        ongoingCount += math.floor(rawNum)
                if ongoingCount >= fileListSize:                # Without this check, if we have more threads than clusters, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                        break
        # Begin the loop
        processing_threads = []
        ongoingCount = 0
        for start, end in chunkPoints:
                build = threading.Thread(target=run_mafft, args=(mafftdir, outputDir, fastaFileList, outputFileNameList, start, end, algorithm))
                processing_threads.append(build)
                build.start()
                ongoingCount += 1

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()

#### USER INPUT SECTION
usage = """%(prog)s will parse the Orthogroups.csv file output by Orthofinder and 
extract and group sequences by orthogroup as individual FASTA files with names 
corresponding to orthogroup IDs e.g., "OG0000000.fasta". File input
can be managed in two ways; 1) A directory path can be provided which contains files
whose names (minus extension) match the column names of the .csv file. Use '.' for
current directory. 2) Files can be listed on the command-line. These files must be
provided in the same order as their corresponding columns are presented in the .csv
file from left to right and there must be the same number of provided files as
there is columns.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="inputLocation", nargs="+",
               help="Input a single directory or multiple file locations")
p.add_argument("-c", "-csv", dest="orthogroups",
               help="Specify the orthogroup file")
p.add_argument("-o", "-output", dest="outputLocation",
               help="Specify the location to write output files to; this location will be populated with orthogroup files")
p.add_argument("-a", dest="align", action = "store_true", default=False,
               help="Provide this tag and below listed arguments if you wish for orthogroups to be aligned.")
p.add_argument("-m", "-mafftdir", dest="mafftdir", type=str, default='',
               help="""Optionally specify the directory where the MAFFT executables are located
               on a Unix system or .bat file on a Windows system is located. If this is already
               in your PATH or you don\'t want to align orthogroups, you can leave this blank.""")
p.add_argument("-t", "-threads", dest="threads", type=int, default=1,
               help="Optionally specify the number of threads to run MAFFT alignment with (default == 1).")

args = p.parse_args()
validate_args(args)

# Parse orthogroup file
orthoDict, header = parse_orthogroups_csv(args.orthogroups)

# Find FASTA files depending on how inputLocation was specified
fastaFiles = []
headerCount = {}
for head in header:
        headerCount[head] = 0
if len(args.inputLocation) == 1:
        # Scan through files and detect our files of interest
        for file in os.listdir(args.inputLocation[0]):
                filePrefix = file.rsplit('.', maxsplit=1)[0]
                if filePrefix in header:
                        fastaFiles.append(os.path.join(args.inputLocation[0], file))
                        # Count this file against its respective header value
                        for head in header:
                                if filePrefix == head:
                                        headerCount[head] += 1
        # Ensure that we found all the files we should have
        if len(fastaFiles) < len(header):
                print('I did not find all the expected files in directory "' + args.inputLocation[0] + '"')
                print('Based on "' + args.orthogroups + '" column headers, I should have found files starting with these values: ' + str(header))
                print('I did not find values matching to these headers:')
                for key, value in headerCount.items():
                        if value == 0:
                                print(key)
                print('Make sure you specified the correct location, or make sure all files are located at this directory. Program will exit now.')
                quit()
        elif len(fastaFiles) > len(header):
                print('I found too many files in directory "' + args.inputLocation[0] + '"')
                print('Based on column headers, I should have found 1 file starting with each of these values: ' + str(header))
                print('I found multiple values matching to these headers:')
                for key, value in headerCount.items():
                        if value > 1:
                                print(key)
                print('Make sure you specified the correct location, or make sure there is only 1 file with the same name as each column header. Program will exit now.')
                quit()
else:
        # Make sure the number of provided arguments is sensible
        if len(args.inputLocation) != len(header):
                print('The amount of files specified on command line (n = ' + str(len(args.inputLocation)) + ') does not match the number of column headers in "' + args.orthogroups + '" (n = ' + str(len(header)) + ')')
                print('Make sure your input orthogroups file matches your FASTA file list and try again.')
                quit()
        # Make sure that the provided files all exist
        for file in args.inputLocation:
                if not os.path.isfile(file):
                        print('Input file "' + file + '" either does not exist or is not a file.')
                        print('Make sure you spelled the file name/location correctly and try again.')
                        quit()
                fastaFiles.append(os.path.abspath(file))

# Load the fasta files and parse their contents
fastaRecords = []
fastaIndices = []
for fasta in fastaFiles:
        fastaRecords.append(Fasta(fasta))
        fastaIndices.append(pyfaidx_long_indexing_orthofinder(fastaRecords[-1]))

# Make the output directory if it doesn't exist
if not os.path.isdir(args.outputLocation):
        os.mkdir(args.outputLocation)

# Main output loop
fileNames = []  # This is for subsequent MAFFT alignment (if it occurs)
for key, value in orthoDict.items():
        if key in header:       # This lets us skip the alternative indexing system values in orthoDict
                continue
        # Extract orthogroups as FASTA
        fileName = os.path.join(args.outputLocation, key + '.fasta')
        with open(fileName, 'w') as ogOut:
                fileNames.append(fileName)
                for i in range(len(header)):
                        if orthoDict[key][header[i]] != []:
                                for seqid in orthoDict[key][header[i]]:
                                        if seqid in fastaIndices[i]:
                                                ogOut.write('>' + header[i].split('_')[0] + '_' + seqid + '\n' + str(fastaRecords[i][fastaIndices[i][seqid]]) + '\n')
                                        elif seqid in fastaRecords[i]:
                                                ogOut.write('>' + header[i].split('_')[0] + '_' + seqid + '\n' + str(fastaRecords[i][seqid]) + '\n')
                                        else:
                                                print('"' + seqid + '" was not found within column\'s corresponding FASTA file ("' + fastaFiles[i] + '")')
                                                print('If you specified files on the command-line, you must make sure these files are ordered to correspond exactly to ' + str(header))
                                                print('Otherwise, your FASTA files are either incorrect or this orthogroups CSV file does not match the files. You need to figure out what\'s going wrong here, sorry.')
                                                print('Problem sequence ID = "' + seqid + '"; note that this is how the ID is presented within your orthogroups file and it\'s possible OrthoFinder substituted characters')
                                                print('Program will exit now.')
                                                quit()

# Perform MAFFT alignment if relevant
if args.align:
        mafft_align_filelist(args.mafftdir, fileNames, args.outputLocation, fileNames, args.threads, 'genafpair')       # Note that we're overwriting the original files by doing this

# Done!
print('Program completed successfully!')
