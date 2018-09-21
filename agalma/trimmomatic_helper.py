#! python3
# trimmomatic_helper.py
# This code will help to automate large Trimmomatic jobs involving numerous
# reads files. By formatting a table with 4 columns [Directory\tFiles\tPrefix\Complete]
# we can help with this process. Files should be space-separated if the reads are paired
# otherwise there should only be a single value here. Directory and files will be
# concatenated to render the file locations; prefix will be provided as -baseout to
# Trimmomatic and we'll keep a count of the number of times this prefix occurs to add
# to this value (e.g., if "prefix" occurs two times, the first occurrence's -baseout
# will be "prefix1", the next will be "prefix2", etc.). Complete will simply check
# to see if a 'Y' is present; if so, we'll skip this (but still count it for the prefix
# occurrence value).

import os, argparse

##### FUNCTION DEFINITION

# Validate arguments
def validate_args(args):
        # Validate input file location
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate that the threads argument is sensible
        if args.numThreads < 1:
                print('-n numThreads argument cannot be less than 1. Fix your input and try again.')
                quit()
        # Validate Trimmomatic file location
        if not os.path.isfile(args.trimmomaticJar):
                print('I am unable to locate the Trimmomatic .jar file (' + args.trimmomaticJar + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Format the trimmomatic commands
        if '{}' in args.trimSECommand:
                formatCount = args.trimSECommand.count('{}')
                trimmoDir = [os.path.dirname(args.trimmomaticJar)]*formatCount  # This lets us handle multiple occurrences of {}; this shouldn't happen but maybe the user has a reason for wanting Trimmomatic's dir to be called more than once
                args.trimSECommand = args.trimSECommand.format(*trimmoDir)
        if '{}' in args.trimPECommand:
                formatCount = args.trimPECommand.count('{}')
                trimmoDir = [os.path.dirname(args.trimmomaticJar)]*formatCount  # This lets us handle multiple occurrences of {}; this shouldn't happen but maybe the user has a reason for wanting Trimmomatic's dir to be called more than once
                args.trimPECommand = args.trimPECommand.format(*trimmoDir)
        # Validate that the adapters file is locatable
        SEAdapterFile = args.trimSECommand.split(':')[1]                          # I think this should always correspond to the adapters file; if it doesn't I'll give informative error text to the user
        if not os.path.isfile(SEAdapterFile):
                print('I am unable to locate the SE adapter file used for Trimmomatic. From analysis of the Trimmomatic command, I think it should be here "' + SEAdapterFile + '".')
                print('Does this seem right to you? If not, the problem likely arises from my expectation that your adapter file is located within the first two instances of the ":" character.')
                print('e.g., for "ILLUMINACLIP:{}/adapters/TruSeq3-SE.fa:2:30:10", {} will be substituted with "/home/user/Trimmomatic" (for example); then, I expect that this file is between "ILLUMINACLIP:" and ":2:30:10"')
                print('Try to fix your input to reflect this format, or just use the default argument.')
                quit()
        PEAdapterFile = args.trimPECommand.split(':')[1]                          # I think this should always correspond to the adapters file; if it doesn't I'll give informative error text to the user
        if not os.path.isfile(PEAdapterFile):
                print('I am unable to locate the PE adapter file used for Trimmomatic. From analysis of the Trimmomatic command, I think it should be here "' + PEAdapterFile + '".')
                print('Does this seem right to you? If not, the problem likely arises from my expectation that your adapter file is located within the first two instances of the ":" character.')
                print('e.g., for "ILLUMINACLIP:{}/adapters/TruSeq3-SE.fa:2:30:10", {} will be substituted with "/home/user/Trimmomatic" (for example); then, I expect that this file is between "ILLUMINACLIP:" and ":2:30:10"')
                print('Try to fix your input to reflect this format, or just use the default argument.')
                quit()
        # Validate that Java runs
        program_execution(os.path.join(args.javadir, 'java -h'))
        return args

# Execute program commands
def program_execution(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '' and not cmderr.decode("utf-8").startswith('Usage') and not 'completed successfully' in cmderr.decode("utf-8").lower():  # Trimmomatic prints 'Completed successfully' when it works
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

# Input table handling
def trimmomatic_reads_table_to_cmds(fileName, javaDir, trimmomaticJar, SECmd, PECmd, threads):
        # Set up
        import os
        prefixes = {}
        cmds = []
        # Main function
        with open(fileName, 'r') as trimmoFile:
                for line in trimmoFile:
                        # Skip comment/blank lines
                        if line.startswith('#') or line == '' or line == '\n' or line == '\r\n':
                                continue
                        sl = line.rstrip('\r\n').split('\t')
                        if set(sl) == {''}:
                                continue
                        # Extract information
                        readsList = sl[1].split(' ')
                        for i in range(len(readsList)):
                                readsList[i] = os.path.join(sl[0], readsList[i])
                        prefix = sl[2]
                        # Make sure our columns have information within them
                        if sl[0] == '' or sl[1] == '' or sl[2] == '':
                                print('trimmomatic_reads_table_to_cmds: There are empty values in the line "' + line)
                                print('These columns should all contain at least some sort of information. Fix this and try again.')
                                quit()
                        # Figure out if we're looking at PE or SE reads
                        if len(readsList) == 1:
                                readType = 'SE'
                        elif len(readsList) == 2:
                                readType = 'PE'
                        else:
                                print('trimmomatic_reads_table_to_cmds: I don\'t recognise the number of reads in this entry "' + sl[1] + '".')
                                print('There should be 1 or 2 values, not ' + str(len(readsList)) + '. Fix this value and try again.')
                                quit()
                        # Figure out which number we are looking at for this prefix
                        if prefix not in prefixes:
                                prefixes[prefix] = 1
                        else:
                                prefixes[prefix] += 1
                        prefixNum = prefixes[prefix]    # If this is the first occurrence of this prefix, we'll give it a suffix of '1'; otherwise this counts upwards
                        # Check to see if we should be skipping this value
                        if sl[3].lower() == 'y':
                                continue
                        # Make sure our reads files can be located
                        for readsFile in readsList:
                                if not os.path.isfile(readsFile):
                                        print(line)
                                        print(sl)
                                        print('trimmomatic_reads_table_to_cmds: I am unable to locate the reads file "' + readsFile + '".')
                                        print('This was formatted by concatenating column 1 and 2 of a line in your input file. Make sure these values are correct and try again.')
                                        quit()
                        # Format our output prefix and cmd, then store it
                        if readType == 'SE':
                                prefix += '_SE' + str(prefixNum) + '.trimmed.fq.gz'
                                cmd = '{} -jar {} {} -threads {} {} {} {}'.format(*[os.path.join(javaDir, 'java'), trimmomaticJar, readType, threads, readsList[0], prefix, SECmd])
                        else:
                                prefix += '_PE' + str(prefixNum) + '.trimmed.fq.gz'
                                cmd = '{} -jar {} {} -threads {} {} {} -baseout {} {}'.format(*[os.path.join(javaDir, 'java'), trimmomaticJar, readType, threads, readsList[0], readsList[1], prefix, PECmd])
                        cmds.append(cmd)
        return cmds

##### USER INPUT SECTION
usage = """%(prog)s will help with running Trimmomatic on a large number of read sets.
This script assumes java is locatable within your PATH.
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-inputTable", dest="inputTable",
                  help="Input tab-delimited information table file name.")
p.add_argument("-j", dest="javadir", type = str, default = '',
                  help="""Specify the directory where the java executable file is located.
                  If this is already in your PATH, you can leave this blank.""")
p.add_argument("-t", dest="trimmomaticJar", type = str,
                  help="Specify the location of the Trimmomatic .jar file.")
p.add_argument("-se", dest="trimSECommand", type = str, default = "ILLUMINACLIP:{}/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25",
                  help="""Optionally change the command to be used when running Trimmomatic (SE).
                  By default we will read in the fasta file within the adapters directory
                  (assumed to be in the same location as the trimmomatic .jar file; {} can
                  be used as a placeholder for this) and use Trinity's recommended settings.
                  Make sure any changes to this are correct.""")
p.add_argument("-pe", dest="trimPECommand", type = str, default = "ILLUMINACLIP:{}/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25",
                  help="""Optionally change the command to be used when running Trimmomatic (PE).
                  By default we will read in the fasta file within the adapters directory
                  (assumed to be in the same location as the trimmomatic .jar file; {} can
                  be used as a placeholder for this) and use Trinity's recommended settings.
                  Make sure any changes to this are correct.""")
p.add_argument("-n", dest="numThreads", type = int, default = 1,
                  help="Specify the number of threads to use when running Trimmomatic.")

args = p.parse_args()
args = validate_args(args)

##### MAIN LOOP

# Parse information table
cmds = trimmomatic_reads_table_to_cmds(args.inputTable, args.javadir, args.trimmomaticJar, args.trimSECommand, args.trimPECommand, args.numThreads)

# Run Trimmomatic commands
for cmd in cmds:
        print(cmd)
        program_execution(cmd)

# Done!
print('Program completed successfully!')
