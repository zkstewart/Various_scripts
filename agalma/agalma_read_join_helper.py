#! python3

import os, argparse

# Define functions for later use
def validate_args(args):
        # Validate input file location
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()

def program_execution(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '' and not cmderr.decode("utf-8").startswith('Usage'):     # Need this extra check for seg since it puts its usage information into stderr rather than stdout
                print('Failed to execute "' + cmd + '".')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

def agalma_info_table_parse_sp_read_dict(fileName):
        # Set up
        speciesReads = {}
        # Main function
        with open(fileName, 'r') as agalmaFile:
                for line in agalmaFile:
                        if line.startswith('#'):
                                continue
                        sl = line.rstrip('\r\n').split('\t')
                        # Check to see if we should be skipping this value
                        if sl[6].lower() == 'y':
                                continue
                        # Extract information
                        readsList = sl[1].split(' ')
                        for i in range(len(readsList)):
                                readsList[i] = [os.path.join(sl[0], readsList[i])]
                        species = sl[3]
                        # Hold onto species reads
                        if species not in speciesReads:
                                speciesReads[species] = readsList
                        else:
                                r1, r2 = readsList
                                speciesReads[species][0] += r1
                                speciesReads[species][1] += r2
        return speciesReads

def file_cat_cmds_by_dict(inputDict):
        # Set up 
        cmds = []
        # Main function
        for key, value in inputDict.items():
                fileExt = ''
                if '.fq' or '.fastq' in value[0]:
                        fileExt += '.fastq'
                if '.gz' or '.gzip' in value[0]:
                        fileExt += '.gz'
                if len(value) > 1:
                        r1Cmd = 'cat ' + ' '.join(value[0]) + ' > ' + key.replace(' ', '_') + '_1' + fileExt
                        r2Cmd = 'cat ' + ' '.join(value[1]) + ' > ' + key.replace(' ', '_') + '_2' + fileExt
                        cmds += [r1Cmd, r2Cmd]
                else:
                        r1Cmd = 'cat ' + ' '.join(value[0]) + ' > ' + key.replace(' ', '_') + '_SE' + fileExt
                        cmds += [r1Cmd]
        return cmds

#### USER INPUT SECTION
usage = """%(prog)s will help with formatting AGALMA catalog inserts in a way that is consistent
and able to be repeated easily. Specifically, this code will facilitate the concatenation
of multiple separate paired-end read sets into a single paired-end set (which is what AGALMA requires).
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-inputTable", dest="inputTable",
                  help="Input tab-delimited information table file name.")

args = p.parse_args()
validate_args(args)

# Parse information table
speciesReads = agalma_info_table_parse_sp_read_dict(args.inputTable)

# Generate cat commands
cmds = file_cat_cmds_by_dict(speciesReads)

# Run cat commands
for cmd in cmds:
        outFile = cmd.split('> ')[1]
        if not os.path.isfile(outFile):
                print(cmd)
                program_execution(cmd)

# Done!
print('Program completed successfully!')
