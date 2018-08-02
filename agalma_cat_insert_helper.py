#! python3

import os, argparse

# Define functions for later use
def validate_args(args):
        # Validate input file location
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate program argument
        if args.agalmaDir == None:
                args.agalmaDir = ''
        program_execution(os.path.join(args.agalmaDir, 'agalma -h'))
        # Validate numeric arguments
        if args.threads < 1:
                print('Threads argument must be at least 1.')
                quit()
        elif args.threads == None:
                print('Threads argument must be specified.')
                quit()
        if args.memory < 1:
                print('Memory argument must be at least 1.')
                quit()
        elif args.memory == None:
                print('Memory argument must be specified.')
                quit()
        return args

def program_execution(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '':     # Need this extra check for seg since it puts its usage information into stderr rather than stdout
                # Check if this is actually failed
                for line in cmderr.decode("utf-8").split('\n'):
                        if line.startswith('biolite.config.parse_env_resources: database='):
                                continue
                else:
                        print('Failed to execute program "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                        print('---')
                        print('stderr is below for debugging purposes.')
                        print(cmderr.decode("utf-8"))
                        print('Program closing now.')
                        quit()

def agalma_info_table_parse(agalmaDir, fileName, threads, mem):
        # Set up
        outCmds = []
        # Main function
        with open(fileName, 'r') as agalmaFile:
                for line in agalmaFile:
                        if line.startswith('#'):
                                continue
                        sl = line.rstrip('\r\n').split('\t')
                        # Extract information
                        species = sl[3]
                        catID = species.replace(' ', '_')
                        # Format reads string
                        readsList = sl[1].split(' ')
                        fileExt = ''
                        if '.fq' or '.fastq' in readsList[0]:
                                fileExt += '.fastq'
                        if '.gz' or '.gzip' in readsList[0]:
                                fileExt += '.gz'
                        if len(readsList) > 1:
                                readsString = catID + '_1' + fileExt + ' ' + catID + '_2' + fileExt
                        else:
                                readsString = catID + '_SE' + fileExt
                        # Format catalog command
                        cmd = os.path.join(agalmaDir, 'agalma') + ' -t ' + str(threads) + ' -m ' + str(mem) + 'G' + ' catalog insert --id ' + catID + ' --paths ' + readsString + ' --species "' + sl[3] + '"'
                        if sl[4] != '.':
                                sl[4] = sl[4].strip(' ')
                                cmd += ' -n ' + sl[4]
                        if sl[5] != '.':
                                sl[5] = sl[5].strip(' ')
                                cmd += ' -d ' + sl[5]
                        outCmds.append(cmd)
        outCmds = list(set(outCmds))
        return outCmds

#### USER INPUT SECTION
usage = """%(prog)s will help with formatting AGALMA catalog inserts in a way that is consistent
and able to be repeated easily
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-inputTable", dest="inputTable",
                  help="Input tab-delimited information table file name.")
p.add_argument("-a", "-agalmaDir", dest="agalmaDir", type = str,
                  help="Specify the directory where AGALMA executables are located.")
p.add_argument("-t", "-threads", dest="threads", type = int,
                  help="Specify the number of threads to provide to AGALMA.")
p.add_argument("-m", "-memory", dest="memory", type = int,
                  help="Specify the amount of memory (in Gb) to provide to AGALMA.")

args = p.parse_args()
args = validate_args(args)

# Parse information table
cmds = agalma_info_table_parse(args.agalmaDir, args.inputTable, args.threads, args.memory)

# Run catalog commands
for cmd in cmds:
        print(cmd)
        program_execution(cmd)

# Done!
print('Program completed successfully!')
