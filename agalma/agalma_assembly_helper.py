#! python3

import os, argparse

# Define functions for later use
def validate_args(args):
        # Validate input file location
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate program argument
        if args.agalmaDir == None:
                args.agalmaDir = ''
        program_execution(os.path.join(args.agalmaDir, 'agalma -h'))
        # Validate numeric arguments
        if args.threads == None:
                print('Threads argument must be specified.')
                quit()
        elif args.threads < 1:
                print('Threads argument must be at least 1.')
                quit()
        if args.memory == None:
                print('Memory argument must be specified.')
                quit()
        elif args.memory < 1:
                print('Memory argument must be at least 1.')
                quit()
        return args

def program_execution(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '':
                print('Program returned stderr; check for any possible debugging purposes.')
                print(cmderr.decode("utf-8"))

def agalma_info_table_parse(agalmaDir, fileName, threads, mem):
        # Set up
        outCmds = []
        # Main function
        with open(fileName, 'r') as agalmaFile:
                for line in agalmaFile:
                        if line.startswith('#'):
                                continue
                        sl = line.rstrip('\r\n').split('\t')
                        # Check to see if we should be skipping this value
                        if sl[7].lower() == 'y':
                                continue
                        # Extract information
                        species = sl[3]
                        catID = species.replace(' ', '_')
                        # Format assembly command
                        cmd = os.path.join(agalmaDir, 'agalma') + ' -t ' + str(threads) + ' -m ' + str(mem) + 'G' + ' transcriptome --id ' + catID
                        if cmd not in outCmds:
                                outCmds.append(cmd)
        return outCmds

#### USER INPUT SECTION
usage = """%(prog)s will help with running AGALMA transcriptome assembly commands.
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
