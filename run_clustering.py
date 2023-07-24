#! python3
# run_mmseqs2.py
# Wrapper script to make it easier to run a search with MMseqs2

import argparse, os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from Function_packages import ZS_ClustIO

# Define functions
def validate_args(args):
    # Validate required parameters
    if not os.path.isfile(args.fastaFile):
        print(f'I am unable to locate the input FASTA file ({args.fastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.threads < 1:
        print("threads must be a positive integer")
        quit()
    
    # Validate "BOTH" parameters
    if not 0 <= args.identity <= 1.0:
        print("identity must be a float in the range 0.0 -> 1.0")
        quit()
    
    # Validate "MMS" parameters
    if not os.path.isdir(args.mmseqsDir):
        print(f'I am unable to locate the MMseqs2 directory ({args.mmseqsDir})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    "--tmpDir is validated by the MM_DB Class"
    if args.evalue < 0:
        print("evalue must be greater than or equal to 0")
        quit()
    if not 0 <= args.coverage <= 1.0:
        print("coverage must be a float in the range 0.0 -> 1.0")
        quit()
    "--mode is controlled by argparse choices"
    
    # Validate "MMS-CASCADE" parameters
    "--sensitivity is controlled by argparse choices"
    if args.steps < 1:
        print("steps must be greater than or equal to 1")
        quit()
    
    # Validate "CD" parameters
    if not os.path.isdir(args.cdhitDir):
        print(f'I am unable to locate the CD-HIT directory ({args.cdhitDir})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    "--molecule is controlled by argparse choices"
    if args.mem < 1000:
        print("mem must be greater than or equal to 1000")
        quit()
    "--align is controlled by argparse choices"
    if not 0 <= args.shorterCov <= 1.0:
        print("shorterCov must be a float in the range 0.0 -> 1.0")
        quit()
    if not 0 <= args.longerCov <= 1.0:
        print("shorterCov must be a float in the range 0.0 -> 1.0")
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()
    if args.program == "cd-hit":
        if os.path.isfile(args.outputFileName + ".clstr"):
            print(f'File already exists at output location ({args.outputFileName + ".clstr"})')
            print('Make sure you specify a unique file name and try again.')
            quit()

# Log file related functions
def temp_file_name_gen(prefix):
    ongoingCount = 1
    while True:
        if not os.path.isfile(prefix):
            return prefix
        elif os.path.isfile(prefix + str(ongoingCount)):
            ongoingCount += 1
        else:
            return prefix + str(ongoingCount)

def log_update(logName, text):
    with open(logName, 'a') as logFile:
        logFile.write(text + '\n')

# Main call
def main():
    usage = """Wrapper script to perform sequence clustering using MMseqs2 or CD-HIT.
    Some parameters are used for all programs, and others are specific depending on which
    algorithm is being used; these are given prefixes in the help text below to help you
    figure out what you need to specify (or leave at the indicated defaults). Note that
    CD-HIT produces the representatives FASTA and cluster TSV file by default. MMseqs
    clustering at this point just generates a cluster TSV file for you to pick
    representatives from yourself."""
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="fastaFile",
                   required=True,
                   help="Input FASTA file for clustering")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for clustering results")
    p.add_argument("-p", dest="program",
                   required=True,
                   choices=["mmseqs-cascade", "mmseqs-linclust", "cd-hit"],
                   help="Specify which algorithm to use for clustering")
    p.add_argument("-t", dest="threads",
                   required=True,
                   type=int,
                   help="Specify how many threads to operate with")
    # Optional - both
    p.add_argument("--identity", dest="identity",
                   required=False,
                   type=float,
                   help="""BOTH: Specify the identity threshold for clustering;
                   default==0.9""")
    # Optional - MMseqs2
    p.add_argument("--mmseqs", dest="mmseqsDir",
                   required=False,
                   help="MMS: Specify the directory containing the MMseqs2 executable")
    p.add_argument("--tmpDir", dest="tmpDir",
                   required=False,
                   help="""MMS: Specify the tmpDir for MMseqs2 running; default='mms2_tmp'
                   in your current working directory""",
                   default="mms2_tmp")
    p.add_argument("--evalue", dest="evalue",
                   required=False,
                   type=float,
                   help="MMS: Specify the evalue threshold for clustering; default==1e-3",
                   default=1e-3)
    p.add_argument("--coverage", dest="coverage",
                   required=False,
                   type=float,
                   help="MMS: Specify the coverage ratio for clustering; default==0.8",
                   default=0.8)
    p.add_argument("--mode", dest="mode",
                   required=False,
                   choices=["set-cover", "connected_component", "greedy"],
                   help="MMS: Specify the clustering mode; default=='set-cover'",
                   default="set-cover")
    p.add_argument("--sensitivity", dest="coverage",
                   required=False,
                   choices=[1,2,3,4,5,5.7,6,7,7.5],
                   help="MMS-CASCADE: Specify the sensitivity value; default==4",
                   default=4)
    p.add_argument("--steps", dest="steps",
                   required=False,
                   type=int,
                   help="""MMS-CASCADE: Specify the number of cascaded clustering steps;
                   default==3""",
                   default=3)
    # Optional - CD-HIT
    p.add_argument("--cdhit", dest="cdhitDir",
                   required=False,
                   help="CD: Specify the directory containing the CD-HIT executables")
    p.add_argument("--molecule", dest="molecule",
                   required=False,
                   choices=["protein", "nucleotide"],
                   help="CD: Specify the molecule type contained in your FASTA file")
    p.add_argument("--mem", dest="mem",
                   required=False,
                   type=int,
                   help="CD: Specify the number of MB of memory to use; default==5000",
                   default=5000)
    p.add_argument("--align", dest="align",
                   required=False,
                   choices=["local", "global"],
                   help="""CD: Specify whether local or global alignment should be used;
                   default == 'global'""")
    p.add_argument("--shorter_coverage", dest="shorterCov",
                   required=False,
                   type=float,
                   help="CD: Specify the coverage ratio for the shorter seq; default==0.0",
                   default=0.0)
    p.add_argument("--longer_coverage", dest="longerCov",
                   required=False,
                   type=float,
                   help="CD: Specify the coverage ratio for the longer seq; default==0.0",
                   default=0.0)
    args = p.parse_args()
    validate_args(args)
    
    # Create log file for ongoing updates of program status
    logName = temp_file_name_gen(os.path.join(args.outputdir, args.output + '.mms2log'))
    open(logName, 'w').close()
    print("# Note that a current quirk of this software means that log file updates " +
          "are made AFTER the command actually finished executing. I might try " + 
          "to fix this at a later date.")
    
    # Branch logic depending on which cluster algorithm is being used
    if args.program == "mmseqs-cascade" or args.program == "mmseqs-linclust":
        main_mmseqs(args, logName)
    elif args.program == "cd-hit":
        main_cdhit(args, logName)
    
    print('Program completed successfully!')

def main_mmseqs(args, logName):
    '''
    Functionally part of the main namespace. Handles the running of MMseqs2 clustering
    with linclust or cascaded clustering.
    '''
    mmDB = ZS_ClustIO.MM_DB(args.fastaFile, args.mmseqsDir, args.tmpDir, args.threads)
    
    # Generate the sequence database
    logText = mmDB.generate()
    log_update(logName, logText)
    
    # Index the sequence database
    logText = mmDB.index()
    log_update(logName, logText)
    
    # Run clustering
    if args.program == "mmseqs-cascade":
        clusterer = ZS_ClustIO.MM_Cascade(
            mmDB, args.evalue, args.identity, args.coverage,
            args.mode, args.threads, args.tmpDir,
            args.sensitivity, args.steps
        )
    else:
        clusterer = ZS_ClustIO.MM_Linclust(
            mmDB, args.evalue, args.identity, args.coverage,
            args.mode, args.threads, args.tmpDir
        )
    
    logText = clusterer.cluster()
    log_update(logName, logText)
    
    # Generate output file
    logText = clusterer.tabulate(args.outputFileName)
    log_update(logName, logText)

def main_cdhit(args, logName):
    '''
    Functionally part of the main namespace. Handles the running of CD-HIT clustering.
    '''
    
    # Parameterise the clusterer object
    cdhit = ZS_ClustIO.CDHIT(args.fastaFile, args.molecule, args.cdhitDir)
    cdhit.identity = args.identity
    cdhit.mem = args.mem
    cdhit.threads = args.threads
    cdhit.clean = False
    cdhit.set_shorter_cov_pct(args.shorterCov)
    cdhit.set_longer_cov_pct(args.longerCov)
    if args.align == "local":
        cdhit.set_local()
    
    # Run CD-HIT
    outputFileName = os.path.abspath(args.outputFileName)
    logText = cdhit.cdhit(args.fastaFile, os.path.dirname(outputFileName), os.path.basename(outputFileName))
    log_update(logName, logText)
    
if __name__ == '__main__':
    main()
