#! python3

import os, argparse, re

# Define functions for later use
## Argument validation
def validate_args(args):
    # Validate input location
    if not os.path.isdir(args.readsDir):
        print('I am unable to locate the input FASTQs directory (' + args.readsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
        quit()

def get_rnaseq_files(readsDir, readsSuffix, isSingleEnd):
    # Locate files from the directory
    forwardReads = []
    reverseReads = []
    for file in os.listdir(readsDir):
        if file.endswith(readsSuffix):
            if isSingleEnd:
                forwardReads.append(os.path.join(readsDir, file))
            else:
                if file.endswith(f"1{readsSuffix}"):
                    forwardReads.append(os.path.join(readsDir, file))
                elif file.endswith(f"2{readsSuffix}"):
                    reverseReads.append(os.path.join(readsDir, file))
                else:
                    raise ValueError(f"{file} ends with the expected suffix '{readsSuffix}' but is not preceeded by a 1 or 2!")
    forwardReads.sort()
    reverseReads.sort()
    
    # Validate that paired files match
    if not isSingleEnd:
        assert len(forwardReads) == len(reverseReads), \
            f"Number of reads don't match for forward ({len(forwardReads)}) and reverse ({len(reverseReads)}) files"
        for i in range(len(forwardReads)):
            prefix = os.path.commonprefix([forwardReads[i], reverseReads[i]])
            assert prefix != "", \
                "forward and reverse read pairs don't have a common prefix?"
            assert forwardReads[i].startswith(f"{prefix}1") and reverseReads[i].startswith(f"{prefix}2"), \
                f"forward and reverse reads don't start with a recognised prefix ({prefix} should preceed a 1 or 2)"
    
    # Return files
    return forwardReads, reverseReads if reverseReads != [] else None

def group_reads_across_lanes(forwardReads, reverseReads, laneRegex):
    '''
    Parameters:
        forwardReads -- a list containing strings pointing to the absolute file paths
                        of forward reads from read pairs
        reverseReads -- a list containing strings pointing to the absolute file paths
                        of reverse reads from read pairs; can be None if handling
                        single end reads
        laneRegex -- a regex that can locate the lane identifier from a string
    Returns:
        forwardReadGroups -- a list containing sublists with structure like:
                             [
                                 [lane1_reads, lane2_reads],
                                 [lane1_reads],
                                 [lane1_reads, lane2_reads],
                                 ...
                             ]
        reverseReadGroups -- likewise with forwardReadGroups, but with the reverse
                             reads; can be None if handling single end reads
    '''
    def _group_reads(readsList, laneRegex):
        # Handle forward reads from pairs or single reads
        readGroups = []
        for readFile in readsList:
            readFileBase = os.path.basename(readFile)
            assert laneRegex.search(readFileBase) != None, \
                f"Lane identifier not found in '{readFileBase}'; can't handle this"
            
            readFileNoLane = laneRegex.sub("", readFileBase)
            
            foundGroup = False
            for readGroup in readGroups:
                groupIdentifier, group = readGroup
                if readFileNoLane == groupIdentifier:
                    group.append(readFile)
                    foundGroup = True
                    break
            if foundGroup == False:
                readGroups.append([readFileNoLane, [readFile]])
        return [ group for groupIdentifier, group in readGroups ]
    
    # Get forward reads group
    forwardReadGroups = _group_reads(forwardReads, laneRegex)
    
    # Handle reverse reads from pairs(if applicable)
    if reverseReads != None:
        reverseReadGroups = _group_reads(reverseReads, laneRegex)
    else:
        reverseReadGroups = None
    
    return forwardReadGroups, reverseReadGroups

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s receives a directory containing FASTQs that been split
    across several lanes. This script will produce a shell script that can be
    qsubbed to concatenate the FASTQs into a single forward/reverse file
    per sample.
    
    Some notes: laneIdentifier by default is "_L00", which means we can expect
    files from different lanes to have values like "_L001" and "_L002" for example.
    Otherwise, these file names should have no differences.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="readsDir",
                   required=True,
                   help="Input directory containing trimmed FASTQs")
    p.add_argument("-rs", dest="readsSuffix",
                   required=True,
                   help="""Suffix which uniquely identifies all relevant read files
                   e.g., 'P.fq.gz' for trimmomatic reads""")
    # Opts
    p.add_argument("-o", dest="outputFileName",
                   required=False,
                   help="Optionally specify shell script name (default==run_read_prep.sh)",
                   default="run_read_prep.sh")
    p.add_argument("--singleEnd", dest="isSingleEnd",
                   required=False,
                   action="store_true",
                   help="Optionally indicate whether the reads are expected to be single-ended rather than paired",
                   default=False)
    p.add_argument("--laneIdentifier", dest="laneIdentifier",
                   required=False,
                   help="""Optionally specify how we tell lanes apart (default==\"_L00\");
                   it is assumed this value will be followed by a single digit that identifies
                   the lane""",
                   default="_L00")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate files from the directory
    forwardReads, reverseReads = get_rnaseq_files(args.readsDir, args.readsSuffix, args.isSingleEnd)
    
    # Group reads from different lanes
    laneRegex = re.compile(args.laneIdentifier + r"\d")
    forwardReadGroups, reverseReadGroups = group_reads_across_lanes(forwardReads, reverseReads, laneRegex)
    
    # Format cat commands to join lane files
    catCmds = []
    for i in range(len(forwardReadGroups)):
        # Format forward read commands
        forwardGroup = forwardReadGroups[i]
        
        if len(forwardGroup) > 1:
            prefix = os.path.commonprefix([
                laneRegex.split(os.path.basename(f))[0].rstrip("_")
                    for f in forwardGroup
            ])
            
            cmd1 = "cat {forwardFiles} > {prefix}{suffix}.fq.gz".format(
                forwardFiles=" ".join(
                    [f"${{TRIMMEDDIR}}/{os.path.basename(f)}" for f in forwardGroup]
                ),
                prefix=prefix,
                suffix="_1" if not args.isSingleEnd else ""
            )
            catCmds.append(cmd1)
        elif len(forwardGroup) == 1:
            prefix = laneRegex.split(os.path.basename(forwardGroup[0]))[0].rstrip("_")
            
            cmd1 = "ln -s {forwardFiles[0]} {prefix}_1.fq.gz".format(
                forwardFiles=[f"${{TRIMMEDDIR}}/{os.path.basename(f)}" for f in forwardGroup],
                prefix=prefix,
                suffix="_1" if not args.isSingleEnd else ""
            )
            catCmds.append(cmd1)
        else:
            print(f"Error: no forward files (_1) found for a group; this shouldn't be possible?")
            quit()
        
        # Format reverse read commands
        reverseGroup = reverseReadGroups[i]
        
        if len(reverseGroup) > 1:
            prefix = os.path.commonprefix([
                laneRegex.split(os.path.basename(f))[0].rstrip("_")
                    for f in reverseGroup
            ])
            
            cmd2 = "cat {reverseFiles} > {prefix}{suffix}.fq.gz".format(
                reverseFiles=" ".join(
                    [f"${{TRIMMEDDIR}}/{os.path.basename(f)}" for f in reverseGroup]
                ),
                prefix=prefix,
                suffix="_2" if not args.isSingleEnd else ""
            )
            catCmds.append(cmd2)
        elif len(reverseGroup) == 1:
            prefix = laneRegex.split(os.path.basename(reverseGroup[0]))[0].rstrip("_")
            
            cmd2 = "ln -s {reverseFiles[0]} {prefix}{suffix}.fq.gz".format(
                reverseFiles=[f"${{TRIMMEDDIR}}/{os.path.basename(f)}" for f in reverseGroup],
                prefix=prefix,
                suffix="_2" if not args.isSingleEnd else ""
            )
            catCmds.append(cmd2)
        else:
            print(f"Error: no reverse files (_2) found for a group; this shouldn't be possible?")
            quit()
    
    # Write the script file
    script = f'''#!/bin/bash -l
#PBS -N prep_reads
#PBS -l walltime=120:00:00
#PBS -l mem=25G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

TRIMMEDDIR={os.path.abspath(args.readsDir)}
'''
    script += "\n".join(catCmds)

    with open(args.outputFileName, "w") as fileOut:
        fileOut.write(script)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
