#! python3
# md5summer.py
# This script automates the generation of a qsub script that will run md5sum on all files
# with a specified suffix in a given directory.

import os, argparse

def validate_args(args):
    # Validate input location
    if not os.path.isdir(args.filesDirectory):
        raise FileNotFoundError(f"Input directory '{args.filesDirectory}' does not exist or is not a directory.")
    args.filesDirectory = os.path.abspath(args.filesDirectory)

def locate_files(filesDirectory, fileSuffix):
    locatedFiles = []
    for file in os.listdir(filesDirectory):
        if file.endswith(fileSuffix):
            locatedFiles.append(os.path.join(filesDirectory, file))
    if locatedFiles == []:
        raise FileNotFoundError(f"No files with suffix '{fileSuffix}' found in directory '{filesDirectory}'.")
    return locatedFiles

def create_md5sum_cmd_file(fileNames, outputFileName="cmd_md5sums.txt"):    
    with open(outputFileName, "w") as fileOut:
        for fileName in fileNames:
            fileOut.write(f"md5sum {fileName} > {fileName}.md5\n")

def create_shell_script(cmdFile, numJobs, outputFileName="run_md5sums.sh", walltime="03:00:00"):
    # Specify hard-coded script features
    jobname = "md5sum"
    mem = "10G"
    
    # Format the PBS -J line
    if numJobs > 1:
        jobLine = f"#PBS -J 1-{numJobs}"
    else:
        jobLine = "PBS_ARRAY_INDEX=1"
    
    # Setup the script's contents
    formatStr = """#!/bin/bash -l
#PBS -N {jobname}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus=1
{jobLine}

cd $PBS_O_WORKDIR

eval $(cat {cmdFile} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1)
"""
    # Write to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(formatStr.format(
            jobname=jobname,
            walltime=walltime,
            mem=mem,
            jobLine=jobLine,
            cmdFile=cmdFile
        ))

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s receives a directory containing files with a specified suffix
    and automates the generation of a qsub script that will run md5sum on all files
    with that suffix.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="filesDirectory",
                   required=True,
                   help="Input directory containing files")
    p.add_argument("-s", dest="filesSuffix",
                   required=True,
                   help="Suffix which uniquely identifies all relevant files")
    # Opts
    p.add_argument("--cmd", dest="cmdFileName",
                   required=False,
                   help="Optionally specify list file name (default==cmd_md5sums.txt)",
                   default="cmd_md5sums.txt")
    p.add_argument("--sh", dest="shellScriptName",
                   required=False,
                   help="Optionally specify shell script name (default==run_md5sums.sh)",
                   default="run_md5sums.sh")
    p.add_argument("--walltime", dest="walltime",
                   required=False,
                   help="Optionally specify the walltime for the job (default==03:00:00)",
                   default="03:00:00")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate files from the directory
    filesForMd5 = locate_files(args.filesDirectory, args.filesSuffix)
    
    # Write the cmd file
    create_md5sum_cmd_file(filesForMd5, outputFileName=args.cmdFileName)
    
    # Write the shell script to run md5sum on the located files
    create_shell_script(args.cmdFileName, len(filesForMd5),
                        outputFileName=args.shellScriptName,
                        walltime=args.walltime)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
