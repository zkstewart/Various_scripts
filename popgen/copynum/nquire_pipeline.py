#! python3
# nquire_pipeline.py
# Script to automatically run nQuire analysis for the
# prediction of genome ploidy values in a set of BAM files
# using a reference genome FASTA.

import os, argparse, sys, platform, subprocess

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find Function_packages
from Function_packages import ZS_Utility

# Define functions
def validate_args(args):
    def _not_specified_error(program):
        print(f"ERROR: {program} not discoverable in your system PATH and was not specified as an argument.")
        quit()
    def _not_found_error(program, path):
        print(f"ERROR: {program} was not found at the location indicated ('{path}')")
        quit()
    
    # Validate input file locations
    if not os.path.isdir(args.bamDirectory):
        print('I am unable to locate the BAM directory (' + args.bamDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate BAM suffix
    foundABAM = False
    for file in os.listdir(args.bamDirectory):
        if file.endswith(args.bamSuffix):
            foundABAM = True
            break
    if not foundABAM:
        print(f'No BAM files with suffix "{args.bamSuffix}" found in directory "{args.bamDirectory}"')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate program discoverability
    if args.nQuire is None:
        args.nQuire = ZS_Utility.wsl_which("nQuire")
        if args.nQuire is None:
            _not_specified_error("nQuire")
    else:
        if not os.path.isfile(args.nQuire):
            _not_found_error("nQuire", args.nQuire)
    
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def run_nQuire_create(bamFile, outputBaseName, nQuirePath):
    '''
    Parameters:
        bamFile -- a string indicating the location of the BAM file to process
        outputBaseName -- a string indicating the basename for nQuire to write
                          the .bin file to (note: nQuire will add the .bin extension)
        nQuirePath -- a string indicating the location of the nQuire executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(nQuirePath)
    cmd += [
        "create", "-b", ZS_Utility.convert_to_wsl_if_not_unix(bamFile),
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputBaseName)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_create = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    createout, createerr = run_create.communicate()
    
    # Check file outputs to see if there was an error
    if (not createout.decode("utf-8") == "") or (not createerr.decode("utf-8") == ""):
        raise Exception(("ERROR: run_nQuire_create encountered an error; have a look " +
                        f'at the stdout ({createout.decode("utf-8")}) and stderr ' + 
                        f'({createerr.decode("utf-8")}) to make sense of this.'))

def run_nQuire_lrdmodel(binFile, outputFileName, nQuirePath):
    '''
    Parameters:
        binFile -- a string indicating the location of the nQuire create .bin file
        outputBaseName -- a string indicating the location to write lrdmodel results to
        nQuirePath -- a string indicating the location of the nQuire executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(nQuirePath)
    cmd += [
        "lrdmodel", ZS_Utility.convert_to_wsl_if_not_unix(binFile)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_lrdmodel = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    lrdmodelout, lrdmodelerr = run_lrdmodel.communicate()
    
    # Check stderr for errors
    if lrdmodelerr.decode("utf-8") != "":
        raise Exception(("ERROR: run_nQuire_lrdmodel encountered an error; have a look " +
                        f'at the stdout ({lrdmodelout.decode("utf-8")}) and stderr ' + 
                        f'({lrdmodelerr.decode("utf-8")}) to make sense of this.'))
    
    # Check that stdout contains the bin file name
    elif not binFile in lrdmodelout.decode("utf-8"):
        raise Exception(("ERROR: run_nQuire_lrdmodel encountered an error; have a look " +
                        f'at the stdout ({lrdmodelout.decode("utf-8")}) and stderr ' + 
                        f'({lrdmodelerr.decode("utf-8")}) to make sense of this.'))
    
    # Write stdout to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(lrdmodelout.decode("utf-8"))

def parse_lrdmodel_file(lrdModelFile):
    '''
    Parameters:
        lrdModelFile -- a string indicating the location of the lrdmodel results file
    Returns:
        ploidyLikelihoodOrder -- a list of strings indicating the ploidy likelihood ordered
                                 from most likely to least likely; strings are "diploid",
                                 "triploid", and "tetraploid"
    '''
    with open(lrdModelFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle first line
            if firstLine:
                assert sl == ["file", "free", "dip", "tri", "tet", "d_dip", "d_tri", "d_tet"], \
                    f"ERROR: lrdmodel file '{lrdModelFile}' is not formatted correctly; header isn't as expected"
                firstLine = False
            # Handle subsequent line
            else:
                f, free, dip, tri, tet, d_dip, d_tri, d_tet = sl
                
                # Get the ordered likelihoods
                ploidyLikelihoodOrder = [[float(dip), "diploid"], [float(tri), "triploid"], [float(tet), "tetraploid"]]
                ploidyLikelihoodOrder.sort(key = lambda x: x[0])
                
                return [x[1] for x in ploidyLikelihoodOrder]

def tabulate_ploidy_estimates(lrdmodelDir, outputFileName):
    '''
    Parameters:
        lrdmodelDir -- a string indicating the location where nQuire lrdmodel stdout
                       results are located
        outputFileName -- a string indicating the location to write the ploidy
                          estimates to
    '''
    with open(outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("sample_id\tlikely_ploidy\tnext_likeliest\tleast_likely\n")
        
        # Iterate through each sample folder
        for lrdResultFile in os.listdir(lrdmodelDir):
            if lrdResultFile.endswith(".lrdmodel"):
                # Parse the lrdmodel results file
                first, second, third = parse_lrdmodel_file(os.path.join(lrdmodelDir, lrdResultFile))
                
                # Write to file
                lrdPrefix = lrdResultFile.rsplit(".", maxsplit=1)[0]
                fileOut.write(f"{lrdPrefix}\t{first}\t{second}\t{third}\n")

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing one or more BAM files. It will run nQuire
    for each BAM file and tabulate the most likely ploidy for each sample.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-b", dest="bamDirectory",
                   required=True,
                   help="Input directory containing BAM file(s)")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (programs)
    p.add_argument("--nQuire", dest="nQuire",
                   required=False,
                   help="""Optionally, specify the nQuire executable file
                   if it is not discoverable in the path""",
                   default=None)
    # Opts (behavioural)
    p.add_argument("--bamSuffix", dest="bamSuffix",
                   required=False,
                   help="""Indicate the suffix of the BAM files to look for;
                   default == '.sorted.bam'""",
                   default=".sorted.bam")
    
    args = p.parse_args()
    validate_args(args)
    
    # Symlink BAM files into working directory
    bamsDir = os.path.join(args.outputDirectory, "bams")
    os.makedirs(bamsDir, exist_ok=True)
    
    if not os.path.exists(os.path.join(args.outputDirectory, "bam_symlink_was_successful.flag")):
        # Locate BAM files
        bamFiles = [
            os.path.join(args.bamDirectory, file)
            for file in os.listdir(args.bamDirectory)
            if file.endswith(args.bamSuffix)
        ]
        
        # Symlink each BAM file
        for bamFile in bamFiles:
            outputBamFile = os.path.join(bamsDir, os.path.basename(bamFile))
            if not os.path.exists(outputBamFile):
                os.symlink(bamFile, outputBamFile)
        open(os.path.join(args.outputDirectory, "bam_symlink_was_successful.flag"), "w").close()
    else:
        print(f"bam symlinking has already been performing; skipping.")
    
    # Run nQuire for each BAM file
    nquireBinDir = os.path.join(args.outputDirectory, "nquire_bin_files")
    os.makedirs(nquireBinDir, exist_ok=True)
    
    if not os.path.exists(os.path.join(args.outputDirectory, "nquire_was_successful.flag")):
        # Locate symlinked BAM files
        bamFiles = [ os.path.join(bamsDir, file) for file in os.listdir(bamsDir) ]
        
        # Iterate through each BAM file
        for bamFile in bamFiles:
            bamBase = os.path.basename(bamFile).replace(args.bamSuffix, "")
            
            # Run nQuire create
            binFileName = os.path.join(nquireBinDir, bamBase + ".bin")
            if not os.path.exists(binFileName):
                run_nQuire_create(bamFile, binFileName.rsplit(".", maxsplit=1)[0], args.nQuire) # remove .bin since nQuire adds it automatically
            else:
                print(f"nQuire create has already been run for '{bamBase}'; skipping.")
            
            # Run nQuire lrdmodel
            lrdmodelFileName = os.path.join(args.outputDirectory, bamBase + ".lrdmodel")
            if not os.path.exists(lrdmodelFileName):
                run_nQuire_lrdmodel(binFileName, lrdmodelFileName, args.nQuire)
            else:
                print(f"nQuire lrdmodel has already been run for '{bamBase}'; skipping.")
        
        open(os.path.join(args.outputDirectory, "nquire_was_successful.flag"), "w").close()
    else:
        print(f"nQuire has already been performed; skipping.")
    
    # Tabulate sample ploidy results
    ploidyTableFile = os.path.join(args.outputDirectory, "ploidy_estimates.tsv")
    if not os.path.exists(ploidyTableFile) or not \
        os.path.exists(os.path.join(args.outputDirectory, "ploidy_tabulation_was_successful.flag")):            
            tabulate_ploidy_estimates(args.outputDirectory, ploidyTableFile)
            open(os.path.join(args.outputDirectory, "ploidy_tabulation_was_successful.flag"), "w").close()
    else:
        print(f"ploidy tabulation file has already been generated; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
