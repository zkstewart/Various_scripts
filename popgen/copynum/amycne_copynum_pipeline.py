#! python3
# amycne_copynum_pipeline.py
# Script to automatically run AMYCNE analysis for the
# prediction of gene copy number variation in a set of BAM files
# using a GFF3 file and a reference genome FASTA file.

import os, argparse, sys, platform, subprocess

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find Function_packages
from Function_packages import ZS_Utility, ZS_GFF3IO

# Define functions
def validate_args(args):
    def _not_specified_error(program):
        print(f"ERROR: {program} not discoverable in your system PATH and was not specified as an argument.")
        quit()
    def _not_found_error(program, path):
        print(f"ERROR: {program} was not found at the location indicated ('{path}')")
        quit()
    
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the genome FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.bamDirectory):
        print('I am unable to locate the BAM directory (' + args.bamDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.amycneDirectory):
        print('I am unable to locate the AMYCNE directory (' + args.amycneDirectory + ')')
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
    for pyScript in ["AMYCNE.py", "Generate_GC_tab.py"]:
        if not os.path.isfile(os.path.join(args.amycneDirectory, pyScript)):
            _not_found_error(pyScript, args.amycneDirectory)
    
    if not os.path.isfile(args.python2):
        _not_found_error("python", args.python2)
    
    if args.tiddit is None:
        args.tiddit = ZS_Utility.wsl_which("tiddit")
        if args.tiddit is None:
            _not_specified_error("tiddit")
    else:
        if not os.path.isfile(args.tiddit):
            _not_found_error("tiddit", args.tiddit)
    
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def generate_gc_file(fastaFile, windowSize, outputFile, python2Exe, amycneDirectory):
    '''
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file to generate the
                     GC content file from
        windowSize -- an integer indicating the window size to calculate GC content over
        outputFile -- a string indicating the location to write the output file to
        python2Exe -- a string indicating the location of the python2 executable
        amycneDirectory -- a string indicating the location of the AMYCNE directory
                           containing Generate_GC_tab.py
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(python2Exe)
    cmd += [
        ZS_Utility.convert_to_wsl_if_not_unix(os.path.join(amycneDirectory, "Generate_GC_tab.py")),
        "--fa", ZS_Utility.convert_to_wsl_if_not_unix(fastaFile),
        "--size", str(windowSize)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_gc_tab = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    gcout, gcerr = run_gc_tab.communicate()
    
    # Check to see if there was an error
    if (gcout.decode("utf-8") == "") or (gcerr.decode("utf-8") != ""):
        raise Exception(("ERROR: generate_gc_file encountered an error; have a look " +
                        f'at the stdout ({gcout.decode("utf-8")}) and stderr ' + 
                        f'({gcerr.decode("utf-8")}) to make sense of this.'))
    
    # Write the stdout to file
    with open(outputFile, "w") as gcFile:
        gcFile.write(gcout.decode("utf-8"))

def generate_tiddit_file(fastaFile, bamFile, windowSize, outputPrefix, tidditPath):
    '''
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file to generate the
                     GC content file from
        
        windowSize -- an integer indicating the window size to calculate GC content over
        outputPrefix -- a string indicating the prefix for output file generation by tiddit
        tidditPath -- a string indicating the location of the tiddit executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(tidditPath)
    cmd += [
        "--cov", "--bam", ZS_Utility.convert_to_wsl_if_not_unix(bamFile),
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputPrefix),
        "-z", str(windowSize), "--ref", ZS_Utility.convert_to_wsl_if_not_unix(fastaFile)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_tiddit = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    tidditout, tidditerr = run_tiddit.communicate()
    
    # Check to see if there was an error
    if (tidditout.decode("utf-8") != "") or (tidditerr.decode("utf-8") != ""):
        raise Exception(("ERROR: generate_tiddit_file encountered an error; have a look " +
                        f'at the stdout ({tidditout.decode("utf-8")}) and stderr ' + 
                        f'({tidditerr.decode("utf-8")}) to make sense of this.'))

def run_amycne_per_gene(gcFile, tidditFile, gff3Obj, outputFileName, python2Exe, amycneDirectory):
    '''
    Parameters:
        gcFile -- a string indicating the location of the gc_content.tab file generated by
                  AMYCNE's Generate_GC_tab.py
        tidditFile -- a string indicating the location of the tiddit file
        gff3Obj -- a ZS_GFF3IO.GFF3 object containing the parsed GFF3 file
        outputFileName -- a string indicating the file name to write AMYCNE copy number
                          predictions to
        python2Exe -- a string indicating the location of the python2 executable
        amycneDirectory -- a string indicating the location of the AMYCNE directory
                           containing AMYCNE.py
    '''
    assert "gene" in gff3Obj.featureTypes, "ERROR: GFF3 file does not contain 'gene' features!?!"
    
    # Iterate over GFF3 gene features
    resultLines = ["gene\tbins\tused_bin_ratio\tref_coverage\tCN_raw\t95%CI\tCN_rounded\tregion_command"]
    for geneFeature in gff3Obj.types["gene"]:
        # Format the gene feature for AMYCNE
        geneRegion = f"{geneFeature.contig}:{geneFeature.start}-{geneFeature.end}"
        
        # Construct the cmd for subprocess
        cmd = ZS_Utility.base_subprocess_cmd(python2Exe)
        cmd += [
            ZS_Utility.convert_to_wsl_if_not_unix(os.path.join(amycneDirectory, "AMYCNE.py")),
            "--genotype", "--gc", ZS_Utility.convert_to_wsl_if_not_unix(gcFile),
            "--coverage", ZS_Utility.convert_to_wsl_if_not_unix(tidditFile),
            "--R", geneRegion
        ]
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
        # Run the command
        run_amycne = subprocess.Popen(cmd, shell = True,
                                    stdout = subprocess.PIPE,
                                    stderr = subprocess.PIPE)
        amycneout, amycneerr = run_amycne.communicate()
        
        # Check to see if there was an error
        if (amycneout.decode("utf-8") == "") or (amycneerr.decode("utf-8") != ""):
            raise Exception(("ERROR: run_amycne_per_gene encountered an error; have a look " +
                            f'at the stdout ({amycneout.decode("utf-8")}) and stderr ' + 
                            f'({amycneerr.decode("utf-8")}) to make sense of this.'))
        
        # Extract information from the AMYCNE output
        amycneResult = amycneout.decode("utf-8").strip("\r\n ").split("\n")
        assert len(amycneResult) == 2, ("ERROR: AMYCNE output is not as expected, " +
                                        f"should have two lines but I see '{amycneResult}'")
        amycneResult = amycneResult[-1]
        
        # Reformat the result and store in the results list
        amycneResult = [geneFeature.ID] + [ value for index, value in enumerate(amycneResult.split("\t")) if index != 0 ]
        resultLines.append("\t".join(amycneResult))
    
    # Write the results to file
    with open(outputFileName, "w") as amycneFile:
        amycneFile.write("\n".join(resultLines))

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3, a genome FASTA, and a directory containing one or more
    BAM files. It will run AMYCNE for each sample/BAM file on each gene in the GFF3 file to
    predict gene copy number variation. The output will be a table of gene copy number estimates
    for each sample. Note: This script is to be run on python3 but AMYCNE is a python2 script,
    hence you need to indicate the location of the python2 executable file.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome FASTA file")
    p.add_argument("-b", dest="bamDirectory",
                   required=True,
                   help="Input directory containing BAM file(s)")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    p.add_argument("-p", dest="python2",
                   required=False,
                   help="""Specify the python2 executable file""",
                   default=None)
    p.add_argument("-a", dest="amycneDirectory",
                   required=True,
                   help="""Specify the location of the AMYCNE directory containing the
                   AMYCNE.py and Generate_GC_tab.py files""")
    # Opts (programs)
    p.add_argument("--tiddit", dest="tiddit",
                   required=False,
                   help="""Optionally, specify the tiddit executable file
                   if it is not discoverable in the path""",
                   default=None)
    # Opts (behavioural)
    p.add_argument("--window", dest="windowSize",
                   required=False,
                   type=int,
                   help="""Optionally specify the window size for GC and tiddit calculations;
                   default == 100000""",
                   default=100000)
    p.add_argument("--bamSuffix", dest="bamSuffix",
                   required=False,
                   help="""Indicate the suffix of the BAM files to look for;
                   default == '.sorted.bam'""",
                   default=".sorted.bam")
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    
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
    
    # Set up the working directory structure
    coverageDir = os.path.abspath(os.path.join(args.outputDirectory, "coverage_files"))
    os.makedirs(coverageDir, exist_ok=True)
    
    amycneDir = os.path.abspath(os.path.join(args.outputDirectory, "amycne_outputs"))
    os.makedirs(amycneDir, exist_ok=True)
    
    # Generate gc_content.tab file
    gcFile = os.path.join(coverageDir, "gc_content.tab")
    if not os.path.exists(gcFile):
        generate_gc_file(args.fastaFile, args.windowSize, gcFile,
                         args.python2, args.amycneDirectory)
    else:
        print(f"gc_content.tab file has already been generated; skipping.")
    
    # Parse the GFF3 file
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, relaxed=args.relaxedParsing)
    
    # Run AMYCNE for each BAM file
    if not os.path.exists(os.path.join(args.outputDirectory, "amycne_was_successful.flag")):
        # Locate symlinked BAM files
        bamFiles = [ os.path.join(bamsDir, file) for file in os.listdir(bamsDir) ]
        
        # Iterate through each BAM file
        for bamFile in bamFiles:
            bamBase = os.path.basename(bamFile).replace(args.bamSuffix, "")
            
            # Run tiddit
            tidditPrefix = os.path.join(coverageDir, bamBase + ".tiddit")
            tidditFileName = tidditPrefix + ".tab"
            if not os.path.exists(tidditFileName):
                generate_tiddit_file(args.fastaFile, bamFile, args.windowSize,
                                     tidditPrefix, args.tiddit)
            else:
                print(f"tiddit has already been run for '{bamBase}'; skipping.")
            
            # Run AMYCNE for each gene
            amycneFileName = os.path.join(amycneDir, bamBase + ".tsv")
            if not os.path.exists(amycneFileName):
                run_amycne_per_gene(gcFile, tidditFileName, gff3Obj, amycneFileName,
                                    args.python2, args.amycneDirectory)
            else:
                print(f"AMYCNE has already been run for '{bamBase}'; skipping.")
        
        open(os.path.join(args.outputDirectory, "amycne_was_successful.flag"), "w").close()
    else:
        print(f"AMYCNE has already been performed; skipping.")
    
    # # Tabulate gene copy number results
    # copynumTableFile = os.path.join(args.outputDirectory, "gene_copy_numbers.tsv")
    # if not os.path.exists(copynumTableFile) or not \
    #     os.path.exists(os.path.join(args.outputDirectory, "copynum_tabulation_was_successful.flag")):            
    #         tabulate_copynum_estimates(freecBaseDir, copynumTableFile) # not implemented yet
    #         open(os.path.join(args.outputDirectory, "copynum_tabulation_was_successful.flag"), "w").close()
    # else:
    #     print(f"copy number table file has already been generated; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
