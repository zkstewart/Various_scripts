#! python3
# cnvnator_copynum_pipeline.py
# Script to automatically run CNVnator analysis for the
# prediction of gene copy number variation in a set of BAM files
# using a GFF3 file and a reference genome FASTA file.

import os, argparse, sys, platform, subprocess, re
from Bio import SeqIO

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
    if args.cnvnator is None:
        args.cnvnator = ZS_Utility.wsl_which("cnvnator")
        if args.cnvnator is None:
            _not_specified_error("cnvnator")
    else:
        if not os.path.isfile(args.cnvnator):
            _not_found_error("cnvnator", args.cnvnator)
    
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def cnvnator_tree(rootFileName, bamFile, cnvnatorFile):
    '''
    Parameters:
        rootFileName -- a string indicating the name of the .root file to create
        bamFile -- a string indicating the location of a BAM file to setup a CNVnator .root file for
        cnvnatorFile -- a string indicating the location of the CNVnator executable file
    '''
    readsPlacedRegex = re.compile(r"Total of (\d+) reads were placed")
    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(cnvnatorFile)
    cmd += [
        "-root", ZS_Utility.convert_to_wsl_if_not_unix(rootFileName),
        "-tree", ZS_Utility.convert_to_wsl_if_not_unix(bamFile)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_cnvnator = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    natorout, natorerr = run_cnvnator.communicate()
    
    # Check to see if there was an error
    numReadsPlaced = readsPlacedRegex.findall(natorout.decode("utf-8"))
    if (natorerr.decode("utf-8") != "") or (numReadsPlaced == []):
        raise Exception(("ERROR: cnvnator_tree encountered an error; have a look " +
                        f'at the stdout ({natorout.decode("utf-8")}) and stderr ' + 
                        f'({natorerr.decode("utf-8")}) to make sense of this.'))
    elif int(numReadsPlaced[0]) == 0:
        raise Exception(f"ERROR: cnvnator_tree placed 0 reads; this is likely an error.")

def cnvnator_his(rootFileName, explosionDir, windowSize, cnvnatorFile):
    '''
    Parameters:
        rootFileName -- a string indicating the name of the .root file to get data histograms for
        explosionDir -- a string indicating the location of the exploded genome FASTA
        windowSize -- an integer indicating the size of the window to use for CNVnator
        cnvnatorFile -- a string indicating the location of the CNVnator executable file
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(cnvnatorFile)
    cmd += [
        "-root", ZS_Utility.convert_to_wsl_if_not_unix(rootFileName),
        "-his", str(windowSize),
        "-d", ZS_Utility.convert_to_wsl_if_not_unix(explosionDir)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_cnvnator = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    natorout, natorerr = run_cnvnator.communicate()
    
    # Check to see if there was an error
    if (natorout.decode("utf-8") == "") or (natorerr.decode("utf-8") != ""):
        raise Exception(("ERROR: cnvnator_his encountered an error; have a look " +
                        f'at the stdout ({natorout.decode("utf-8")}) and stderr ' + 
                        f'({natorerr.decode("utf-8")}) to make sense of this.'))

def cnvnator_stat(rootFileName, windowSize, cnvnatorFile):
    '''
    Parameters:
        rootFileName -- a string indicating the name of the .root file to get data histograms for
        windowSize -- an integer indicating the size of the window to use for CNVnator
        cnvnatorFile -- a string indicating the location of the CNVnator executable file
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(cnvnatorFile)
    cmd += [
        "-root", ZS_Utility.convert_to_wsl_if_not_unix(rootFileName),
        "-stat", str(windowSize)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_cnvnator = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    natorout, natorerr = run_cnvnator.communicate()
    
    # Check to see if there was an error
    if not "Average RD per bin" in natorout.decode("utf-8"):
        raise Exception(("ERROR: cnvnator_stat encountered an error; have a look " +
                        f'at the stdout ({natorout.decode("utf-8")}) and stderr ' + 
                        f'({natorerr.decode("utf-8")}) to make sense of this.'))

def cnvnator_partition(rootFileName, windowSize, cnvnatorFile):
    '''
    Parameters:
        rootFileName -- a string indicating the name of the .root file to get data histograms for
        windowSize -- an integer indicating the size of the window to use for CNVnator
        cnvnatorFile -- a string indicating the location of the CNVnator executable file
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(cnvnatorFile)
    cmd += [
        "-root", ZS_Utility.convert_to_wsl_if_not_unix(rootFileName),
        "-partition", str(windowSize)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_cnvnator = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    natorout, natorerr = run_cnvnator.communicate()
    
    # Check to see if there was an error
    if (not "Partitioning RD signal for" in natorout.decode("utf-8")) or (natorerr.decode("utf-8") != ""):
        raise Exception(("ERROR: cnvnator_partition encountered an error; have a look " +
                        f'at the stdout ({natorout.decode("utf-8")}) and stderr ' + 
                        f'({natorerr.decode("utf-8")}) to make sense of this.'))

def cnvnator_call(rootFileName, windowSize, callFileName, cnvnatorFile):
    '''
    Parameters:
        rootFileName -- a string indicating the name of the .root file to get data histograms for
        windowSize -- an integer indicating the size of the window to use for CNVnator
        callFileName -- a string indicating the output file name for CNV calls to be written to
        cnvnatorFile -- a string indicating the location of the CNVnator executable file
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(cnvnatorFile)
    cmd += [
        "-root", ZS_Utility.convert_to_wsl_if_not_unix(rootFileName),
        "-call", str(windowSize)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_cnvnator = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    natorout, natorerr = run_cnvnator.communicate()
    
    # Check to see if there was an error
    if natorerr.decode("utf-8") != "":
        raise Exception(("ERROR: cnvnator_call encountered an error; have a look " +
                        f'at the stdout ({natorout.decode("utf-8")}) and stderr ' + 
                        f'({natorerr.decode("utf-8")}) to make sense of this.'))
    
    # Write call stdout to file
    with open(callFileName, "w") as fileOut:
        fileOut.write(natorout.decode("utf-8"))

def parse_cnvnator_calls(callFile, evalueCutoff=1e-5):
    '''
    Parameters:
        callFile -- a string indicating the location of a CNVnator call file
        evalueCutoff -- a float indicating the E-value significance of the T-test
                        for a CNV to be called; default == 1e-5
    Returns:
        callsList -- a list of lists with format like:
                     [
                         [chromosome, start, end, normalisedRD],
                         ...,
                     ]
    '''
    MAGIC_NUM = 0.1
    callsList = []
    
    with open(callFile, "r") as fileIn:
        for line in fileIn:
            cnvType, coords, cnvSize, normalisedRD, eval1, \
                eval2, eval3, eval4, q0 = line.rstrip("\r\n ").split("\t")
            
            # Skip if E-value 1 (T-test statistic) does not meet significance
            if float(eval1) > evalueCutoff:
                continue
            
            # Extract start and end from coords
            chromosome, coords = coords.split(":")
            start, end = coords.split("-")
            
            # Skip if RD doesn't diverge enough from 1
            normalisedRD = float(normalisedRD)
            if ((round((normalisedRD + MAGIC_NUM)*2) / 2) == 1) or ((round((normalisedRD - MAGIC_NUM)*2) / 2) == 1):
                continue
            
            # Store values
            callsList.append([chromosome, int(start), int(end), normalisedRD])
    
    return callsList

def find_genes_within_cnvs(callsList, gff3Obj, overlapProportion=0.5):
    '''
    Parameters:
        callsList -- a list of lists with format like:
                     [
                         [chromosome, start, end, normalisedRD],
                         ...,
                     ]
        gff3Obj -- a ZS_GFF3IO.GFF3 object with NCLS indexing applied
        overlapProportion -- a float indicating the proportion of a gene that must
                             be within a CNV to be considered; default == 0.5
    Returns:
        callsList -- a list of lists with format like:
                     [
                         [chromosome, start, end, normalisedRD],
                         ...,
                     ]
    '''
    cnvGenes = []
    for chromosome, start, end, normalisedRD in callsList:
        # Find genes overlapping this region via NCLS
        overlappingGenes = gff3Obj.ncls_finder(start, end, "contig", chromosome)
        if overlappingGenes == []:
            continue
        
        # Calculate the proportion of the gene(s) that are within the CNV
        foundGenes = []
        for geneFeature in overlappingGenes:
            geneStart, geneEnd = geneFeature.coords
            if geneStart <= end and geneEnd >= start:
                thisOverlap = min(end, geneEnd) - max(start, geneStart) + 1
                thisProportion = thisOverlap / (end - start + 1)
                if thisProportion >= overlapProportion:
                    foundGenes.append(geneFeature.ID)
        
        # Store any identified genes with their CNV information
        for foundGeneID in foundGenes:
            cnvGenes.append([foundGeneID, normalisedRD])
    
    return cnvGenes

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3, a directory containing an exploded FASTA file, and a directory
    containing one or more BAM files. It will run CNVnator for each sample/BAM file to prevent CNV
    genomic segments. Using the GFF3 it will attempt to identify genes within the CNV segments and
    estimate how many gene copies there are.
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
    # Opts (programs)
    p.add_argument("--cnvnator", dest="cnvnator",
                   required=False,
                   help="""Optionally, specify the cnvnator executable file
                   if it is not discoverable in the path""",
                   default=None)
    # Opts (behavioural)
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
        print(f"bam symlinking has already been performed; skipping.")
    
    # Explode FASTA file into working directory
    explosionDir = os.path.join(args.outputDirectory, "exploded_fasta")
    os.makedirs(explosionDir, exist_ok=True)
    
    if not os.path.exists(os.path.join(args.outputDirectory, "explosion_was_successful.flag")):
        # Load FASTA file in as SeqIO records
        with open(args.fastaFile, "r") as fileIn:
            records = SeqIO.parse(fileIn, "fasta")
            
            # Iterate over and create file for each contig/chromosome/'record'
            for record in records:
                seqFile = os.path.join(explosionDir, record.id + ".fa")
                if not os.path.exists(seqFile):
                    with open(seqFile, "w") as fileOut:
                        fileOut.write(record.format("fasta"))
        
        open(os.path.join(args.outputDirectory, "explosion_was_successful.flag"), "w").close()
    else:
        print(f"FASTA explosion has already been performed; skipping.")
    
    # Set up the working directory structure
    rootDir = os.path.abspath(os.path.join(args.outputDirectory, "root_files"))
    os.makedirs(rootDir, exist_ok=True)
    
    callsDir = os.path.abspath(os.path.join(args.outputDirectory, "calls"))
    os.makedirs(callsDir, exist_ok=True)
    
    # Run CNVnator for each BAM file
    WINDOW_SIZES = [10000, 100000]
    if not os.path.exists(os.path.join(args.outputDirectory, "cnvnator_was_successful.flag")):
        # Locate symlinked BAM files
        bamFiles = [ os.path.join(bamsDir, file) for file in os.listdir(bamsDir) ]
        
        # Iterate through each BAM file
        for bamFile in bamFiles:
            bamBase = os.path.basename(bamFile).replace(args.bamSuffix, "")
            rootFile = os.path.join(rootDir, bamBase + ".root")
            
            # Run tree
            treeFlag = os.path.join(rootDir, bamBase + ".tree_was_successful.flag")
            if not os.path.exists(treeFlag):
                cnvnator_tree(rootFile, bamFile, args.cnvnator)
                open(treeFlag, "w").close()
            else:
                print(f"cnvnator 'tree' already been run for '{bamBase}'; skipping.")
            
            # Run his
            hisFlag = os.path.join(rootDir, bamBase + ".his_was_successful.flag")
            if not os.path.exists(hisFlag):
                for windowSize in WINDOW_SIZES:
                    cnvnator_his(rootFile, explosionDir, windowSize, args.cnvnator)
                open(hisFlag, "w").close()
            else:
                print(f"cnvnator 'his' already been run for '{bamBase}'; skipping.")
            
            # Run stat
            statFlag = os.path.join(rootDir, bamBase + ".stat_was_successful.flag")
            if not os.path.exists(statFlag):
                for windowSize in WINDOW_SIZES:
                    cnvnator_stat(rootFile, windowSize, args.cnvnator)
                open(statFlag, "w").close()
            else:
                print(f"cnvnator 'stat' already been run for '{bamBase}'; skipping.")
            
            # Run partition
            partitionFlag = os.path.join(rootDir, bamBase + ".partition_was_successful.flag")
            if not os.path.exists(partitionFlag):
                for windowSize in WINDOW_SIZES:
                    cnvnator_partition(rootFile, windowSize, args.cnvnator)
                open(partitionFlag, "w").close()
            else:
                print(f"cnvnator 'partition' already been run for '{bamBase}'; skipping.")
            
            # Run call
            callFlag = os.path.join(rootDir, bamBase + ".call_was_successful.flag")
            if not os.path.exists(callFlag):
                for windowSize in WINDOW_SIZES:
                    callFile = os.path.join(callsDir, f"{bamBase}.calls.{windowSize}.tsv")
                    if not os.path.exists(callFile):
                        cnvnator_call(rootFile, windowSize, callFile, args.cnvnator)
                open(callFlag, "w").close()
            else:
                print(f"cnvnator 'call' already been run for '{bamBase}'; skipping.")
        
        open(os.path.join(args.outputDirectory, "cnvnator_was_successful.flag"), "w").close()
    else:
        print(f"cnvnator has already been performed; skipping.")
    
    # Parse the GFF3 file with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse = not args.relaxedParsing)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Tabulate gene copy number results
    CHOSEN_WINDOW_SIZE = 10000
    
    copynumTableFile = os.path.join(args.outputDirectory, "gene_copy_numbers.tsv")
    if not os.path.exists(os.path.join(args.outputDirectory, "tabulation_was_successful.flag")):
        # Locate symlinked BAM files
        bamFiles = [ os.path.join(bamsDir, file) for file in os.listdir(bamsDir) ]
        
        # Begin generating result tabulation
        with open(copynumTableFile, "w") as fileOut:
            # Write header line
            fileOut.write("sampleid\tgeneid\tmultiplication\n")
            
            # Iterate through each sample's call file
            for bamFile in bamFiles:
                samplePrefix = os.path.basename(bamFile).replace(args.bamSuffix, "")
                
                # Parse the call file for CNV regions
                callFile = os.path.join(callsDir, f"{samplePrefix}.calls.{CHOSEN_WINDOW_SIZE}.tsv")
                callsList = parse_cnvnator_calls(callFile)
                
                # Find genes contained within a call region
                cnvGenes = find_genes_within_cnvs(callsList, gff3Obj)
                
                # Write any results to file
                for geneID, normalisedRD in cnvGenes:
                    multiplicationFactor = round(normalisedRD*2) / 2
                    fileOut.write(f"{samplePrefix}\t{geneID}\t{multiplicationFactor}\n")
        
        open(os.path.join(args.outputDirectory, "tabulation_was_successful.flag"), "w").close()
    else:
        print(f"tabulation has already been performed; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
