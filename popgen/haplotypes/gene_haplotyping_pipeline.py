#! python3
# gene_haplotyping_pipeline.py
# Script to take in a GFF3 file and a directory containing
# sorted, indexed BAM files to 

import os, argparse, sys, platform, subprocess, shutil
from intervaltree import IntervalTree
from Bio import SeqIO
from pyfaidx import Fasta

sys.path.append(os.path.dirname(os.path.abspath(__file__))) # the current dir is where we find haplotypes
from predict_variant_effects import convert_vcf_snps_to_cds_snps

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_GFF3IO, ZS_Utility, ZS_VCFIO

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
    if args.samtools is None:
        args.samtools = ZS_Utility.wsl_which("samtools")
        if args.samtools is None:
            _not_specified_error("samtools")
    else:
        if not os.path.isfile(args.samtools):
            _not_found_error("samtools", args.samtools)
    
    if args.bcftools is None:
        args.bcftools = ZS_Utility.wsl_which("bcftools")
        if args.bcftools is None:
            _not_specified_error("bcftools")
    else:
        if not os.path.isfile(args.bcftools):
            _not_found_error("bcftools", args.bcftools)
    
    if args.picard is None:
        args.picard = ZS_Utility.wsl_which("picard")
        if args.picard is None:
            _not_specified_error("picard")
    else:
        if not os.path.isfile(args.picard):
            _not_found_error("picard", args.picard)
    
    if args.whatshap is None:
        args.whatshap = ZS_Utility.wsl_which("whatshap")
        if args.whatshap is None:
            _not_specified_error("whatshap")
    else:
        if not os.path.isfile(args.whatshap):
            _not_found_error("whatshap", args.whatshap)
    
    if args.bgzip is None:
        args.bgzip = ZS_Utility.wsl_which("bgzip")
        if args.bgzip is None:
            _not_specified_error("bgzip")
    else:
        if not os.path.isfile(args.bgzip):
            _not_found_error("bgzip", args.bgzip)
    
    if args.tabix is None:
        args.tabix = ZS_Utility.wsl_which("tabix")
        if args.tabix is None:
            _not_specified_error("tabix")
    else:
        if not os.path.isfile(args.tabix):
            _not_found_error("tabix", args.tabix)
    # Validate numeric inputs
    if 0 > args.plusMinus:
        print("--plusMinus should be 0 or greater")
        quit()
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def gff3_to_coordinatesDict(gff3Obj, plusMinus):
    '''
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object
        plusMinus -- an integer indicating how far to the left and right
                     of a gene region to subset BAMs to
    Returns:
        coordinatesDict -- a dictionary with structure like:
                           {
                                "contig1": [
                                    [start1, end1],
                                    [start2, end2],
                                    ...
                                ],
                                "contig2": [
                                    [start1, end1],
                                    [start2, end2],
                                    ...
                                ],
                                ...
                           }
    '''
    coordinatesDict = {}
    for parentType in gff3Obj.parentTypes:
        for parentFeature in gff3Obj.types[parentType]:
            coordinatesDict.setdefault(parentFeature.contig, [])
            
            # Get the coordinates for subsetting
            start = parentFeature.start - plusMinus
            if start < 1:
                start = 1
            end = parentFeature.end + plusMinus
            
            # Store
            coordinatesDict[parentFeature.contig].append([start, end])
    return coordinatesDict

def merge_coordinates_list(coordinatesList):
    '''
    Parameters:
        coordinatesList -- a list of lists, each containing two integers
                           indicating the start and end of a region
    Returns:
        mergedCoordinatesList -- a list of lists, each containing two integers
                                 indicating the start and end of a region, with
                                 overlapping regions merged
    '''
    tree = IntervalTree.from_tuples(coordinatesList)
    mergedCoordinatesList = []
    for interval in tree:
        mergedCoordinatesList.append([interval.begin, interval.end])
    mergedCoordinatesList.sort()
    return mergedCoordinatesList

def merge_overlapping_coordinates(coordinatesDict):
    '''
    Parameters:
        coordinatesDict -- a dictionary with structure like:
                           {
                                "contig1": [
                                    [start1, end1],
                                    [start2, end2],
                                    ...
                                ],
                                "contig2": [
                                    [start1, end1],
                                    [start2, end2],
                                    ...
                                ],
                                ...
                           }
    Returns:
        correctedCoordinatesDict -- a dictionary with the same structure as the
                                    input dictionary, but with overlapping coordinates
                                    merged.
    '''
    correctedCoordinatesDict = {}
    for contig, coordsList in coordinatesDict.items():
        correctedCoordinatesDict[contig] = merge_coordinates_list(coordsList)
    return correctedCoordinatesDict

def format_samtools_regions(coordinatesDict):
    '''
    Parameters:
        coordinatesDict -- a dictionary with structure like:
                           {
                                "contig1": [
                                    [start1, end1],
                                    [start2, end2],
                                    ...
                                ],
                                "contig2": [
                                    [start1, end1],
                                    [start2, end2],
                                    ...
                                ],
                                ...
                           }
    Returns:
        samtoolsRegions -- a list of strings formatted for samtools
                           like:
                           [
                               "contig1:start1-end1",
                               "contig1:start2-end2",
                               ...
                               "contig2:start1-end1",
                               "contig2:start2-end2",
                               ...
                           ]
    '''
    samtoolsRegions = []
    for contig, coordsList in coordinatesDict.items():
        for coords in coordsList:
            samtoolsRegions.append(f"{contig}:{coords[0]}-{coords[1]}")
    return samtoolsRegions

def samtools_subset_bam(bamFile, outputFile, samtoolsRegions, samtoolsPath):
    '''
    Parameters:
        bamFile -- a string indicating the location of the BAM file to subset
        outputFile -- a string indicating the location to write the subset BAM to
        samtoolsRegions -- a list of strings formatted for samtools
                           like:
                           [
                               "contig1:start1-end1",
                               "contig1:start2-end2",
                               ...
                               "contig2:start1-end1",
                               "contig2:start2-end2",
                               ...
                           ]
        samtoolsPath -- a string indicating the location of the samtools executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(samtoolsPath)
    cmd += [
        "view", "-b", "-h",
        ZS_Utility.convert_to_wsl_if_not_unix(bamFile),
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFile),
        *samtoolsRegions
    ]
    
    # Run the command
    run_samtools_subset = subprocess.Popen(cmd, shell = True,
                                           stdout = subprocess.PIPE,
                                           stderr = subprocess.PIPE)
    subsetout, subseterr = run_samtools_subset.communicate()
    if subsetout.decode("utf-8") != "" and subseterr.decode("utf-8") == "":
        print("WARNING: samtools_subset_bam may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({subsetout.decode("utf-8")})')
    elif subseterr.decode("utf-8") != "":
        raise Exception(("ERROR: samtools_subset_bam encountered an error; have a look " +
                         f'at the stdout > {subsetout.decode("utf-8")} < and stderr ' + 
                         f'> {subseterr.decode("utf-8")} < to make sense of this.'))

def picard_addreplace_readgroups(inputBamFile, outputBamFile, rgID, rgLB, rgPL, rgPU, rgSM, picardPath):
    '''
    Parameters:
        inputBamFile -- a string indicating the location of the input BAM file
        outputBamFile -- a string indicating the location to write the output BAM to
        rgID -- a string indicating the read group ID
        rgLB -- a string indicating the read group library
        rgPL -- a string indicating the read group platform
        rgPU -- a string indicating the read group platform unit
        rgSM -- a string indicating the read group sample
        picardPath -- a string indicating the location of the picard executable
    '''
    # Construct the cmd for subprocess
    if platform.system() == "Windows":
        cmd = ZS_Utility.base_subprocess_cmd("bash") # necessary for executable-ised picard
        cmd.append(ZS_Utility.convert_to_wsl_if_not_unix(picardPath))
    else:
        cmd = ZS_Utility.base_subprocess_cmd(picardPath)
    
    cmd += [
        "AddOrReplaceReadGroups",
        "I=" + ZS_Utility.convert_to_wsl_if_not_unix(inputBamFile),
        "O=" + ZS_Utility.convert_to_wsl_if_not_unix(outputBamFile),
        "RGID=" + rgID,
        "RGLB=" + rgLB,
        "RGPL=" + rgPL,
        "RGPU=" + rgPU,
        "RGSM=" + rgSM
    ]
    
    # Run the command
    run_picard_addreplace = subprocess.Popen(cmd, shell = True,
                                             stdout = subprocess.PIPE,
                                             stderr = subprocess.PIPE)
    addreplaceout, addreplaceerr = run_picard_addreplace.communicate()
    if addreplaceout.decode("utf-8") != "" and addreplaceerr.decode("utf-8") == "":
        print("WARNING: picard_addreplace_readgroups may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({addreplaceout.decode("utf-8")})')
    elif "done" not in addreplaceerr.decode("utf-8"):
        raise Exception(("ERROR: picard_addreplace_readgroups has probably encountered an error; have a look " +
                         f'at the stdout ({addreplaceout.decode("utf-8")}) and stderr ' + 
                         f'({addreplaceerr.decode("utf-8")}) to make sense of this.'))

def samtools_faidx(fastaFile, samtoolsPath):
    '''
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file to index
        samtoolsPath -- a string indicating the location of the samtools executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(samtoolsPath)
    cmd += [
        "faidx",
        ZS_Utility.convert_to_wsl_if_not_unix(fastaFile)
    ]
    
    # Run the command
    run_samtools_index = subprocess.Popen(cmd, shell = True,
                                          stdout = subprocess.PIPE,
                                          stderr = subprocess.PIPE)
    indexout, indexerr = run_samtools_index.communicate()
    if indexout.decode("utf-8") != "" and indexerr.decode("utf-8") == "":
        print("WARNING: samtools_faidx may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({indexout.decode("utf-8")})')
    elif indexerr.decode("utf-8") != "":
        raise Exception(("ERROR: samtools_faidx encountered an error; have a look " +
                         f'at the stdout > {indexout.decode("utf-8")} < and stderr ' + 
                         f'> {indexerr.decode("utf-8")} < to make sense of this.'))

def samtools_index_bam(bamFile, samtoolsPath):
    '''
    Parameters:
        bamFile -- a string indicating the location of the BAM file to subset
        samtoolsPath -- a string indicating the location of the samtools executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(samtoolsPath)
    cmd += [
        "index",
        ZS_Utility.convert_to_wsl_if_not_unix(bamFile)
    ]
    
    # Run the command
    run_samtools_index = subprocess.Popen(cmd, shell = True,
                                          stdout = subprocess.PIPE,
                                          stderr = subprocess.PIPE)
    indexout, indexerr = run_samtools_index.communicate()
    if indexout.decode("utf-8") != "" and indexerr.decode("utf-8") == "":
        print("WARNING: samtools_index_bam may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({indexout.decode("utf-8")})')
    elif indexerr.decode("utf-8") != "":
        raise Exception(("ERROR: samtools_index_bam encountered an error; have a look " +
                         f'at the stdout > {indexout.decode("utf-8")} < and stderr ' + 
                         f'> {indexerr.decode("utf-8")} < to make sense of this.'))

def bcftools_mpileup(bamFiles, fastaFile, outputFile, bcftoolsPath, regions=None):
    '''
    Parameters:
        bamFiles -- a list containing one or more strings indicating the location of the
                    BAM files to mpileup
        fastaFile -- a string indicating the location of the genome FASTA file
        outputFile -- a string indicating the location to write the mpileup result to
        bcftoolsPath -- a string indicating the location of the bcftools executable
        regions -- optionally, a list of strings formatted for bcftools like:
                    [
                        "contig1:start1-end1",
                        "contig1:start2-end2",
                        ...
                        "contig2:start1-end1",
                        "contig2:start2-end2",
                        ...
                    ]
    '''
    BAD_WORDS = ["failed", "error", "warning", "abort", "exception", "fatal", "fail", "unrecogni"]
    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(bcftoolsPath)
    cmd += [
        "mpileup", "-q", "10", "-Q", "20", "-a", "AD", "--threads", "2",
        "-f", ZS_Utility.convert_to_wsl_if_not_unix(fastaFile),
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFile)
    ]
    
    if regions != None:
        cmd += ["-r", ",".join(regions)]
    
    cmd += [ ZS_Utility.convert_to_wsl_if_not_unix(x) for x in bamFiles ]
    
    # Run the command
    run_bcftools_mpileup = subprocess.Popen(cmd, shell = True,
                                            stdout = subprocess.PIPE,
                                            stderr = subprocess.PIPE)
    mpileupout, mpileuperr = run_bcftools_mpileup.communicate()
    if mpileupout.decode("utf-8") != "" and (not any([ bw in mpileuperr.decode("utf-8").lower() for bw in BAD_WORDS ])):
        print("WARNING: bcftools_mpileup may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({mpileupout.decode("utf-8")})')
    elif any([ bw in mpileuperr.decode("utf-8").lower() for bw in BAD_WORDS ]):
        raise Exception(("ERROR: bcftools_mpileup encountered an error; have a look " +
                         f'at the stdout ({mpileupout.decode("utf-8")}) and stderr ' + 
                         f'({mpileuperr.decode("utf-8")}) to make sense of this.'))

def bcftools_call(mpileupFile, outputFile, bcftoolsPath):
    '''
    Parameters:
        mpileupFile -- a string indicating the result of bcftools mpileup
        outputFile -- a string indicating the location to write the VCF result to
        bcftoolsPath -- a string indicating the location of the bcftools executable
    '''
    BAD_WORDS = ["failed", "error", "warning", "abort", "exception", "fatal", "fail", "unrecogni"]
    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(bcftoolsPath)
    cmd += [
        "call", "-m", "-v", "-Oz", "--write-index",
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFile),
        ZS_Utility.convert_to_wsl_if_not_unix(mpileupFile)
    ]
    
    # Run the command
    run_bcftools_call = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    callout, callerr = run_bcftools_call.communicate()
    if callout.decode("utf-8") != "" and (not any([ bw in callerr.decode("utf-8").lower() for bw in BAD_WORDS ])):
        print("WARNING: bcftools_call may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({callout.decode("utf-8")})')
    elif any([ bw in callerr.decode("utf-8").lower() for bw in BAD_WORDS ]):
        raise Exception(("ERROR: bcftools_call encountered an error; have a look " +
                         f'at the stdout ({callout.decode("utf-8")}) and stderr ' + 
                         f'({callerr.decode("utf-8")}) to make sense of this.'))

def bcftools_norm_multiallelics(vcfFile, outputFile, bcftoolsPath):
    '''
    Parameters:
        vcfFile -- a string indicating the location of the VCF file to multiallele normalise
        outputFile -- a string indicating the location to write the VCF result to
        bcftoolsPath -- a string indicating the location of the bcftools executable
    '''
    BAD_WORDS = ["failed", "error", "warning", "abort", "exception", "fatal", "fail", "unrecogni"]
    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(bcftoolsPath)
    cmd += [
        "norm", "-m", "+any",
        "-Oz", "--write-index",
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFile),
        ZS_Utility.convert_to_wsl_if_not_unix(vcfFile)
    ]
    
    # Run the command
    run_bcftools_norm_m = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    normout, normerr = run_bcftools_norm_m.communicate()
    if normout.decode("utf-8") != "" and (not any([ bw in normerr.decode("utf-8") for bw in BAD_WORDS ])):
        print("WARNING: bcftools_norm_multiallelics may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({normout.decode("utf-8")})')
    elif any([ bw in normerr.decode("utf-8") for bw in BAD_WORDS ]):
        raise Exception(("ERROR: bcftools_norm_multiallelics encountered an error; have a look " +
                         f'at the stdout ({normout.decode("utf-8")}) and stderr ' + 
                         f'({normerr.decode("utf-8")}) to make sense of this.'))

def bcftools_norm_leftalign(vcfFile, fastaFile, outputFile, bcftoolsPath):
    '''
    Parameters:
        vcfFile -- a string indicating the location of the VCF file to left-align normalise
        fastaFile -- a string indicating the location of the genome FASTA file
        outputFile -- a string indicating the location to write the VCF result to
        bcftoolsPath -- a string indicating the location of the bcftools executable
    '''
    BAD_WORDS = ["failed", "error", "warning", "abort", "exception", "fatal", "fail", "unrecogni"]
    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(bcftoolsPath)
    cmd += [
        "norm", "-f", ZS_Utility.convert_to_wsl_if_not_unix(fastaFile),
        "-Oz", "--write-index",
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFile),
        ZS_Utility.convert_to_wsl_if_not_unix(vcfFile)
    ]
    
    # Run the command
    run_bcftools_norm_left = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    normout, normerr = run_bcftools_norm_left.communicate()
    if normout.decode("utf-8") != "" and (not any([ bw in normerr.decode("utf-8").lower() for bw in BAD_WORDS ])):
        print("WARNING: run_bcftools_norm_left may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({normout.decode("utf-8")})')
    elif any([ bw in normerr.decode("utf-8").lower().lower() for bw in BAD_WORDS ]):
        raise Exception(("ERROR: run_bcftools_norm_left encountered an error; have a look " +
                         f'at the stdout ({normout.decode("utf-8")}) and stderr ' + 
                         f'({normerr.decode("utf-8")}) to make sense of this.'))

def whatshap_phase(vcfFile, fastaFile, bamFiles, outputFile, whatshapPath):
    '''
    Parameters:
        vcfFile -- a string indicating the location of the VCF file to phase
        fastaFile -- a string indicating the location of the genome FASTA file
        bamFiles -- a list containing one or more strings indicating the location of the
                    BAM files to mpileup
        outputFile -- a string indicating the location to write the whatshap result to
        whatshapPath -- a string indicating the location of the whatshap executable
    '''    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(whatshapPath)
    cmd += [
        "phase", "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFile),
        f"--reference={ZS_Utility.convert_to_wsl_if_not_unix(fastaFile)}",
        ZS_Utility.convert_to_wsl_if_not_unix(vcfFile),
        *[ ZS_Utility.convert_to_wsl_if_not_unix(x) for x in bamFiles ]
    ]
    
    # Run the command
    run_whatshap_phase = subprocess.Popen(cmd, shell = True,
                                          stdout = subprocess.PIPE,
                                          stderr = subprocess.PIPE)
    phaseout, phaseerr = run_whatshap_phase.communicate()
    if phaseout.decode("utf-8") != "" and (not "Total elapsed time" in phaseerr.decode("utf-8")):
        print("WARNING: whatshap_phase may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({phaseout.decode("utf-8")})')
    elif not "Total elapsed time" in phaseerr.decode("utf-8"):
        raise Exception(("ERROR: whatshap_phase encountered an error; have a look " +
                         f'at the stdout ({phaseout.decode("utf-8")}) and stderr ' + 
                         f'({phaseerr.decode("utf-8")}) to make sense of this.'))

def whatshap_polyphase(vcfFile, fastaFile, bamFiles, outputFile, whatshapPath, ploidy=2):
    '''
    Parameters:
        vcfFile -- a string indicating the location of the VCF file to phase
        fastaFile -- a string indicating the location of the genome FASTA file
        bamFiles -- a list containing one or more strings indicating the location of the
                    BAM files to mpileup
        outputFile -- a string indicating the location to write the whatshap result to
        whatshapPath -- a string indicating the location of the whatshap executable
    '''    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(whatshapPath)
    cmd += [
        "polyphase", "--ploidy", str(ploidy), "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFile),
        f"--reference={ZS_Utility.convert_to_wsl_if_not_unix(fastaFile)}",
        ZS_Utility.convert_to_wsl_if_not_unix(vcfFile),
        *[ ZS_Utility.convert_to_wsl_if_not_unix(x) for x in bamFiles ]
    ]
    
    # Run the command
    run_whatshap_poly = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    polyout, polyerr = run_whatshap_poly.communicate()
    if polyout.decode("utf-8") != "" and (not "Total elapsed time" in polyerr.decode("utf-8")):
        print("WARNING: whatshap_polyphase may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({polyout.decode("utf-8")})')
    elif not "Total elapsed time" in polyerr.decode("utf-8"):
        raise Exception(("ERROR: whatshap_polyphase encountered an error; have a look " +
                         f'at the stdout ({polyout.decode("utf-8")}) and stderr ' + 
                         f'({polyerr.decode("utf-8")}) to make sense of this.'))

def merge_whatshap_vcfs(phaseFile, polyFile, outputFile):
    '''
    This function will merge the 'whatshap phase' and 'whatshap polyphase' VCF files
    into a single VCF file. This allows the benefit of 'polyphase' (multiallelic phasing)
    to be combined with additional strength of 'phase' for diploid phasing some variants that
    'polyphase' fails to phase.
    
    Parameters:
        phaseFile -- a string indicating the location of the 'whatshap phase' VCF file
        polyFile -- a string indicating the location of the 'whatshap polyphase' VCF file
        outputFile -- a string indicating the location to write the whatshap result to
    '''
    with open(phaseFile, "r") as phaseIn, open(polyFile, "r") as polyIn, open(outputFile, "w") as fileOut:
        # Iterate through phase file, writing header to output
        while True:
            phaseLine = phaseIn.readline()
            if phaseLine.startswith("#"):
                fileOut.write(phaseLine)
            if phaseLine.startswith("#CHROM"):
                break
        
        # Iterate through poly file, skipping over header
        while True:
            polyLine = polyIn.readline()
            if polyLine.startswith("#"):
                pass
            if polyLine.startswith("#CHROM"):
                break
        
        # Iterate line-by-line through both files, writing the best phased result to file
        while True:
            phaseLine = phaseIn.readline()
            polyLine = polyIn.readline()
            
            # Exit condition
            if phaseLine == "" or polyLine == "":
                assert phaseLine == "" and polyLine == "", \
                    "ERROR: phase and polyphase files have different numbers of lines!"
                break
            
            # Extract details from each line
            phaseSl = phaseLine.rstrip("\r\n ").split("\t")
            polySl = polyLine.rstrip("\r\n ").split("\t")
            
            # Get the GT field
            phaseGTindex = phaseSl[8].split(":").index("GT")
            polyGTindex = polySl[8].split(":").index("GT")
            
            phaseGT = phaseSl[9].split(":")[phaseGTindex]
            polyGT = polySl[9].split(":")[polyGTindex]
            
            # Write the best phased line to the output
            if "|" in phaseGT: # 'whatshap phase' is assumed to be better if it phased something
                fileOut.write(phaseLine)
            elif "|" in polyGT: # 'whatshap polyphase', if phase failed, is assumed to be better than nothing
                fileOut.write(polyLine)
            else: # If both failed to phase, write the 'whatshap phase' line
                fileOut.write(phaseLine)

def bgzip_file(fileName, bgzipPath):
    '''
    Parameters:
        fileName -- a string indicating the location of a file to index
        bgzipPath -- a string indicating the location of the bgzip executable
    '''
    BAD_WORDS = ["failed", "error", "warning", "abort", "exception", "fatal", "fail", "unrecogni"]
    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(bgzipPath)
    cmd.append(ZS_Utility.convert_to_wsl_if_not_unix(fileName))
    
    # Run the command
    run_bgzip = subprocess.Popen(cmd, shell = True,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE)
    bgzipout, bgziperr = run_bgzip.communicate()
    if bgzipout.decode("utf-8") != "" and (not any([ bw in bgziperr.decode("utf-8").lower() for bw in BAD_WORDS ])):
        print("WARNING: bgzip_file may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({bgzipout.decode("utf-8")})')
    elif any([ bw in bgziperr.decode("utf-8").lower() for bw in BAD_WORDS ]):
        raise Exception(("ERROR: bgzip_file encountered an error; have a look " +
                         f'at the stdout ({bgzipout.decode("utf-8")}) and stderr ' + 
                         f'({bgziperr.decode("utf-8")}) to make sense of this.'))

def tabix_file(fileName, tabixPath):
    '''
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file to index
        tabixPath -- a string indicating the location of the tabix executable
    '''
    BAD_WORDS = ["failed", "error", "warning", "abort", "exception", "fatal", "fail", "unrecogni"]
    
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(tabixPath)
    cmd.append(ZS_Utility.convert_to_wsl_if_not_unix(fileName))
    
    # Run the command
    run_tabix = subprocess.Popen(cmd, shell = True,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE)
    tabixout, tabixerr = run_tabix.communicate()
    if tabixout.decode("utf-8") != "" and (not any([ bw in tabixerr.decode("utf-8").lower() for bw in BAD_WORDS ])):
        print("WARNING: tabix_file may have encountered an error, since the stdout is not empty as expected. " +
              f'Please check the stdout for more information ({tabixout.decode("utf-8")})')
    elif any([ bw in tabixerr.decode("utf-8").lower() for bw in BAD_WORDS ]):
        raise Exception(("ERROR: tabix_file encountered an error; have a look " +
                         f'at the stdout ({tabixout.decode("utf-8")}) and stderr ' + 
                         f'({tabixerr.decode("utf-8")}) to make sense of this.'))

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3 and a directory containing one or more BAM files.
    It will subset the BAM files to the parent feature regions contained in the GFF3 file,
    so you should filter the GFF3 file to only contain the features you're interested in.
    Reads aligning to parent features +/- the --plusMinus distance will be mpileup'd into
    a combined VCF file. Note that this code will operate in WSL where possible if you are
    running on Windows. Also note that picard is expected to be "executable-ised" e.g., like
    with a file in your local/bin called picard containing 'java -jar /mnt/c/bio/picard/picard.jar $@'.
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
    p.add_argument("--samtools", dest="samtools",
                   required=False,
                   help="""Optionally, specify the samtools executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--bcftools", dest="bcftools",
                   required=False,
                   help="""Optionally, specify the bcftools executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--picard", dest="picard",
                   required=False,
                   help="""Optionally, specify the picard executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--whatshap", dest="whatshap",
                   required=False,
                   help="""Optionally, specify the whatshap executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--bgzip", dest="bgzip",
                   required=False,
                   help="""Optionally, specify the bgzip executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--tabix", dest="tabix",
                   required=False,
                   help="""Optionally, specify the tabix executable file
                   if it is not discoverable in the path""",
                   default=None)
    # Opts (behavioural)
    p.add_argument("--plusMinus", dest="plusMinus",
                   type=int,
                   required=False,
                   help="""Optionally, specify how far to the left and right
                   of a gene region to subset BAMs to; default == 0""",
                   default=0)
    p.add_argument("--bamSuffix", dest="bamSuffix",
                   required=False,
                   help="""Indicate the suffix of the BAM files to look for;
                   default == '.sorted.bam'""",
                   default=".sorted.bam")
    p.add_argument("--setRG", dest="setRG",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should re-set read groups for each BAM;
                   will use the BAM file name sans bamSuffix as the read group ID.""",
                   default=False)
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    # Opts (read groups)
    p.add_argument("--rgLB", dest="rgLB",
                   required=False,
                   help="""If you are setting read groups, optionally specify the
                   library name to use; default == 'lib1'""",
                   default="lib1")
    p.add_argument("--rgPL", dest="rgPL",
                   required=False,
                   help="""If you are setting read groups, optionally specify the
                   platform name to use; default == 'ILLUMINA'""",
                   default="ILLUMINA")
    p.add_argument("--rgPU", dest="rgPU",
                   required=False,
                   help="""If you are setting read groups, optionally specify the
                   unit name to use; default == 'unit1'""",
                   default="unit1")
    
    args = p.parse_args()
    validate_args(args)
    
    # Index the FASTA file
    if not os.path.exists(f"{args.fastaFile}.fai"):
        samtools_faidx(args.fastaFile, args.samtools)
    
    # Parse GFF3
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse = not args.relaxedParsing)
    
    # Locate BAM files
    bamFiles = [
        os.path.join(args.bamDirectory, file)
        for file in os.listdir(args.bamDirectory)
        if file.endswith(args.bamSuffix)
    ]
    
    # Locate regions to filter BAMs to
    coordinatesDict = gff3_to_coordinatesDict(gff3Obj, args.plusMinus)
    
    # Merge overlapping ranges
    coordinatesDict = merge_overlapping_coordinates(coordinatesDict)
    
    # Format regions for samtools
    samtoolsRegions = format_samtools_regions(coordinatesDict)
    
    # Subset BAMs using samtools
    subsetBamFiles = []
    for bamFile in bamFiles:
        outputFile = os.path.join(args.outputDirectory, os.path.basename(bamFile))
        if not os.path.exists(outputFile):
            samtools_subset_bam(bamFile, outputFile, samtoolsRegions, args.samtools)
        else:
            print(f"Subset BAM file '{outputFile}' already exists; skipping.")
        subsetBamFiles.append(outputFile)
    
    # Update read groups if necessary
    if args.setRG:
        if not os.path.exists(os.path.join(args.outputDirectory, "readgroups_were_set.flag")):
            for bamFile in subsetBamFiles:
                assert os.path.exists(bamFile), f"Subset BAM file '{bamFile}' not found; this is unexpected."
                
                try:
                    # Move the original file
                    tmpFile = f"{bamFile}.tmp"
                    shutil.move(bamFile, tmpFile)
                    
                    # Figure out what our read group ID/sample name should be
                    rgID = os.path.basename(bamFile).replace(args.bamSuffix, "")
                    rgSM = rgID
                    
                    # Update read group using picard
                    picard_addreplace_readgroups(tmpFile, bamFile,
                                                rgID, args.rgLB, args.rgPL,
                                                args.rgPU, rgSM, args.picard)
                    
                    # Clean up tmp file
                    os.remove(tmpFile)
                except Exception as e:
                    # Clean up
                    if os.path.exists(tmpFile):
                        if not os.path.exists(bamFile):
                            shutil.move(tmpFile, bamFile)
                        else:
                            os.remove(tmpFile)
                    
                    # Propagate error
                    raise e
            
            # Set a flag to indicate that we've updated read groups
            open(os.path.join(args.outputDirectory, "readgroups_were_set.flag"), "w").close()
        else:
            print("Read groups have already been set for these BAM files; skipping.")
    
    # Index BAMs
    for bamFile in subsetBamFiles:
        if not os.path.exists(f"{bamFile}.bai"):
            samtools_index_bam(bamFile, args.samtools)
        else:
            print(f"Indexed BAM file '{outputFile}' already exists; skipping.")
    
    # Run bcftools mpileup
    mpileupFileName = os.path.join(args.outputDirectory, "bcftools.mpileup")
    if not os.path.exists(mpileupFileName) or not \
        os.path.exists(os.path.join(args.outputDirectory, "mpileup_was_successful.flag")):
            bcftools_mpileup(subsetBamFiles, args.fastaFile, mpileupFileName, args.bcftools, samtoolsRegions)
            open(os.path.join(args.outputDirectory, "mpileup_was_successful.flag"), "w").close()
    else:
        print(f"bcftools mpileup file has already been generated; skipping.")
    
    # Run bcftools call
    vcfFileName = os.path.join(args.outputDirectory, "bcftools.vcf.gz")
    if not os.path.exists(vcfFileName) or not \
        os.path.exists(os.path.join(args.outputDirectory, "call_was_successful.flag")):
            bcftools_call(mpileupFileName, vcfFileName, args.bcftools)
            open(os.path.join(args.outputDirectory, "call_was_successful.flag"), "w").close()
    else:
        print(f"bcftools call file has already been generated; skipping.")
    
    # Normalise multiallelic variants
    multialleleFileName = os.path.join(args.outputDirectory, "bcftools.multiallele.norm.vcf.gz")
    if not os.path.exists(vcfFileName) or not \
        os.path.exists(os.path.join(args.outputDirectory, "norm_multiallele_was_successful.flag")):            
            bcftools_norm_multiallelics(vcfFileName, multialleleFileName, args.bcftools)
            open(os.path.join(args.outputDirectory, "norm_multiallele_was_successful.flag"), "w").close()
    else:
        print(f"bcftools norm -m file has already been generated; skipping.")
    
    # Left-align and produce finalised variants
    finalVcfFileName = os.path.join(args.outputDirectory, "bcftools.final.vcf.gz")
    if not os.path.exists(vcfFileName) or not \
        os.path.exists(os.path.join(args.outputDirectory, "finalise_variants_was_successful.flag")):            
            bcftools_norm_leftalign(multialleleFileName, args.fastaFile, finalVcfFileName, args.bcftools)
            open(os.path.join(args.outputDirectory, "finalise_variants_was_successful.flag"), "w").close()
    else:
        print(f"bcftools norm -f file has already been generated; skipping.")
    
    # Run WhatsHap phase
    mergedWhatsHapName = os.path.join(args.outputDirectory, "whatshap.merged.vcf")
    compressedWhatsHapName = f"{mergedWhatsHapName}.gz"
    
    phasedFileName = os.path.join(args.outputDirectory, "whatshap.phase.vcf")
    if (not os.path.exists(phasedFileName) and not os.path.exists(compressedWhatsHapName)) and not \
        os.path.exists(os.path.join(args.outputDirectory, "phase_was_successful.flag")):
            whatshap_phase(finalVcfFileName, args.fastaFile, subsetBamFiles, phasedFileName, args.whatshap)
            open(os.path.join(args.outputDirectory, "phase_was_successful.flag"), "w").close()
    else:
        print(f"whatshap phase file has already been generated; skipping.")
    
    # Run WhatsHap polyphase
    polyFileName = os.path.join(args.outputDirectory, "whatshap.poly.vcf")
    if (not os.path.exists(polyFileName) and not os.path.exists(compressedWhatsHapName)) and not \
        os.path.exists(os.path.join(args.outputDirectory, "poly_was_successful.flag")):
            whatshap_polyphase(finalVcfFileName, args.fastaFile, subsetBamFiles, polyFileName, args.whatshap)
            open(os.path.join(args.outputDirectory, "poly_was_successful.flag"), "w").close()
    else:
        print(f"whatshap polyphase file has already been generated; skipping.")
    
    # Merge phase and polyphase files
    if (not os.path.exists(mergedWhatsHapName) and not os.path.exists(compressedWhatsHapName)) and not \
        os.path.exists(os.path.join(args.outputDirectory, "merge_was_successful.flag")):
            merge_whatshap_vcfs(phasedFileName, polyFileName, mergedWhatsHapName)
            open(os.path.join(args.outputDirectory, "merge_was_successful.flag"), "w").close()
    else:
        print(f"whatshap phase+polyphase merge file has already been generated; skipping.")
    
    # Compress the phased VCF
    if not os.path.exists(compressedWhatsHapName) or not \
        os.path.exists(os.path.join(args.outputDirectory, "compression_was_successful.flag")):
            bgzip_file(mergedWhatsHapName, args.bgzip)
            open(os.path.join(args.outputDirectory, "compression_was_successful.flag"), "w").close()
    else:
        print(f"whatshap VCF has already been compressed; skipping.")
    
    # Index the phased VCF
    compressedIndexName = f"{compressedWhatsHapName}.tbi"
    if not os.path.exists(compressedIndexName) or not \
        os.path.exists(os.path.join(args.outputDirectory, "tabix_was_successful.flag")):
            tabix_file(compressedWhatsHapName, args.tabix)
            open(os.path.join(args.outputDirectory, "tabix_was_successful.flag"), "w").close()
    else:
        print(f"whatshap VCF has already been tabix indexed; skipping.")
    
    # Parse the VCF
    phasedVCF_obj = ZS_VCFIO.PhasedVCF(compressedWhatsHapName)
    phasedVCF_obj.parse_whatshap_vcf()
    
    # Generate phased gene sequences for each sample
    FASTA_obj = Fasta(args.fastaFile)
    
    for haplotypeNum in range(1, 3):
        haplotypeFastaFile = os.path.join(args.outputDirectory, f"haplotype_{haplotypeNum}.fasta")
        # Skip if we've generated this file already
        if not os.path.exists(haplotypeFastaFile) or not \
            os.path.exists(os.path.join(args.outputDirectory, f"fasta_phasing_{haplotypeNum}_was_successful.flag")):
                # Start writing phased sequences for this haplotype
                with open(haplotypeFastaFile, "w") as fileOut:
                    for parentType in gff3Obj.parentTypes:
                        
                        # Iterate through parent features
                        for parentFeature in gff3Obj.types[parentType]:
                            mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(parentFeature)
                            
                            # Skip if there are no SNPs on this contig
                            if not parentFeature.contig in phasedVCF_obj:
                                fileOut.write(f"# {parentFeature.ID} is reference type\n")
                                continue
                            
                            # Get the exon regions for this gene
                            exonCoords = [ exonFeature.coords for exonFeature in mrnaFeature.exon ]
                            
                            # Find any SNPs located within this gene
                            positionsInGene = [
                                pos
                                for pos in phasedVCF_obj[parentFeature.contig]
                                for start, end in exonCoords
                                if start <= pos <= end
                            ]
                            
                            # Skip if there are no SNPs in this gene
                            if len(positionsInGene) == 0:
                                fileOut.write(f"# {parentFeature.ID} is reference type\n")
                                continue
                            
                            # Get the mRNA feature sequence
                            exon_FastASeq_obj, exon_featureType, exon_startingFrame = \
                                gff3Obj.retrieve_sequence_from_FASTA(FASTA_obj, mrnaFeature.ID, "exon")

                            #if parentFeature.start < 4743973 < parentFeature.end:
                            #    stophere
                            
                            # Iterate through samples and generate haplotypes
                            sampleHaplotypes = {}
                            for sampleID in phasedVCF_obj.samples:
                                sampleHaplotypes[sampleID] = []
                                
                                # Get any SNPs located within this gene for this sample
                                sampleSNPs = {
                                    pos: phasedVCF_obj[parentFeature.contig][pos][sampleID]
                                    for pos in positionsInGene
                                }
                                
                                # Extract variants within exon regions, and modify positions to be relative to the exon
                                sampleSNPs = convert_vcf_snps_to_cds_snps(mrnaFeature, sampleSNPs, featureType="exon")
                                
                                # Get the ordered SNP positions
                                orderedPositions = list(sampleSNPs.keys())
                                orderedPositions.sort()
                                
                                # Format haplotype
                                haplotypes = [[], []]
                                for pos in orderedPositions:
                                    genotypeDict = sampleSNPs[pos]
                                    ## TBD: Figure out how to link phase groups using PS tag and such
                                    haplotypes[0].append([genotypeDict["ref_alt"][0], genotypeDict["GT"][0]])
                                    haplotypes[1].append([genotypeDict["ref_alt"][0], genotypeDict["GT"][1]])
                            
                            # Get just the unique haplotypes so we can avoid redundant sequence extraction
                            # uniqueHaplotypes = list(set(
                            #     tuple(subgenotype)
                            #     for genotype in haploCodeDict.values()
                            #     for subgenotype in genotype
                            # ))
                            
                            # Extract haplotype sequences for all unique haplotypes
                            # haplotypeSequences = []
                            # for haplotype in uniqueHaplotypes:
                            #     haplotypeSequence = edit_reference_to_haplotype_sequence(cds_FastASeq_obj.seq[:], haplotype, orderedPositions, newSnpDict, mrnaFeature.strand)
                            #     haplotypeSequences.append(haplotypeSequence)
                            
                            # Write haplotype sequences to file
                            ## TBD...
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
