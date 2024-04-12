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
from Function_packages import ZS_GFF3IO, ZS_Utility, ZS_VCFIO, ZS_SeqIO

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
    tree.merge_overlaps(strict=False)
    
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

def parse_whatshap_for_merge(vcfFile):
    '''
    This function will parse a phased WhatsHap result file, produced using either
    'phase' or 'polyphase', into a data structure that can be used to merge the
    two VCF files.
    
    Parameters:
        vcfFile -- a string indicating the location of a WhatsHap phased VCF file
    Returns:
        phaseDict -- a dictionary with structure like:
                    {
                        "contig1": {
                            position1: {
                                "GT": "0|1",
                                "PS": "1"
                            },
                            position2: {
                                "GT": "1|0",
                                "PS": "1"
                            },
                            ...
                        },
                        "contig2": {
                            position1: {
                                "GT": "0|1",
                                "PS": "1"
                            },
                            position2: {
                                "GT": "1|0",
                                "PS": "1"
                            },
                            ...
                        },
                        ...
                    }
    '''
    phaseDict = {}
    with open(vcfFile, "r") as phaseIn:
        for line in phaseIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = sl[9:]
            
            # Handle content lines
            else:
                chrom, pos, id, ref, alt, \
                    qual, filter, info, _format = sl[0:9]
                pos = int(pos)
                formatList = _format.split(":")
                
                GTindex = formatList.index("GT")
                PSindex = formatList.index("PS") if "PS" in formatList else -1
                
                # Set up phaseDict structure
                phaseDict.setdefault(chrom, {})
                phaseDict[chrom].setdefault(pos, {})
                
                # Extract GT and PS fields per sample
                for i, sample in enumerate(samples):
                    sampleDetails = sl[9 + i].split(":")
                    phaseDict[chrom][pos][sample] = {
                        "GT": sampleDetails[GTindex],
                        "PS": sampleDetails[PSindex] if PSindex != -1 else None
                    }
    return phaseDict

def get_positions_from_list(orderedList, start, end):
    '''
    Parameters:
        orderedList -- a list of integers in ascending order
        start -- the start of the range to search for
        end -- the end of the range to search for
    Returns:
        a list of integers in the range [start, end] that are in the orderedList
    '''
    return [ x for x in orderedList if x >= start and x <= end ]

def find_proximal_snps(pos, chrom, sample, phaseDict, polyDict, contigPositions, psBoundaries):
    '''
    This function will find the most proximal phased SNPs before and after a given position
    in a given sample.
    
    Parameters:
        pos -- an integer indicating the position to search around
        chrom -- a string indicating the contig to search on
        sample -- a string indicating the sample to search in
        phaseDict -- a dictionary resulting from parse_whatshap_for_merge run on the 'phase' VCF
        polyDict -- a dictionary resulting from parse_whatshap_for_merge run on the 'polyphase' VCF
        contigPositions -- a dictionary of lists of integers, indicating the variant positions on each contig
        psBoundaries -- a dictionary of lists of lists, indicating the start and end boundaries of each
                        PS group, indexed by the sample and contig
    Returns:
        prevPhase -- a string indicating the phased GT of the SNP before the given position
        prevPoly -- a string indicating the polyphase GT of the SNP before the given position
        nextPhase -- a string indicating the phased GT of the SNP after the given position
        nextPoly -- a string indicating the polyphase GT of the SNP after the given position
    '''
    boundariesList = psBoundaries[chrom][sample]
    psStart, psEnd = [ x for x in boundariesList if pos >= x[0] and pos <= x[1] ][0]
    
    psPositions = get_positions_from_list(contigPositions[chrom], psStart, psEnd)
    currentPosIndex = psPositions.index(pos)
    
    found = 0
    prevPhase, prevPoly = None, None
    for i in range(currentPosIndex - 1, -1, -1):
        prevPos = psPositions[i]
        phaseGT = phaseDict[chrom][prevPos][sample]["GT"]
        polyGT = polyDict[chrom][prevPos][sample]["GT"]
        
        if "|" in phaseGT and prevPhase == None:
            prevPhase = phaseGT
            found += 1
        if "|" in polyGT and prevPoly == None:
            prevPoly = polyGT
            found += 1
        if found == 2:
            break
    
    found = 0
    nextPhase, nextPoly = None, None
    for i in range(currentPosIndex + 1, len(psPositions)):
        nextPos = psPositions[i]
        phaseGT = phaseDict[chrom][nextPos][sample]["GT"]
        polyGT = polyDict[chrom][nextPos][sample]["GT"]
        
        if "|" in phaseGT and nextPhase == None:
            nextPhase = phaseGT
            found += 1
        if "|" in polyGT and nextPoly == None:
            nextPoly = polyGT
            found += 1
        if found == 2:
            break
    
    return prevPhase, prevPoly, nextPhase, nextPoly

def merge_whatshap_vcfs(vcfFile, phaseDict, polyDict, outputFile):
    '''
    This function will merge the 'whatshap phase' and 'whatshap polyphase' VCF files
    into a single VCF file. This allows the benefit of 'polyphase' (multiallelic phasing)
    to be combined with additional strength of 'phase' for diploid phasing some variants that
    'polyphase' fails to phase. In order to do this, you must first parse the two VCFs
    with parse_whatshap_for_merge() and provide their results as input.
    
    Parameters:
        vcfFile -- a string indicating the location of a WhatsHap phased VCF file
        phaseDict -- a dictionary resulting from parse_whatshap_for_merge run on the 'phase' VCF
        polyDict -- a dictionary resulting from parse_whatshap_for_merge run on the 'polyphase' VCF
        outputFile -- a string indicating the location to write the merged whatshap result to
    '''
    assert set(phaseDict.keys()) == set(polyDict.keys()), "Contigs do not match in phase and polyphase dicts"
    for chrom in phaseDict.keys():
        assert set(phaseDict[chrom].keys()) == set(polyDict[chrom].keys()), \
            "Positions do not match in phase and polyphase dicts"
        for pos in phaseDict[chrom].keys():
            assert set(phaseDict[chrom][pos].keys()) == set(polyDict[chrom][pos].keys()), \
                "Samples do not match in phase and polyphase dicts"
    
    # Establish data structures for easy searching of coordinates
    contigPositions = {
        chrom: sorted(phaseDict[chrom].keys())
        for chrom in phaseDict.keys()
    }
    
    # Find maximal boundaries of each PS group
    psCoords = {}
    for chrom, posList in contigPositions.items():
        psCoords.setdefault(chrom, {})
        for pos in posList:
            phaseSampleDict = phaseDict[chrom][pos]
            for sampleID, phaseValueDict in phaseSampleDict.items():
                psCoords[chrom].setdefault(sampleID, {})
                polyValueDict = polyDict[chrom][pos][sampleID]
                
                if phaseValueDict["PS"] != None:
                    psCoords[chrom][sampleID].setdefault(phaseValueDict["PS"], [pos, pos]) # init start and end as same for PS group
                    psCoords[chrom][sampleID][phaseValueDict["PS"]][1] = pos # update the end of the PS group with most recent pos
                if polyValueDict["PS"] != None:
                    psCoords[chrom][sampleID].setdefault(polyValueDict["PS"], [pos, pos])
                    psCoords[chrom][sampleID][polyValueDict["PS"]][1] = pos
    
    # Collapse PS group boundaries across phase and poly
    psBoundaries = {}
    for chrom, sampleDict in psCoords.items():
        psBoundaries.setdefault(chrom, {})
        for sampleID, psDict in sampleDict.items():
            thisCoords = [ x for x in psDict.values() if x[0] != x[1] ] # single position PS groups don't make sense and cause interval tree errors
            psBoundaries[chrom][sampleID] = merge_coordinates_list(thisCoords)
    
    # Produce the merged VCF
    with open(vcfFile, "r") as vcfIn, open(outputFile, "w") as vcfOut:
        for line in vcfIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Write the header
            if line.startswith("#CHROM"):
                samples = sl[9:]
            if line.startswith("#"):
                vcfOut.write(line)
                continue
            
            # Extract details from the body
            chrom, pos, id, ref, alt, \
                qual, filter, info, format = sl[0:9]
            pos = int(pos)
            formatList = format.split(":")
            
            gtIndex = formatList.index("GT")
            psIndex = formatList.index("PS") if "PS" in formatList else -1
            
            # Iterate through the samples
            rowSampleData = sl[9:]
            phasedSomething = False
            for i, sample in enumerate(samples):
                # Extract details of the sample
                sampleData = rowSampleData[i].split(":")
                phaseGT = sampleData[gtIndex]
                
                # Update PS group for all phased lines
                "Our phase set start site might change by merging 'polyphase' results in"
                if "|" in phaseGT:
                    if chrom in psBoundaries:
                        for start, end in psBoundaries[chrom][sample]:
                            if pos >= start and pos <= end:
                                sampleData[psIndex] = str(start)
                                break
                
                # Handle unphased 'phase' lines
                if not "|" in phaseGT:
                    # Get the 'polyphase' data
                    polyData = polyDict[chrom][pos][sample]
                    polyGT = polyData["GT"]
                    
                    # Skip if the 'polyphase' line is also unphased
                    if not "|" in polyGT:
                        continue
                    
                    # See if we need to flip the 'polyphase' GT
                    prevPhase, prevPoly, nextPhase, nextPoly = find_proximal_snps(pos, chrom, sample, \
                        phaseDict, polyDict, contigPositions, psBoundaries)
                    
                    flip = 0
                    dontflip = 0
                    if prevPhase != None and prevPoly != None:
                        if prevPhase != prevPoly:
                            flip += 1
                        else:
                            dontflip += 1
                    if nextPhase != None and nextPoly != None:
                        if nextPhase != nextPoly:
                            flip += 1
                        else:
                            dontflip += 1
                    
                    if flip > dontflip:
                        polyGT = "|".join(polyGT.split("|")[::-1])
                    elif dontflip > flip:
                        pass
                    else:
                        """If we're here, we have a tie in the flip/don't flip decision which means the phasing is
                        unreliable and there's no benefit to merging in 'polyphase's data. We'll skip this line."""
                        continue
                    
                    # Specifically set PS group if we're merging a 'polyphase' line in where the original was unphased
                    if chrom in psBoundaries:
                        for start, end in psBoundaries[chrom][sample]:
                            if pos >= start and pos <= end:
                                sampleData.append(str(start)) # add PS group to the sampleData list
                                break
                    
                    # Update the sample details
                    sampleData[gtIndex] = polyGT
                    phasedSomething = True
                
                # Update the sample in the line
                rowSampleData[i] = ":".join(sampleData)
            
            # Fix the line up if PS did not exist in 'poly'
            if psIndex == -1 and phasedSomething:
                # Add a blank PS if we didn't phase this sample
                for i, sample in enumerate(samples):
                    if rowSampleData[i].count(":") == len(formatList) - 1:
                        rowSampleData[i] += ":."
                
                # Add PS to the format list
                formatList.append("PS")
            
            # Write the line to the output
            newLine = "\t".join(sl[:8] + [":".join(formatList)] + rowSampleData) + "\n"
            vcfOut.write(newLine)

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

def edit_reference_to_haplotype_sequence(referenceSeq, haplotypeEditList,
                                         strand, VALIDATE_STRICT=True): # turn it to False when testing/dev is done
    '''
    This function will take in a reference nucleotide sequence, typically representing
    a CDS for a gene in either +ve or -ve strand, and generates the haplotype version
    of that sequence.
    
    Parameters:
        referenceSeq -- the sequence as a string prior to any editing
        haplotypeEditList -- a list with structure like:
                             [
                                 [pos1, "ref", "allele"],
                                 [pos2, "ref", "allele"],
                                 ...
                             ]; at each position, the reference allele is replaced with the
                             allele allele
        strand -- a string equal to "+" or "-" indicating the strandedness of this gene
        VALIDATE_STRICT -- a boolean indicating whether sequence validation should occur;
                           this is useful during testing and development but it causes issues
                           for things that are very hard to validate (i.e., the changes this
                           code make are correct, but nested variants can mess with the
                           validation).
    Returns:
        editedSeq -- an edited version of the input referenceSeq with all variations made
    '''
    # referenceSeq, haplotypeEditList, strand = exon_FastASeq_obj.seq[:], haplotypeExemplar, mrnaFeature.strand
    # VALIDATE_STRICT=True
    
    # Make sure haplotypeEditList is sorted by position
    assert [ x[0] for x in haplotypeEditList ] == sorted([ x[0] for x in haplotypeEditList ]), \
        "haplotypeEditList must be sorted by position!"
    
    # Perform the editing operation
    for i in range(len(haplotypeEditList)-1, -1, -1): # iterate backwards through positions and variants
        pos, ref, variant = haplotypeEditList[i]
        
        # Skip if variant is reference type
        if ref == variant:
            continue
        
        # Validate that our position is correct and edit the sequence
        if VALIDATE_STRICT:
            # Prevent errors with nested alleles
            if i != len(haplotypeEditList)-1: # if it's not the last variant, meaning we could've edited part of this variant
                if pos + len(ref) < haplotypeEditList[i+1][0]: # if the next variant doesn't overlap this one, meaning it's not nested
                    assert referenceSeq[pos:pos+len(ref)].upper() == ref.upper(), \
                        "Zac, you need to fix your +ve haplotype positioning code!"
        referenceSeq = referenceSeq[:pos] + variant + referenceSeq[pos+len(ref):]
    
    # Reverse complement the sequence if it's on the -ve strand
    editedSeq = referenceSeq if strand == "+" else ZS_SeqIO.FastASeq.get_reverse_complement(self=None, staticSeq=referenceSeq)
    return editedSeq.upper() # we'd like all sequences to be upper cased

def localise_vcf_snps_to_feature(mrnaFeature, snpDict, featureType="CDS"):
    '''
    Receives a mRNA feature and a dictionary indicating SNP locations as interpreted
    from a VCF, and alters the positions to point to locations in the CDS where edits
    should be made. This function handles +ve and -ve stranded mRNA features differenly
    to give the appropriate coordinates in a 5' -> 3' reading direction.
    
    Parameters:
        mrnaFeature -- a ZS_GFF3IO.Feature object representing a mRNA
        snpDict -- a dictionary with structure like:
                   {
                       pos1: { ... } # contents of dictionary don't matter
                   }
        embedOriginalPos -- optional; boolean to indicate whether the output dict
                            should also index the original SNP index. Provided as
                            optional since I'm adding this functionality into
                            legacy code and I don't want to break stuff.
        featureType -- a string indicating the type of feature we're localising to
                       e.g., a child feature of the mRNA like "CDS" or "exon".
    '''
    #mrnaFeature, snpDict, featureType = mrnaFeature, sampleSNPs, "exon"
    
    assert hasattr(mrnaFeature, featureType), \
        f"ERROR: mrnaFeature '{mrnaFeature.ID}' does not have a {featureType} attribute!"
    
    # Get the order of snp positions to iterate through
    "Ordering is necessary for the ongoingCount value to be correct"
    orderedPositions = sorted(snpDict.keys(), reverse = True)
    
    newSnpDict = {}
    ongoingCount = 0
    for childFeature in sorted(mrnaFeature.__dict__[featureType], key = lambda x: x.start):        
        # Check each position to see if we need to localise it
        for pos in orderedPositions:
            # If the position overlaps this CDS section
            if pos >= childFeature.start and pos <= childFeature.end:
                genotypeDict = snpDict[pos]
                # Get the adjusted position
                newPos = pos - childFeature.start + ongoingCount
                
                # Handle splice site variants
                refAllele = genotypeDict["ref_alt"][0]
                skipThisPos = False
                if (pos + len(refAllele) - 1) > childFeature.end:
                    # Modify the ref allele to be contained within the exon
                    allowedAlleleLength = childFeature.end - pos + 1
                    newRefAllele = refAllele[0:allowedAlleleLength]
                    # Modify alt allele(s)
                    newAltAlleles = [allele[0:allowedAlleleLength] for allele in genotypeDict["ref_alt"][1:]]
                    # If an alt allele is identical to our reference, eliminate it now
                    deleteIndices = []
                    for x in range(len(newAltAlleles)):
                        if newAltAlleles[x] == newRefAllele:
                            deleteIndices.append(x)
                            for sampleID, genotype in genotypeDict.items():
                                if sampleID != "ref_alt":
                                    for z in range(len(genotype)):
                                        if genotype[z] == x+1: # x+1 gives the index of our alt allele in ["ref_alt"]
                                            genotype[z] = 0 # set it to the ref allele index
                                    genotypeDict[sampleID] = genotype
                    for index in deleteIndices[::-1]:
                        del newAltAlleles[index] # remove it from our alt alleles values
                    # If we no longer have any variants contained within the CDS region, eliminate this variant position
                    if newAltAlleles == []:
                        skipThisPos = True
                    # Otherwise, update the ref_alt allele in our dictionary
                    else:
                        genotypeDict["ref_alt"] = [newRefAllele, *newAltAlleles]
                
                # Handle normal scenarios / index the modified alleles if relevant
                if skipThisPos is False:
                    newSnpDict[newPos] = genotypeDict
        ongoingCount += childFeature.end - childFeature.start + 1 # feature coords are 1-based inclusive, so 1->1 is a valid coord
    return newSnpDict

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
            phaseDict = parse_whatshap_for_merge(phasedFileName)
            polyDict = parse_whatshap_for_merge(polyFileName)
            
            merge_whatshap_vcfs(phasedFileName, phaseDict, polyDict, mergedWhatsHapName)
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
    
    # Skip if we've completed this already
    if os.path.exists(os.path.join(args.outputDirectory, f"fasta_phasing_was_successful.flag")):
        print("Phased FASTA files have already been generated; skipping.")
    
    # Otherwise, generate the phased FASTA files for each gene
    else:
        # Make the output directory
        phasingOutputDir = os.path.join(args.outputDirectory, "phasing")
        os.makedirs(phasingOutputDir, exist_ok=True)
        
        # Perform main phasing loop
        with open(os.path.join(args.outputDirectory, "haplotyping_results.tsv"), "w") as logOut:
            # Write log file header
            logOut.write("GeneID\thas_snps\t{0}\n".format("\t".join( [ f"{sid}_haplotype1\t{sid}_haplotype2" for sid in phasedVCF_obj.samples ])))
                
            # Iterate through genes
            for parentType in gff3Obj.parentTypes:
                for parentFeature in gff3Obj.types[parentType]:
                    mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(parentFeature)
                    
                    # Get the output file name and skip if it exists
                    fastaFileName = os.path.join(phasingOutputDir, f"{parentFeature.ID}.fasta")
                    if os.path.isfile(fastaFileName):
                        print(f"'{parentFeature.ID}' FASTA file has already been generated; skipping.")
                        continue
                    
                    # Skip if there are no SNPs on this contig
                    if not parentFeature.contig in phasedVCF_obj:
                        logOut.write("{0}\tno\t{1}\n".format(parentFeature.ID, "\t".join(["."] * len(phasedVCF_obj.samples) * 2)))
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
                        logOut.write("{0}\tno\t{1}\n".format(parentFeature.ID, "\t".join(["."] * len(phasedVCF_obj.samples) * 2)))
                        continue
                    
                    # Start writing the FASTA file
                    with open(fastaFileName, "w") as fileOut:
                        # Get the mRNA feature sequence
                        exon_FastASeq_obj, exon_featureType, exon_startingFrame = \
                            gff3Obj.retrieve_sequence_from_FASTA(FASTA_obj, mrnaFeature.ID, "exon", skipComplement=True)
                        
                        # Iterate through samples and generate haplotype codes
                        sampleHaplotypes = {}
                        sampleSwitchPoints = {sampleID: False for sampleID in phasedVCF_obj.samples}
                        for sampleID in phasedVCF_obj.samples:
                            sampleHaplotypes[sampleID] = []
                            
                            # Get any SNPs located within this gene for this sample
                            sampleSNPs = {
                                pos: phasedVCF_obj[parentFeature.contig][pos][sampleID]
                                for pos in positionsInGene
                            }
                            
                            # Extract variants within exon regions, and modify positions to be relative to the exon
                            sampleSNPs = localise_vcf_snps_to_feature(mrnaFeature, sampleSNPs, featureType="exon")
                            
                            # Get the ordered SNP positions
                            orderedPositions = list(sampleSNPs.keys())
                            orderedPositions.sort()
                            
                            # Format haplotype noting if there is a potential switch point
                            hap1, hap2 = [], [] # contains [pos, ref, allele]
                            for i, pos in enumerate(orderedPositions):
                                genotypeDict = sampleSNPs[pos]
                                
                                # Determine if this is a potential switch point
                                if i != 0:
                                    if not genotypeDict["phased"]:
                                        sampleSwitchPoints[sampleID] = True
                                    elif "PS" in genotypeDict and lastPS != None and genotypeDict["PS"] != lastPS:
                                        sampleSwitchPoints[sampleID] = True
                                
                                # Store the haplotype
                                hap1.append([pos, genotypeDict["ref_alt"][0], genotypeDict["GT"][0]])
                                hap2.append([pos, genotypeDict["ref_alt"][0], genotypeDict["GT"][1]])
                                
                                lastPS = genotypeDict["PS"] if "PS" in genotypeDict else None

                            # Store the haplotypes for this sample
                            sampleHaplotypes[sampleID] = [hap1, hap2]
                        
                        # Get just the unique haplotypes so we can avoid redundant sequence extraction
                        "Done for efficiency purposes"
                        sampleCodes = {
                            sampleID : [
                                tuple( x[-1] for x in haplotypeList[0] ),
                                tuple( x[-1] for x in haplotypeList[1] )
                            ]
                            for sampleID, haplotypeList in sampleHaplotypes.items()
                        }
                        
                        uniqueHaplotypes = list(set(
                            haplotypeCode
                            for haplotypeCodeList in sampleCodes.values()
                            for haplotypeCode in haplotypeCodeList
                        ))
                        
                        # Extract haplotype sequences for all unique haplotypes
                        haplotypeSequences = {}
                        for haplotypeCode in uniqueHaplotypes:
                            # Find an example sample with this haplotype
                            "The sample value is richer and contains data our codes lack"
                            haplotypeExemplar = [
                                haplotype
                                for haplotypesList in sampleHaplotypes.values()
                                for haplotype in haplotypesList
                                if tuple( x[-1] for x in haplotype ) == haplotypeCode
                            ][0]
                            
                            # Get the edited sequence
                            haplotypeSequence = edit_reference_to_haplotype_sequence(exon_FastASeq_obj.seq[:], 
                                                    haplotypeExemplar, mrnaFeature.strand)
                            haplotypeSequences[haplotypeCode] = haplotypeSequence
                        
                        # Write haplotype sequences to file
                        sampleLoggingInfo = []
                        for sampleID, haplotypeCodeList in sampleCodes.items():
                            sampleSeq1 = haplotypeSequences[haplotypeCodeList[0]]
                            sampleSeq2 = haplotypeSequences[haplotypeCodeList[1]]
                            
                            fileOut.write(f">{sampleID}_haplotype_1" + \
                                    " has_switch_error={0}\n{1}\n".format(
                                        "yes" if sampleSwitchPoints[sampleID] else "no",
                                        sampleSeq1
                                    ))
                            fileOut.write(f">{sampleID}_haplotype_2" + \
                                    " has_switch_error={0}\n{1}\n".format(
                                        "yes" if sampleSwitchPoints[sampleID] else "no",
                                        sampleSeq2
                                    ))
                            
                            # Add info to logging
                            sampleLoggingInfo.append(("*" if sampleSwitchPoints[sampleID] else "") + "-".join(haplotypeCodeList[0]))
                            sampleLoggingInfo.append(("*" if sampleSwitchPoints[sampleID] else "") + "-".join(haplotypeCodeList[1]))
                        
                        # Write logging information
                        logOut.write("{0}\tyes\t{1}\n".format(parentFeature.ID, "\t".join(sampleLoggingInfo)))
        
        open(os.path.join(args.outputDirectory, "fasta_phasing_was_successful.flag"), "w").close()
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
