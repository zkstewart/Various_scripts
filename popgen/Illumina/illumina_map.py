#! python3
# illumina_map.py
# Script to handle Illumina read mapping via BWA-MEM.
# It follows the process of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/

import os, argparse, re

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"I am unable to locate the genome FASTA file ({args.fastaFile})")
    if not os.path.isfile(args.metadataFile):
        raise FileNotFoundError(f"I am unable to locate the metadata file ({args.metadataFile})")
    if not os.path.isfile(args.bwa):
        raise FileNotFoundError(f"I am unable to locate the bwa executable file ({args.bwa})")
    if not os.path.isdir(args.fastqDirectory):
        raise FileNotFoundError(f"I am unable to locate the FASTQ files directory ({args.fastqDirectory})")
    # Validate numeric inputs
    if args.cpus < 1:
        raise ValueError(f"--cpus value '{args.cpus}' should be a positive integer (greater than zero)")

def parse_metadata_csv(metadataFile, prefixCol="prefix", idCol="id", smCol="sm"):
    '''
    Parameters:
        metadataFile -- a string indicating the location of a comma or tab-delimited file containing
                        (at least) the three columns indicated by the other parameters.
        prefixCol -- a string indicating the header for the column containing file prefixes
                     which uniquely identify single or paired end reads
        idCol -- a string indicating the header for the column containing unique IDs for samples
        smCol -- a string indicating the header for the column containing genotype IDs for samples
                 that should be grouped together (if not unique)
    Returns:
        prefixes -- a list of the file prefixes
        IDs -- a list of the sample IDs (ID tag)
        SMs -- a list of the sample genotypes (SM tag)
    '''
    prefixes, IDs, SMs = [], [], []
    
    firstLine = True
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if "," in line:
                sl = l.split(",")
            elif "\t" in line:
                sl = l.split("\t")
            else:
                raise ValueError(f"Unable to parse line '{l}' in metadata CSV; please use tabs or commas to separate columns")
            
            # Handle header line
            if firstLine == True:
                if not prefixCol in sl:
                    raise ValueError(f"--prefix '{prefixCol}' not found in metadata CSV")
                elif not idCol in sl:
                    raise ValueError(f"--id '{idCol}' not found in metadata CSV")
                elif not smCol in sl:
                    raise ValueError(f"--sample '{smCol}' not found in metadata CSV")
                
                prefixIndex = sl.index(prefixCol)
                idIndex = sl.index(idCol)
                smIndex = sl.index(smCol)
                firstLine = False
            # Handle content lines
            else:
                prefixValue, idValue, smValue = sl[prefixIndex], sl[idIndex], sl[smIndex]
                if prefixValue == "" or idValue == "" or smValue == "":
                    raise ValueError(f"Empty values found in metadata CSV; please fill in all fields")
                
                if prefixValue in prefixes:
                    raise ValueError(f"Prefix tag '{prefixValue}' isn't unique!")
                if idValue in IDs:
                    raise ValueError(f"ID tag '{idValue}' isn't unique!")
                
                prefixes.append(prefixValue)
                IDs.append(idValue)
                SMs.append(smValue)
    
    # Alert user to how the metadata has been interpreted
    if len(SMs) != len(set(SMs)):
        print("# WARNING: Some samples have the same --sample ID; they may be pooled during variant calling")
    
    return prefixes, IDs, SMs

def get_files_from_prefix(fastqDirectory, prefixes):
    '''
    Parameters:
        fastqDirectory -- a string indicating the location where files (as indicated by their prefix)
                          can be found.
        prefixes -- a list containing one or more file prefixes which point to single or paired end
                    reads.
    Returns:
        fastqFiles -- a list of lists containing one value (if single end) or two values (if paired end)
                      for each sample
    '''
    fileSortRegex = re.compile(r"\d+")
    numFoundWithPrefix = []
    fastqFiles = []
    
    # Get all files in the directory
    files = os.listdir(fastqDirectory)
    for filePrefix in prefixes:
        filesWithPrefix = [os.path.join(fastqDirectory, f) for f in files if f.startswith(filePrefix)]
        
        if len(filesWithPrefix) == 0:
            raise FileNotFoundError(f"Unable to find a file in '{fastqDirectory}' that starts with prefix '{filePrefix}'")
        elif len(filesWithPrefix) > 2:
            raise ValueError(f"Found more than two matches in '{fastqDirectory}' for files that start with prefix '{filePrefix}'")
        
        # Make sure files are ordered appropriately
        filesWithPrefix.sort(key = lambda x: list(map(int, fileSortRegex.findall(x))))
        
        fastqFiles.append(filesWithPrefix)
        numFoundWithPrefix.append(len(filesWithPrefix))
    
    # Check if all files are of the same type (single or paired end)
    if len(set(numFoundWithPrefix)) != 1:
        raise ValueError(f"'{fastqDirectory}' seems to contain a mix of single and paired end reads; " + 
                         "make sure your prefixes uniquely identify sample file(s) or ensure that " +
                         "all samples have their corresponding paired end files.")
    
    # Alert user to the paired end nature of the files
    if len(fastqFiles[0]) == 2:
        print("# Reads are PAIRED end")
    else:
        print("# Reads are SINGLE end")
    
    # Alert user to how many samples were found
    print(f"# Found {len(fastqFiles)} samples")
    
    return fastqFiles

def process_readgroups(speciesIds, genotypes, platform, library, unit):
    readgroups=[]
    for i in range(len(genotypes)):
        readgroups.append("@RG\\tID:{0}\\tSM:{1}\\tPL:{2}\\tLB:{3}\\tPU:{4}".format(speciesIds[i], genotypes[i], platform, library, unit))
    return readgroups

def create_cmd_file(fastqFiles, speciesIds, readgroups, genomeFile, bwa, outputFileName="cmd_illumina_map.txt", cpus=1):    
    # Validations
    if os.path.isfile(outputFileName):
        raise FileExistsError(f"'{outputFileName}' already exists; delete, move, or rename then try again")
    assert len(fastqFiles) == len(readgroups)
    
    # Write to file
    with open(outputFileName, "w") as fileOut:
        for i in range(len(fastqFiles)):
            sid = speciesIds[i]
            fq = " ".join(fastqFiles[i]) # join in case it's paired end
            rg = readgroups[i]
            
            fileOut.write("{bwa} mem -t {cpus} -R '{rg}' {genomeFile} {fq} > {sampleID}.sam\n".format(
                bwa=bwa,
                cpus=cpus,
                rg=rg,
                genomeFile=genomeFile,
                fq=fq,
                sampleID=sid
            ))

def create_shell_script(cmdFile, numJobs, outputFileName="run_illumina_map.sh", cpus=1):
    # Specify hard-coded script features
    jobname = "illumina_map"
    walltime = "12:00:00"
    mem = "25G"
    
    # Setup the script's contents
    formatStr = """#!/bin/bash -l
#PBS -N {jobname}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -J 1-{numJobs}

cd $PBS_O_WORKDIR

eval $(cat {cmdFile} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1)
"""
    # Write to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(formatStr.format(
            jobname=jobname,
            walltime=walltime,
            mem=mem,
            cpus=cpus,
            numJobs=numJobs,
            cmdFile=cmdFile
        ))

def main():
    # User input
    usage = """%(prog)s enables easier running of BWA-MEM with readgroup
    tags for downstream popgen analysis. It will output a shell script amenable
    to batched job submission via PBS.
    
    It will automatically infer whether your samples are single or paired end
    based on the prefixes given in the metadata file.
    
    The --prefix column must contain a unique value that identifies each sample's read files
    (forward and reverse). This is used solely to locate the files in the -d directory.
    
    The --id column must contain unique values for each sample. This will dictate the output
    SAM files and their embedded read group identifiers, which will dictate what the samples
    are called in your VCF files.
    
    During variant calling, if two samples have the same --sample value, they may be pooled together.
    You probably want this to be unique for each sample and be the same as your --id value unless
    you have a good reason to pool samples.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="fastqDirectory",
                   required=True,
                   help="Input directory where fastq files reside")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome fasta file")
    p.add_argument("-csv", dest="metadataFile",
                   required=True,
                   help="Input Metadata CSV or TSV file")
    p.add_argument("-b", dest="bwa",
                   required=True,
                   help="Input the full path to the bwa executable")
    # Optional
    p.add_argument("--prefix", dest="prefixCol",
                   required=False,
                   help="Column name where file prefix is located (default == 'prefix')",
                   default="prefix")
    p.add_argument("--id", dest="idCol",
                   required=False,
                   help="Column name where species ID is located (default == 'id')",
                   default="id")
    p.add_argument("--sample", dest="sampleCol",
                   required=False,
                   help="Column name where sample ID is located (default == 'sm')",
                   default="sm")
    p.add_argument("--platform", dest="platform",
                   required=False,
                   choices=["illumina", "pacbio", "nanopore"], #incomplete list
                   help="String to use for readgroup platform e.g., 'illumina' by default",
                   default="illumina")
    p.add_argument("--library", dest="library",
                   required=False,
                   help="String to use for library e.g., 'lib1' by default",
                   default="lib1")
    p.add_argument("--unit", dest="unit",
                   required=False,
                   help="String to use for unit e.g., 'unit1' by default",
                   default="unit1")
    p.add_argument("--cpus", dest="cpus", type=int,
                   required=False,
                   help="Optionally specify the number of CPUs to run for each file (default == 1)",
                   default=1)
    args = p.parse_args()
    validate_args(args)
    
    # Parse CSV file columns
    prefixes, ids, samples = parse_metadata_csv(args.metadataFile, args.prefixCol, args.idCol, args.sampleCol)
    
    # Validate that files exist and get their locations
    fastqFiles = get_files_from_prefix(args.fastqDirectory, prefixes)
    
    # Create readgroups
    readgroups = process_readgroups(ids, samples, args.platform, args.library, args.unit)
    
    # Create cmd file
    cmdFileName = "cmd_illumina_map.txt"
    create_cmd_file(fastqFiles, ids, readgroups, args.fastaFile, args.bwa, cmdFileName, args.cpus)
    
    # Create shell script
    create_shell_script(cmdFileName, len(readgroups), cpus=args.cpus)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
