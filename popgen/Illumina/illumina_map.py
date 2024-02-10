#! python3
# illumina_map.py
# Script to handle Illumina read mapping via BWA-MEM.
# It follows the process of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/

import os, argparse, re

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the genome FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataCsv):
        print('I am unable to locate the metadata CSV file (' + args.metadataCsv + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.bwa):
        print('I am unable to locate the bwa executable file (' + args.bwa + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.fastqDirectory):
        print('I am unable to locate the FASTQ files directory (' + args.fastqDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if args.cpus < 1:
        print("cpus should be a positive integer (greater than zero)")
        quit()

def parse_metadata_csv(metadataCsv, prefixCol="prefix", idCol="id", smCol="sm"):
    '''
    Parameters:
        metadataCsv -- a string indicating the file location of a .csv file containing
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
    with open(metadataCsv, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split(",")
            # Handle header line
            if firstLine == True:
                if any([col not in sl for col in [prefixCol, idCol, smCol]]):
                    print("Provided columns are not found in the metadata CSV")
                    print("Make sure --prefix, --id, and --sample are specified correctly")
                    quit()
                prefixIndex = sl.index(prefixCol)
                idIndex = sl.index(idCol)
                smIndex = sl.index(smCol)
                firstLine = False
            # Handle content lines
            else:
                prefix, id, sm = sl[prefixIndex], sl[idIndex], sl[smIndex]
                assert id not in IDs, "ID tag '{0}' isn't unique!".format(id)
                
                prefixes.append(prefix)
                IDs.append(id)
                SMs.append(sm)
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
    
    files = os.listdir(fastqDirectory)
    for filePrefix in prefixes:
        filesWithPrefix = [os.path.join(fastqDirectory, f) for f in files if f.startswith(filePrefix)]
        
        if len(filesWithPrefix) == 0:
            print(f"Error: Unable to find a file in '{fastqDirectory}' that starts with prefix '{filePrefix}'")
            quit()
        elif len(filesWithPrefix) > 2:
            print(f"Error: Found more than two matches in '{fastqDirectory}' for files that start with prefix '{filePrefix}'")
            quit()
        
        # Make sure files are ordered appropriately
        filesWithPrefix.sort(key = lambda x: list(map(int, fileSortRegex.findall(x))))
        
        fastqFiles.append(filesWithPrefix)
        numFoundWithPrefix.append(len(filesWithPrefix))
    
    if len(set(numFoundWithPrefix)) != 1:
        print("Error: There seems to be a mix of single and paired end reads.")
        print("Make sure your prefixes uniquely identify only one type of file and try again.")
        quit()
    
    return fastqFiles

def process_readgroups(speciesIds, genotypes, platform, library, unit):
    readgroups=[]
    for i in range(len(genotypes)):
        readgroups.append("@RG\\tID:{0}\\tSM:{1}\\tPL:{2}\\tLB:{3}\\tPU:{4}".format(speciesIds[i], genotypes[i], platform, library, unit))
    return readgroups

def create_cmd_file(fastqFiles, speciesIds, readgroups, genomeFile, bwa, outputFileName="cmd_illumina_map.txt", cpus=1):    
    # Validations
    if os.path.isfile(outputFileName):
        raise FileExistsError("File name <{0}> already exists".format(outputFileName))
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
    formatStr = """
#!/bin/bash -l
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
    
    This script is suitable for use with single end and paired reads. It will automatically
    infer this based on the prefixes given in the metadata file.
    
    The --prefix column should contain a unique value that identifies each sample's read files
    (forward and reverse).
    
    The --id column should contain unique values for each sample.
    
    The --sample column should contain genotype values which group together
    samples that should be treated as coming from the same material; usually this
    should be unique unless you want these samples to be (essentially) concatenated.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="fastqDirectory",
                   required=True,
                   help="Input directory where fastq files reside")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome fasta file")
    p.add_argument("-csv", dest="metadataCsv",
                   required=True,
                   help="Input Metadata CSV file")
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
                   choices=["illumina", "pacbio"], #incomplete list
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
    prefixes, ids, samples = parse_metadata_csv(args.metadataCsv, args.prefixCol, args.idCol, args.sampleCol)
    
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
