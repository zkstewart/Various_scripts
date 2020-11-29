#! python3
# dge_pipe.py
# Pipeline script for performing differential gene expression analysis

import os, argparse, subprocess, shutil, math
from pathlib import Path

# Define classes
class DGE_Meta:
    def __init__(self, metadataFile):
        self.names = []
        self.reads = []
        self.comparisons = [] # Will be overwritten with meta_open()
        self.meta_open(metadataFile)
    
    def meta_open(self, metadataFile):
        with open(metadataFile, "r") as fileIn:
            header = None
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                # Parse header line
                if header == None:
                    header = sl
                    header_codes = ["" for i in range(len(header))] # For body line parsing
                    found = [False, False, 0] # For header validation
                    for i in range(len(header)):
                        colname = header[i]
                        if colname == "sample":
                            if found[0] == False:
                                header_codes[i] = "s"
                                found[0] = True
                            else:
                                print("Multiple sample fields found in metadata file.")
                                print("Make sure this file conforms to the format and try again.")
                                quit()
                        elif colname.startswith("techrep"):
                            header_codes[i] = "t"
                            found[1] = True
                        elif colname.startswith("comparison"):
                            header_codes[i] = "c"
                            found[2] += 1
                        else:
                            print("Metadata file has unknown column titled {0}".format(colname))
                            print("Make sure this file conforms to the format and try again.")
                            quit()
                    if not (found[0] == True and found[1] == True and found[2] > 0):
                        print("Metadata file is missing necessary fields i.e., sample name, \
                            technical replicate, and/or comparisons.")
                        print("Make sure this file conforms to the format and try again.")
                        quit()
                    # Setup comparisons storage
                    self.comparisons = [[] for i in range(found[2])]
                # Parse body lines
                else:
                    comparisonCount = 0
                    for i in range(len(sl)):
                        if header_codes[i] == "s":
                            self.names.append(sl[i])
                            self.reads.append([])
                        elif header_codes[i] == "t":
                            self.reads[-1].append(sl[i].replace("\\", "/")) # Potential problem with paths
                        elif header_codes[i] == "c":
                            self.comparisons[comparisonCount].append(sl[i])
                            comparisonCount += 1
    
    def format_reads_names(self, rnaseqDirectory):
        for i in range(len(self.reads)):
            self.reads[i][0] = Path(rnaseqDirectory, self.reads[i][0]).as_posix()
            self.reads[i][1] = Path(rnaseqDirectory, self.reads[i][1]).as_posix()

class Locations:
    def __init__(self, genomeLocation, rnaseqLocation, trimWorkDir, starWorkDir, countWorkDir, rnaseqFiles):
        self.genomeLocation = genomeLocation
        self.rnaseqLocation = rnaseqLocation
        self.trimWorkDir = trimWorkDir
        self.starWorkDir = starWorkDir
        self.countWorkDir = countWorkDir
        self.rnaseqFiles = rnaseqFiles

# Define functions
def validate_args(args):
    def fail_validation():
        print("Make sure you\'ve typed the file name or location correctly and try again.")
        quit()
    
    # Validate that all arguments are provided
    for key, value in vars(args).items():
        if value == None:
            print(key + " argument was not specified. Fix this and try again.")
            quit()
    
    # Validate input file locations
    if not os.path.isfile(args.genomeFile):
        print("I am unable to locate the genome FASTA file ({0})".format(args.genomeFile))
        fail_validation()
    if not os.path.isfile(args.annotationFile):
        print("I am unable to locate the GFF3 annotation file ({0})".format(args.annotationFile))
        fail_validation()
    if not os.path.isdir(args.rnaseqDirectory):
        print("I am unable to locate the directory containing RNAseq reads ({0})".format(
            args.rnaseqDirectory))
        fail_validation()
    if not os.path.isdir(args.trimmomaticDir):
        print("I am unable to locate the directory containing the Trimmomatic .jar file ({0})".format(
            args.trimmomaticDir))
        fail_validation()
    if not os.path.isfile(Path(args.trimmomaticDir, args.trimmomaticJar).as_posix()):
        print("I am unable to locate the Trimmomatic .jar file in the specified directory ({0})".format(
            Path(args.trimmomaticDir, args.trimmomaticJar).as_posix()))
        fail_validation()
    if not os.path.isdir(args.starDir):
        print("I am unable to locate the directory containing the STAR executable file ({0})".format(
            args.starDir))
        fail_validation()
    if not os.path.isfile(Path(args.starDir, args.starExe).as_posix()):
        print("I am unable to locate the STAR executable file in the specified directory ({0})".format(
            Path(args.starDir, args.starExe).as_posix()))
        fail_validation()
    if not os.path.isdir(args.python2Dir):
        print("I am unable to locate the directory containing the Python 2 executable file ({0})".format(
            args.python2Dir))
        fail_validation()
    if os.path.isfile(Path(args.python2Dir, "python").as_posix()):
        args.python2Exe = "python"
    elif os.path.isfile(Path(args.python2Dir, "python2.7").as_posix()):
        args.python2Exe = "python2.7"
    elif os.path.isfile(Path(args.python2Dir, "python2").as_posix()):
        args.python2Exe = "python2"
    else:
        print("I am unable to locate the Python 2 executable file in the specified directory ({0})".format(
            Path(args.starDir).as_posix()))
        print("I expect there to be a file called \"python\", \"python2.7\", or \"python2\"")
        fail_validation()

def setup_working_directory(baseDir, species, genomeFile, metadata):
    # Create data location dirs
    genomeLocation = Path(baseDir, "genome")
    genomeLocation.mkdir(exist_ok=True)
    
    rnaseqLocation = Path(baseDir, "rnaseq_reads")
    rnaseqLocation.mkdir(exist_ok=True)

    # Symbolic link files to data locations
    newGenomeName = Path(genomeLocation, "{0}.fasta".format(species))
    if not os.path.islink(newGenomeName):
        os.symlink(genomeFile, newGenomeName)
    
    rnaseqFiles = []
    readFilesSuffix = os.path.basename(metadata.reads[0][0]).split(".",  maxsplit=1)[1] # Program assumes all files are same
    for i in range(len(metadata.reads)):
        # Derive new file names with standard format
        newR1Name = Path(rnaseqLocation, "{0}_R1.{1}".format(
                metadata.names[i], readFilesSuffix)).as_posix()
        newR2Name = Path(rnaseqLocation, "{0}_R2.{1}".format(
            metadata.names[i], readFilesSuffix)).as_posix()
        rnaseqFiles.append(newR1Name)
        rnaseqFiles.append(newR2Name)

        # Concatenate when multiple technical replicates exist
        if len(metadata.reads[i]) > 1:
            # Concatenate R1 file
            if not os.path.exists(newR1Name):
                with open(newR1Name, "wb") as r1FileOut:
                    for readFile in metadata.reads[i]:
                        with open(readFile.format(1), "rb") as fileIn:
                            shutil.copyfileobj(fileIn, r1FileOut)
            # Concatenate R2 file
            if not os.path.exists(newR2Name):
                with open(newR2Name, "wb") as r2FileOut:
                    for readFile in metadata.reads[i]:
                        with open(readFile.format(2), "rb") as fileIn:
                            shutil.copyfileobj(fileIn, r2FileOut)
        # Symbolic link when no technical replicates exist
        else:
            readFile = metadata.reads[i][0]
            if not os.path.islink(newR1Name):
                os.symlink(readFile.format(1), newR1Name) # The metadata SHOULD have {} placeholders already present
            if not os.path.islink(newR2Name):
                os.symlink(readFile.format(2), newR2Name)

    # Create working dirs
    trimWorkDir = Path(baseDir, "trimmomatic")
    starWorkDir = Path(baseDir, "star_map")
    countWorkDir = Path(baseDir, "htseq_count")

    trimWorkDir.mkdir(exist_ok=True)
    starWorkDir.mkdir(exist_ok=True)
    countWorkDir.mkdir(exist_ok=True)

    locations = Locations(genomeLocation, rnaseqLocation, trimWorkDir, starWorkDir, countWorkDir, rnaseqFiles)

    return locations

def testing_setup_working_directory(baseDir, species, genomeFile, metadata):
    # Create data location dirs
    genomeLocation = Path(baseDir, "genome")
    genomeLocation.mkdir(exist_ok=True)
    
    rnaseqLocation = Path(baseDir, "rnaseq_reads")
    rnaseqLocation.mkdir(exist_ok=True)

    # Symbolic link files to data locations
    newGenomeName = Path(genomeLocation, "{0}.fasta".format(species))

    
    rnaseqFiles = []
    readFilesSuffix = os.path.basename(metadata.reads[0][0]).split(".",  maxsplit=1)[1] # Program assumes all files are same
    for i in range(len(metadata.reads)):
        # Derive new file names with standard format
        newR1Name = Path(rnaseqLocation, "{0}_R1.{1}".format(
                metadata.names[i], readFilesSuffix)).as_posix()
        newR2Name = Path(rnaseqLocation, "{0}_R2.{1}".format(
            metadata.names[i], readFilesSuffix)).as_posix()
        rnaseqFiles.append(newR1Name)
        rnaseqFiles.append(newR2Name)

        # Concatenate when multiple technical replicates exist
        if len(metadata.reads[i]) > 1:
            pass
        # Symbolic link when no technical replicates exist
        else:
            pass

    # Create working dirs
    trimWorkDir = Path(baseDir, "trimmomatic")
    starWorkDir = Path(baseDir, "star_map")
    countWorkDir = Path(baseDir, "htseq_count")

    trimWorkDir.mkdir(exist_ok=True)
    starWorkDir.mkdir(exist_ok=True)
    countWorkDir.mkdir(exist_ok=True)

    locations = Locations(genomeLocation, rnaseqLocation, trimWorkDir, starWorkDir, countWorkDir, rnaseqFiles)

    return locations

def qsub(scriptName):
    p = subprocess.Popen(["qsub", scriptName])
    stdout, stderr = p.communicate()
    return stdout

def generate_trim_script(scriptName, locations, species, trimDir, trimJar, threads, walltime=24, mem="25G"):
    # Format reads for Trimmomatic purposes
    trimReads = []
    suffix = locations.rnaseqFiles[0].split("_R1.")[1]
    for name in metadata.names:
        prefix = Path(locations.rnaseqLocation, name).as_posix()
        trimReads.append(prefix)

    trimScript = r"""#!/bin/bash -l
#PBS -N trim_{species}
#PBS -l walltime={walltime}:00:00
#PBS -l mem={mem}
#PBS -l ncpus={threads}
#PBS -J 1-{numJobs}

cd {workDir}

### MANUAL SETUP BELOW
## SETUP: Load modules
module load java/1.8.0_92

## SETUP: Specify trimmomatic location
TRIMDIR={trimDir}
TRIMJAR={trimJar}

## SETUP: Specify file prefixes
SPECIES={species}

## SETUP: Specify RNAseq files
RNAFILES="{rnaseqFiles}"
ARRAY=($RNAFILES)
# Note: FILEPREFIXES assumes a file name like ${{PREFIX}}_1.fastq with paired ${{PREFIX}}_2.fastq
### MANUAL SETUP END

### AUTOMATIC SETUP START
## SETUP: Paired-end trimmomatic settings
COMMAND="ILLUMINACLIP:{trimDir}/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"
### AUTOMATIC SETUP END

### RUN PROGRAM
ARRAY_INDEX=$((${{PBS_ARRAY_INDEX}}-1))
FILE=${{ARRAY[${{ARRAY_INDEX}}]}}
BASENAME=$(basename ${{FILE}})

## STEP 1: Run Trimmomatic
java -jar $TRIMDIR/$TRIMJAR PE -threads {threads} -trimlog $SPECIES.logfile ${{FILE}}_R1.{suffix} ${{FILE}}_R2.{suffix} -baseout {workDir}/${{BASENAME}}.trimmed.fq.gz ${{COMMAND}}

## STEP 2: Unzip files
gunzip {workDir}/${{BASENAME}}.trimmed_1P.fq.gz {workDir}/${{BASENAME}}.trimmed_2P.fq
    """.format(species=species, walltime=walltime, mem=mem, threads=threads, workDir=locations.trimWorkDir,
            trimDir=trimDir, trimJar=trimJar, rnaseqFiles=" ".join(trimReads), suffix=suffix,
            numJobs=len(trimReads))
    
    with open(scriptName, "w") as fileOut:
        fileOut.write(trimScript)

def generate_star_script(scriptName, starDir, starExe, locations, metadata, species, threads, previousJob, walltime=12, mem="50G"):
    # Format reads for STAR purposes
    starReads = []
    for name in metadata.names:
        # Derive Trimmomatic output name
        prefix = "{0}/{1}.trimmed_".format(locations.trimWorkDir, name)
        starReads.append(prefix)
    
    starScript = r"""#!/bin/bash -l
#PBS -N star_{species}
#PBS -l walltime={walltime}:00:00
#PBS -l mem={mem}
#PBS -l ncpus={threads}
#PBS -W depend=afterok:{previousJob}
#PBS -J 1-{numJobs}

cd {workDir}

## MANUAL SETUP BELOW
# SETUP: Specify STAR location
STARDIR={starDir}

# SETUP: Specify prefix
SPECIES={species}

# SETUP: Specify computational resources
CPUS={threads}

# SETUP: Genome location
GENDIR={genomeLocation}
GENFILE={species}.fasta

## SETUP: Specify RNAseq files
RNAFILES="{rnaseqFiles}"
ARRAY=($RNAFILES)
## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Copy genome here. Need to do this since STAR can only tolerate 1 index per directory...
cp $GENDIR/$GENFILE .

# STEP 2: Generate index
$STARDIR/source/{starExe} --runThreadN $CPUS --runMode genomeGenerate --genomeDir {workDir} --genomeFastaFiles $GENFILE

# STEP 3: Run 2-pass procedure for each sample
ARRAY_INDEX=$((${{PBS_ARRAY_INDEX}}-1))
FILE=${{ARRAY[${{ARRAY_INDEX}}]}}
BASENAME=$(basename ${{FILE}} .trimmed_)
mkdir $BASENAME
cd $BASENAME
$STARDIR/source/{starExe} --runThreadN $CPUS --genomeDir {workDir} --readFilesIn ${{FILE}}1P.fq ${{FILE}}2P.fq --twopassMode Basic
cd ..
""".format(species=species, walltime=walltime, mem=mem, threads=threads, workDir=locations.starWorkDir, 
            starDir=starDir, genomeLocation=locations.genomeLocation, rnaseqFiles=" ".join(starReads),
            previousJob=previousJob, starExe=starExe, numJobs=len(starReads))

    with open(scriptName, "w") as fileOut:
        fileOut.write(starScript)

def generate_htseq_script(scriptName, locations, species, py2Dir, annotationFile, samFile, previousJob, walltime=10, mem="50G"):
    htseqScript = r"""#!/bin/bash -l
#PBS -N htseq_{species}
#PBS -l walltime={walltime}:00:00
#PBS -l mem={mem}
#PBS -l ncpus=1
#PBS -W depend=afterok:{previousJob}

cd {workDir}

## MANUAL SETUP BELOW
# SETUP: Specify mapped reads SAM file
SAMFILE={samFile}

# SETUP: Specify Python2 location
PY2DIR={py2Dir}

# SETUP: GFF annotation file
GFF={annotationFile}

# SETUP: Output file name
OUT={species}.htseq.counts
## MANUAL SETUP END

## RUN PROGRAM
$PY2DIR/python2.7 -m HTSeq.scripts.count -r name -a 0 -t gene -i ID $SAMFILE $GFF > $OUT
    """.format(species=species, walltime=walltime, mem=mem, workDir=locations.countWorkDir, py2Dir=py2Dir,
            samFile=samFile, annotationFile=annotationFile, previousJob=previousJob)

    with open(scriptName, "w") as fileOut:
        fileOut.write(htseqScript)

# Argument input
usage = """%(prog)s does things.
It will establish a working environment within the directory where this script is located. As such, you should
perform a single DGE analysis within the present directory.

Also note that skipping steps should only be done if they were done previously by this script, or if you have set
up the working directory to emulate how it would be laid out were this script run previously.
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", dest="genomeFile",
                help="Specify the location and name of the genome file")
p.add_argument("-a", dest="annotationFile",
                help="Specify the location and name of the gene annotation GFF3 file")
p.add_argument("-m", dest="metadataFile",
                help="Specify the location and name of the sample metadata file")
p.add_argument("-c", dest="cpus", type=int,
                help="Specify the number of cpus/threads to use",
                default = 1)
p.add_argument("-rd", dest="rnaseqDirectory",
                help="Specify the base directory for the RNAseq reads file(s). \
                    If the full path is listed in the metadata file, leave this blank",
                default = "")
p.add_argument("-sp", dest="species",
                help="Provide a short, descriptive name for the species being analysed")
p.add_argument("-td", dest="trimmomaticDir",
                help="Specify the location of the Trimmomatic jar file")
p.add_argument("-te", dest="trimmomaticJar",
                help="Optionally specify the name of the Trimmomatic jar file",
                default = "trimmomatic-0.36.jar")
p.add_argument("-sd", dest="starDir",
                help="Specify the location of the STAR executable file")
p.add_argument("-se", dest="starExe",
                help="Optionally specify the name of the STAR executable file",
                default = "STAR")
p.add_argument("-p2", dest="python2Dir",
                help="Specify the location of the Python 2 directory which has HTSeq installed",
                default = "STAR")
p.add_argument("-skip_trim", dest="skip_trim", action="store_true",
                help="Optionally skip trimming if already performed")
p.add_argument("-skip_map", dest="skip_map", action="store_true",
                help="Optionally skip mapping if already performed")
p.add_argument("-skip_count", dest="skip_count", action="store_true",
                help="Optionally skip counting if already performed")
#args = p.parse_args()
#validate_args(args)

## Hard-coded testing
class Arg:
    def __init__(self):
        self.genomeFile = r"F:\plant_rnaseq\genomes\alm.fasta"
        self.annotationFile = r"F:\plant_rnaseq\annotations\Prudul26A.chromosomes.gff3"
        self.metadataFile = r"F:\plant_rnaseq\sample_metadata\alm_meta.txt"
        self.cpus=1
        self.rnaseqDirectory = r"F:\plant_rnaseq\reads"
        self.species = "alm"
        self.trimmomaticDir = r"D:\Bioinformatics\Protein_analysis\Trimmomatic-0.36"
        self.trimmomaticJar = "trimmomatic-0.36.jar"
        self.starDir = r"D:\Bioinformatics\Protein_analysis\STAR-2.7.6a\bin\Linux_x86_64"
        self.starExe = "STAR"
        self.python2Dir = r"D:\Bioinformatics\Anaconda_2"
        self.skip_trim = False
        self.skip_map = False
        self.skip_count = False

args = Arg()

# Parse metadata
metadata = DGE_Meta(args.metadataFile)
metadata.format_reads_names(args.rnaseqDirectory)

# Setup directory
baseDir = os.getcwd()
#locations = setup_working_directory(baseDir, args.species, args.genomeFile, metadata)
locations = testing_setup_working_directory(baseDir, args.species, args.genomeFile, metadata)

# Trim reads
if not args.skip_trim:
    trimScriptName = Path(locations.trimWorkDir, "run_trim.sh").as_posix()
    generate_trim_script(trimScriptName, locations, args.species, args.trimmomaticDir,
        args.trimmomaticJar, math.ceil(args.cpus / 4)) # Arbitrary division by 4 for batch submission
    trimJob = qsub(trimScriptName)
else:
    trimJob = ""

# Map reads
if not args.skip_map:
    starScriptName = Path(locations.starWorkDir, "run_star.sh").as_posix()
    generate_star_script(starScriptName, args.starDir, args.starExe, locations, metadata,
        args.species, math.ceil(args.cpus / 4), trimJob) # Arbitrary division by 4 for batch submission
    starJob = qsub(starScriptName)
else:
    starJob = ""

# Count mapped reads
if not args.skip_count:
    htseqScriptName = Path(locations.countWorkDir, "run_htseq.sh").as_posix()
    generate_htseq_script(htseqScriptName, locations, args.species, args.python2Dir, args.annotationFile,
        Path(locations.starWorkDir, "Aligned.out.sam").as_posix(), starJob)
    countJob = qsub(htseqScriptName)
else:
    countJob = ""

# DGE analysis

# DGE visualisation

# GOseq analysis

# GOseq visualisation

# Visualisation

