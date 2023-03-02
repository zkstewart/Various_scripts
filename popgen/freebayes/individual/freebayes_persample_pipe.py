#! python3
# freebayes_persample_pipe.py
# Script to automate the pipeline needed to predict variants using
# Freebayes using a per-sample variant call followed by a per-sample
# fill-in where variant calling is done again with VCF input.
# The need to automate this arises from a bug present in Freebayes
# which complicates things greatly.

import os, argparse, shutil, subprocess

# Define validation functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.bamDirectory):
        print(f'I am unable to locate the BAM directory ({args.bamDirectory})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFile):
        print(f'I am unable to locate the genome FASTA file ({args.genomeFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print(f'I am unable to locate the metadata text file ({args.metadataFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.freebayesExe):
        print(f'I am unable to locate the Freebayes executable file ({args.freebayesExe})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.varScriptDir):
        print(f'I am unable to locate the Various_scripts directory ({args.varScriptDir})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.moduleFile != None:
        if not os.path.isfile(args.moduleFile):
            print(f'I am unable to locate the nodules text file ({args.moduleFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate numeric inputs
    if args.sampleDP < 0:
        print("sampleDP should be an integer greater than or equal to zero")
        quit()
    if not 0 <= args.popMissing <= 1:
        print("popMissing should be a float in the range of 0->1 inclusive")
        quit()

def validate_programs_locatable(programList):
    '''
    Raises an error if any program is not locatable.
    
    Parameters:
        programList -- a list containing one or more strings corresponding to
                       program executable names which are expected to be locatable
                       via the system PATH
    '''
    for program in programList:
        which = shutil.which(program)
        if which == None:
            raise FileNotFoundError(
                f"{program} not found in system PATH!"
            )

def validate_metadata_against_bams(bamDirectory, bamSuffix, metadataFile):
    '''
    Raises an error if any BAMS specified in the metadata file are absent
    (or if any BAMs found in the directory are not specified in the metadata file)
    
    Returns:
        sampleIDs -- a set which provides the sample IDs for this analysis.
                     Useful for knowing how many subjobs to use with batch submission!
    '''
    # Parse metadata file
    metadataSet = set()
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l == "":
                continue
            else:
                sampleID, pop = l.split("\t")
                metadataSet.add(sampleID)
    
    # Parse BAM directory
    bamSet = set()
    for file in os.listdir(bamDirectory):
        if file.endswith(bamSuffix):
            bamSet.add(file.split(bamSuffix)[0])
    
    # Compare sets and end program if needed
    if metadataSet != bamSet:
        print("Metadata and BAM directories do not match up!")
        print("# Samples in metadata not found in BAM directory:")
        for sampleID in metadataSet:
            if not sampleID in bamSet:
                print(sampleID)
        print("# Samples in BAM directory not found in metadata:")
        for sampleID in bamSet:
            if not sampleID in metadataSet:
                print(sampleID)
        quit()
    
    return metadataSet

# Parameter handling containers
class Container:
    def __init__(self, paramsDict):
        for key, value in paramsDict.items():
            self.__dict__[key] = value

class Parameters:
    def __init__(self, size):
        self.VALID_SIZES = ["low", "medium", "high"]
        
        assert size in self.VALID_SIZES, \
            f"{size} not recognised as being a valid size ({self.VALID_SIZES})"
        
        self.set_params(size)
    
    def set_params(self, size):
        self.freebayes_r1_time = ["04:00:00", "12:00:00", "72:00:00"][self.VALID_SIZES.index(size)]
        self.freebayes_r1_mem = ["10G", "15G", "20G"][self.VALID_SIZES.index(size)]
        self.freebayes_r1_cpu = "1"
        
        self.freebayes_r2_time = ["08:00:00", "24:00:00", "90:00:00"][self.VALID_SIZES.index(size)]
        self.freebayes_r2_mem = ["20G", "25G", "30G"][self.VALID_SIZES.index(size)]
        self.freebayes_r2_cpu = "1"
        
        self.normalise_r1_time = ["04:00:00", "08:00:00", "12:00:00"][self.VALID_SIZES.index(size)]
        self.normalise_r1_mem = "10G"
        self.normalise_r1_cpu = "1"
        
        self.normalise_r2_time = ["08:00:00", "24:00:00", "48:00:00"][self.VALID_SIZES.index(size)]
        self.normalise_r2_mem = ["10G", "15G", "20G"][self.VALID_SIZES.index(size)]
        self.normalise_r2_cpu = "1"
        
        self.merge_time = ["08:00:00", "24:00:00", "90:00:00"][self.VALID_SIZES.index(size)]
        self.merge_mem = ["10G", "25G", "40G"][self.VALID_SIZES.index(size)]
        self.merge_cpu = "1"
        
        self.filter_time = ["04:00:00", "08:00:00", "24:00:00"][self.VALID_SIZES.index(size)]
        self.filter_mem = ["70G", "150G", "250G"][self.VALID_SIZES.index(size)]
        self.filter_cpu = "1"
        
        self.split_time = "04:00:00"
        self.split_mem = "10G"
        self.split_cpu = "1"

# Script running and creation functions
def parse_modules(moduleFile):
    '''
    Very simply parses a text file which SHOULD contain module load
    commands but this function doesn't care if it doesn't. That's a
    you problem.
    '''
    modules = []
    with open(moduleFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l != "":
                modules.append(l)
    return "\n".join(modules)

def qsub(scriptFileName):
    qsubProcess = subprocess.Popen(f"qsub {scriptFileName}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    jobID, stderr = qsubProcess.communicate()
    jobID, stderr = jobID.decode(), stderr.decode()
    if stderr == "":
        return jobID.strip(" \r\n")
    else:
        raise Exception(f"qsub died with stderr == {stderr}")

def make_freebayes_r1_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N fbr1_{prefix}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -J 1-{numJobs}
{waitingLine}

cd {workingDir}
{modules}
# Specify the location of Freebayes
FBDIR={fbDir}

# Specify the location of the genome FASTA
GENOMEDIR={genomeDir}
GENOME={genomeFile}

# Specify the location of the mapped BAM files
MAPDIR={bamDir}

# Specify the suffix that identifies mapped BAM files
SUFFIX={suffix}

####

# > STEP 1: Get our file list
declare -a BAMFILES
i=0
for f in ${{MAPDIR}}/*${{SUFFIX}}; do
    BAMFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# > STEP 2: Get our array index
declare -i index
index=${{PBS_ARRAY_INDEX}}-1

# > STEP 3: Get our file for analysis
INPUTFILE=${{BAMFILES[${{index}}]}}

# > STEP 4: Get our output file prefix
PREFIX=$(basename ${{INPUTFILE}} ${{SUFFIX}})

# > STEP 5: Run freebayes
if [[ ! -f ${{PREFIX}}.vcf ]]; then
    ${{FBDIR}}/{fbExe} -f ${{GENOMEDIR}}/${{GENOME}} ${{INPUTFILE}} > ${{PREFIX}}.vcf;
fi
""".format(
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    numJobs=argsContainer.numJobs,
    waitingLine="#PBS -W depend=afterok:{0}".format(argsContainer.prevJobs) if argsContainer.prevJobs != None
        else "",
    
    walltime=argsContainer.walltime,
    mem=argsContainer.mem,
    cpus=argsContainer.cpus,
    
    modules=f"\n####\n{argsContainer.modules}\n####\n" if argsContainer.modules != None
        else "",
    
    fbDir=os.path.dirname(argsContainer.freebayes),
    fbExe=os.path.basename(argsContainer.freebayes),
    
    genomeDir=os.path.dirname(argsContainer.genome),
    genomeFile=os.path.basename(argsContainer.genome),
    
    bamDir=argsContainer.bamDir,
    suffix=argsContainer.suffix
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_normalisation_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N norm_{prefix}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -J 1-{numJobs}
#PBS -W depend=afterok:{prevJobs}

cd {workingDir}

# Specify the location of the genome FASTA
GENOMEDIR={genomeDir}
GENOME={genomeFile}

# Specify the location of the mapped BAM files
MAPDIR={bamDir}

# Specify the suffix that identifies mapped BAM files
SUFFIX={suffix}

#################################

# > STEP 1: Get our file list
declare -a BAMFILES
i=0
for f in ${{MAPDIR}}/*${{SUFFIX}}; do
    BAMFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# > STEP 2: Get our array index
declare -i index
index=${{PBS_ARRAY_INDEX}}-1

# > STEP 3: Get our file prefix
INPUTFILE=${{BAMFILES[${{index}}]}}
PREFIX=$(basename ${{INPUTFILE}} ${{SUFFIX}})

#################################

# > STEP 1: Drop sites which have no coverage
if [[ ! -f ${{PREFIX}}.dp.vcf  ]]; then
    vcffilter -f "DP > 0" ${{PREFIX}}.vcf > ${{PREFIX}}.dp.vcf
fi
if [[ ! -f ${{PREFIX}}.dp.vcf.gz  ]]; then
    bgzip -c ${{PREFIX}}.dp.vcf > ${{PREFIX}}.dp.vcf.gz
fi
if [[ ! -f ${{PREFIX}}.dp.vcf.gz.csi  ]]; then
    bcftools index ${{PREFIX}}.dp.vcf.gz
fi

# > STEP 2: Split multiallelic records to biallelic
if [[ ! -f ${{PREFIX}}.split.vcf.gz  ]]; then
    bcftools norm -m- -Oz -o ${{PREFIX}}.split.vcf.gz -N ${{PREFIX}}.dp.vcf.gz
fi

# > STEP 3: Rejoin biallic sites into multiallelic sites
if [[ ! -f ${{PREFIX}}.rejoin.vcf.gz  ]]; then
    bcftools norm -m+ -Oz -o ${{PREFIX}}.rejoin.vcf.gz -N ${{PREFIX}}.split.vcf.gz
fi

# > STEP 4: Left-align and normalise everything
if [[ ! -f ${{PREFIX}}.normalised.vcf  ]]; then
    bcftools norm -f ${{GENOMEDIR}}/${{GENOME}} -Ov -o ${{PREFIX}}.normalised.vcf ${{PREFIX}}.rejoin.vcf.gz
fi

# > STEP 5: vt decompose SNPs
if [[ ! -f ${{PREFIX}}.decomposed.vcf  ]]; then
    vt decompose_blocksub ${{PREFIX}}.normalised.vcf > ${{PREFIX}}.decomposed.vcf
fi

# > STEP 6: Make file ready for next steps
if [[ ! -f ${{PREFIX}}.decomposed.vcf.gz  ]]; then
    bgzip -c ${{PREFIX}}.decomposed.vcf > ${{PREFIX}}.decomposed.vcf.gz
fi
if [[ ! -f ${{PREFIX}}.decomposed.vcf.gz.csi  ]]; then
    bcftools index ${{PREFIX}}.decomposed.vcf.gz
fi

""".format(
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    numJobs=argsContainer.numJobs,
    prevJobs=argsContainer.prevJobs,
    
    walltime=argsContainer.walltime,
    mem=argsContainer.mem,
    cpus=argsContainer.cpus,
    
    genomeDir=os.path.dirname(argsContainer.genome),
    genomeFile=os.path.basename(argsContainer.genome),
    
    bamDir=argsContainer.bamDir,
    suffix=argsContainer.suffix
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_merge_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N merge_{prefix}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -W depend=afterok:{prevJobs}

cd {workingDir}

# Specify the location of the Various_scripts directory
VARSCRIPTDIR={varScriptDir}

# Specify prefix for the output file
PREFIX={prefix}

#################################

# > STEP 1: Get our file list
declare -a VCFFILES
i=0
for f in *.decomposed.vcf.gz; do
    VCFFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# > STEP 2: Get our input files argument
SEPARATOR=" "
VCFTOOLS_ARG="$( printf "${{SEPARATOR}}%s" "${{VCFFILES[@]}}" )"

# > STEP 3: Merge individual VCFs
if [[ ! -f ${{PREFIX}}.merged.vcf.gz  ]]; then
    bcftools merge -Oz -o ${{PREFIX}}.merged.vcf.gz ${{VCFTOOLS_ARG}}
fi

# > STEP 4: Index VCF
if [[ ! -f ${{PREFIX}}.merged.vcf.gz.csi  ]]; then
    tabix ${{PREFIX}}.merged.vcf.gz;
    tabix -C ${{PREFIX}}.merged.vcf.gz;
fi
""".format(
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    prevJobs=argsContainer.prevJobs,
    
    walltime=argsContainer.walltime,
    mem=argsContainer.mem,
    cpus=argsContainer.cpus,
    
    varScriptDir=argsContainer.varScriptDir
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_vcf_split_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N vcfsplit_{prefix}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -W depend=afterok:{prevJobs}

cd {workingDir}

# Specify the location of the VCF
VCFDIR={vcfDir}

# Specify prefix for input and output files
PREFIX={prefix}

#################################

# > STEP 1: Derive our VCF file name
VCFFILE=${{VCFDIR}}/${{PREFIX}}.merged.vcf.gz

# > STEP 2: Stream VCF contig IDs through bcftools to split
bcftools index -s ${{VCFFILE}} | cut -f 1 | while read C; do
    bcftools view -Oz -o ${{C}}.vcf.gz ${{VCFFILE}} "${{C}}";
    tabix ${{C}}.vcf.gz;
    tabix -C ${{C}}.vcf.gz;
done
""".format(
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    prevJobs=argsContainer.prevJobs,
    
    walltime=argsContainer.walltime,
    mem=argsContainer.mem,
    cpus=argsContainer.cpus,
    
    vcfDir=argsContainer.vcfDir
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_freebayes_r2_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N fbr2_{prefix}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -J 1-{numJobs}
{waitingLine}

cd {workingDir}
{modules}
# Specify the location of Freebayes
FBDIR={fbDir}

# Specify the location of the genome FASTA
GENOMEDIR={genomeDir}
GENOME={genomeFile}

# Specify the location of the prepped, split VCF files
VCFDIR={vcfDir}
VCFSUFFIX=.vcf.gz

# Specify the location of the mapped BAM files
MAPDIR={bamDir}

# Specify the suffix that identifies mapped BAM files
SUFFIX={suffix}

####

# > STEP 1: Get our file list
declare -a BAMFILES
i=0
for f in ${{MAPDIR}}/*${{SUFFIX}}; do
    BAMFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# > STEP 2: Get our array index
declare -i index
index=${{PBS_ARRAY_INDEX}}-1

# > STEP 3: Get our file for analysis
INPUTFILE=${{BAMFILES[${{index}}]}}

# > STEP 4: Get our output file prefix
PREFIX=$(basename ${{INPUTFILE}} ${{SUFFIX}})

# > STEP 5: Set up directory for containing component files
mkdir -p ${{PREFIX}}
cd ${{PREFIX}}

# > STEP 6: Run freebayes per contig with BAM data being streamed in
for VCF in ${{VCFDIR}}/*${{VCFSUFFIX}}; do
    CONTIG=$(basename ${{VCF}} ${{VCFSUFFIX}});
    if [[ ! -f ${{CONTIG}}.vcf ]]; then
        samtools view -b ${{INPUTFILE}} ${{CONTIG}} | \\
            ${{FBDIR}}/{fbExe} -f ${{GENOMEDIR}}/${{GENOME}} -@ ${{VCF}} --only-use-input-alleles --stdin > ${{CONTIG}}.vcf;
    fi
done

# > STEP 7: Merge per contig VCF outputs together
# >> 7.1: Get a list of all VCFs in this directory
declare -a VCFFILES
i=0
for f in *.vcf; do
    VCFFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# >> 7.2: Get our input files argument
SEPARATOR=" "
VCFTOOLS_ARG="$( printf "${{SEPARATOR}}%s" "${{VCFFILES[@]}}" )"

# >> 7.3: Merge individual VCFs
if [[ ! -f ${{PREFIX}}.vcf  ]]; then
    bcftools merge -Ov -o ${{PREFIX}}.vcf ${{VCFTOOLS_ARG}}
fi

# > STEP 8: Make the merged VCF accessible outside of this working directory
cd ..
if [[ ! -L ${{PREFIX}}.vcf  ]]; then
    ln -s ${{PREFIX}}/${{PREFIX}}.vcf .;
fi
""".format(
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    numJobs=argsContainer.numJobs,
    waitingLine="#PBS -W depend=afterok:{0}".format(argsContainer.prevJobs) if argsContainer.prevJobs != None
        else "",
    
    walltime=argsContainer.walltime,
    mem=argsContainer.mem,
    cpus=argsContainer.cpus,
    
    modules=f"\n####\n{argsContainer.modules}\n####\n" if argsContainer.modules != None
        else "",
    
    fbDir=os.path.dirname(argsContainer.freebayes),
    fbExe=os.path.basename(argsContainer.freebayes),
    
    genomeDir=os.path.dirname(argsContainer.genome),
    genomeFile=os.path.basename(argsContainer.genome),
    
    vcfDir=argsContainer.vcfDir,
    
    bamDir=argsContainer.bamDir,
    suffix=argsContainer.suffix
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_filter_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N filter_{prefix}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -W depend=afterok:{prevJobs}

cd {workingDir}

# Specify the location of the Various_scripts directory
VARSCRIPTDIR={varScriptDir}

# Specify the location of the pop map file
POPSFILE={metadataFile}

# Specify prefix for the output file
PREFIX={prefix}

#################################

# > STEP 1: Get our file list
declare -a VCFFILES
i=0
for f in *.decomposed.vcf.gz; do
    VCFFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# > STEP 2: Get our input files argument
SEPARATOR=" "
VCFTOOLS_ARG="$( printf "${{SEPARATOR}}%s" "${{VCFFILES[@]}}" )"

# > STEP 3: Merge individual VCFs
if [[ ! -f ${{PREFIX}}.merged.vcf  ]]; then
    bcftools merge -Ov -o ${{PREFIX}}.merged.vcf ${{VCFTOOLS_ARG}}
fi

# > STEP 4: Filter vcfs according to publication Pete emailed to pro_zac
MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30, whatever that means
MAC=1 ## minor allele count must be >= 1
MINDEPTH=3 ## minimum depth of 3 for a genotype call
MAF=0.05 ## minor allele frequency greater than or equal to 0.05
if [[ ! -f ${{PREFIX}}.filtered.vcf  ]]; then
    vcftools --vcf ${{PREFIX}}.merged.vcf --max-missing ${{MISSING}} --mac ${{MAC}} --minQ ${{MINQ}} --min-meanDP ${{MINDEPTH}} --remove-filtered-all --recode --recode-INFO-all --maf ${{MAF}} --out ${{PREFIX}}.filtered.vcf;
    mv ${{PREFIX}}.filtered.vcf.recode.vcf ${{PREFIX}}.filtered.vcf;
fi

# > STEP 5: Remove indels and keep only biallelic sites
if [[ ! -f ${{PREFIX}}.filtered.noindels.vcf  ]]; then
    bcftools view --max-alleles 2 --exclude-types indels -Ov -o ${{PREFIX}}.filtered.noindels.vcf ${{PREFIX}}.filtered.vcf
fi

# > STEP 6: More custom filtering
POPMISSING=0.5 ## remove a site if each population does not have at least this percentage of called genotypes
SAMPLEDP=3 ## set a site to be ambiguous if it does not have at least this much DP
if [[ ! -f ${{PREFIX}}.final.vcf  ]]; then
    python ${{VARSCRIPTDIR}}/popgen/VCF/filter_vcf.py -v ${{PREFIX}}.filtered.noindels.vcf -p ${{POPSFILE}} -o ${{PREFIX}}.final.vcf --mpp ${{POPMISSING}} --dp ${{SAMPLEDP}}
fi
""".format(
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    prevJobs=argsContainer.prevJobs,
    
    walltime=argsContainer.walltime,
    mem=argsContainer.mem,
    cpus=argsContainer.cpus,
    
    varScriptDir=argsContainer.varScriptDir,
    metadataFile=argsContainer.metadataFile
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def main():
    # User input
    usage = """%(prog)s automates the two-round Freebayes prediction process whereby
    per-sample variants are called, with the merged variants VCF used for a second
    calling procedure to obtain reference allele calls in all samples.
    
    Note: It's expected that vcftools, bcftools, and vt are locatable in your
    system PATH. This script will attempt to verify that before running.
    """
    # Required (file locations)
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-b", dest="bamDirectory",
                   required=True,
                   help="Input directory where BAM files reside")
    p.add_argument("-g", dest="genomeFile",
                   required=True,
                   help="Input genome fasta file")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Input metadata text file")
    # Required (program locations)
    p.add_argument("-f", dest="freebayesExe",
                   required=True,
                   help="Input the full location of the freebayes executable")
    p.add_argument("-v", dest="varScriptDir",
                   required=True,
                   help="Input the location of the Various_scripts directory")
    # Required (program behaviour)
    p.add_argument("-p", dest="prefix",
                   required=True,
                   help="Specify the prefix for output files")
    p.add_argument("-s", dest="bamSuffix",
                   required=True,
                   help="Specify the suffix which identifies BAM files to be used")
    # Optional
    p.add_argument("--modules", dest="moduleFile",
                   required=False,
                   help="""Optionally provide a text file containing module load
                   commands to embed within Freebayes scripts""")
    p.add_argument("--popMissing", dest="popMissing", type=float,
                   required=False,
                   help="""Optionally, during filtering, drop variant rows if each population
                   does not have at least this proportion of individuals genotyped
                   (default == 0.5)""",
                   default=0.5)
    p.add_argument("--sampleDP", dest="sampleDP", type=int,
                   required=False,
                   help="""Optionally, during filtering, delete a sample's variant call
                   if it's depth is lower than this value (default == 1)""",
                   default=1)
    p.add_argument("--size", dest="size",
                   required=False,
                   choices=["low", "medium", "high"],
                   help="""Optionally, specify how much resources you think the jobs will
                   need on a scale from low -> high (default == "medium")""",
                   default="medium")
    args = p.parse_args()
    validate_args(args)
    
    # Validate that PATH programs are locatable
    validate_programs_locatable(["vcftools", "bcftools", "vt"])
    
    # Validate that metadata matches BAM files
    sampleSet = validate_metadata_against_bams(args.bamDirectory, args.bamSuffix, args.metadataFile)
    
    # Get our parameters
    params = Parameters(args.size)
    
    # Load in modules (if applicable)
    if args.moduleFile != None:
        modules = parse_modules(args.moduleFile)
    else:
        modules = None
    
    #### ROUND 1 ####
    
    # Create directory where round 1 will be run
    round1Dir = os.path.join(os.getcwd(), "round1")
    os.makedirs(round1Dir, exist_ok=True)
    
    # Format and qsub freebayes
    fbRound1ScriptName = os.path.join(round1Dir, "run_freebayes_r1.sh")
    make_freebayes_r1_script(Container({
        "prefix": args.prefix,
        "workingDir": round1Dir,
        "prevJobs": None,
        "numJobs": len(sampleSet),
        "walltime": params.freebayes_r1_time,
        "mem": params.freebayes_r1_mem,
        "cpus": params.freebayes_r1_cpu,
        "modules": modules,
        "freebayes": args.freebayesExe,
        "genome": args.genomeFile,
        "bamDir": args.bamDirectory,
        "suffix": args.bamSuffix,
        "outputFileName": fbRound1ScriptName
    }))
    fbRound1JobID = qsub(fbRound1ScriptName)
    
    # Format and qsub normalisation
    normRound1ScriptName = os.path.join(round1Dir, "run_normalise_r1.sh")
    make_normalisation_script(Container({
        "prefix": "r1_" + args.prefix,
        "workingDir": round1Dir,
        "numJobs": len(sampleSet),
        "prevJobs": fbRound1JobID,
        "walltime": params.normalise_r1_time,
        "mem": params.normalise_r1_mem,
        "cpus": params.normalise_r1_cpu,
        "genome": args.genomeFile,
        "bamDir": args.bamDirectory,
        "suffix": args.bamSuffix,
        "outputFileName": normRound1ScriptName,
    }))
    normRound1JobID = qsub(normRound1ScriptName)
    
    # Format and qsub merging
    mergeRound1ScriptName = os.path.join(round1Dir, "run_merge_r1.sh")
    make_merge_script(Container({
        "prefix": args.prefix,
        "workingDir": round1Dir,
        "prevJobs": normRound1JobID,
        "walltime": params.merge_time,
        "mem": params.merge_mem,
        "cpus": params.merge_cpu,
        "varScriptDir": args.varScriptDir,
        "outputFileName": mergeRound1ScriptName
    }))
    mergeRound1JobID = qsub(mergeRound1ScriptName)
    
    #### PREP ROUND ####
    
    # Create directory where prepping will be performed
    prepDir = os.path.join(os.getcwd(), "prep")
    os.makedirs(prepDir, exist_ok=True)
    
    # Split VCF by contig
    vcfSplitScriptName = os.path.join(prepDir, "run_vcf_split.sh")
    make_vcf_split_script(Container({
        "prefix": args.prefix,
        "workingDir": prepDir,
        "prevJobs": mergeRound1JobID,
        "walltime": params.split_time,
        "mem": params.split_mem,
        "cpus": params.split_cpu,
        "vcfDir": round1Dir,
        "outputFileName": vcfSplitScriptName
    }))
    vcfSplitJobID = qsub(vcfSplitScriptName)
    
    #### ROUND 2 ####
    
    # Create directory where round 1 will be run
    round2Dir = os.path.join(os.getcwd(), "round2")
    os.makedirs(round2Dir, exist_ok=True)
    
    # Format and qsub freebayes
    fbRound2ScriptName = os.path.join(round2Dir, "run_freebayes_r2.sh")
    make_freebayes_r2_script(Container({
        "prefix": args.prefix,
        "workingDir": round2Dir,
        "numJobs": len(sampleSet),
        "prevJobs": vcfSplitJobID,
        "walltime": params.freebayes_r2_time,
        "mem": params.freebayes_r2_mem,
        "cpus": params.freebayes_r2_cpu,
        "modules": modules,
        "freebayes": args.freebayesExe,
        "genome": args.genomeFile,
        "bamDir": args.bamDirectory,
        "vcfDir": prepDir,
        "suffix": args.bamSuffix,
        "outputFileName": fbRound2ScriptName
    }))
    fbRound2JobID = qsub(fbRound2ScriptName)
    
    # Format and qsub normalisation
    normRound2ScriptName = os.path.join(round2Dir, "run_normalise_r2.sh")
    make_normalisation_script(Container({
        "prefix": "r2_" + args.prefix,
        "workingDir": round2Dir,
        "numJobs": len(sampleSet),
        "prevJobs": fbRound2JobID,
        "walltime": params.normalise_r2_time,
        "mem": params.normalise_r2_mem,
        "cpus": params.normalise_r2_cpu,
        "genome": args.genomeFile,
        "bamDir": args.bamDirectory,
        "suffix": args.bamSuffix,
        "outputFileName": normRound2ScriptName,
    }))
    normRound2JobID = qsub(normRound2ScriptName)
    
    # Format and qsub merging
    mergeRound2ScriptName = os.path.join(round2Dir, "run_merge_r2.sh")
    make_merge_script(Container({
        "prefix": args.prefix,
        "workingDir": round2Dir,
        "prevJobs": normRound2JobID,
        "walltime": params.merge_time,
        "mem": params.merge_mem,
        "cpus": params.merge_cpu,
        "varScriptDir": args.varScriptDir,
        "outputFileName": mergeRound2ScriptName
    }))
    mergeRound2JobID = qsub(mergeRound2ScriptName)
    
    # Format and qsub filtering
    filterScriptName = os.path.join(round2Dir, "run_filter.sh")
    make_filter_script(Container({
        "prefix": args.prefix,
        "workingDir": round2Dir,
        "prevJobs": mergeRound2JobID,
        "walltime": params.filter_time,
        "mem": params.filter_mem,
        "cpus": params.filter_cpu,
        "varScriptDir": args.varScriptDir,
        "metadataFile": args.metadataFile,
        "outputFileName": filterScriptName
    }))
    filterJobID = qsub(filterScriptName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
