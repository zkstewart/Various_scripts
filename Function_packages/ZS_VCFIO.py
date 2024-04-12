#! python3
# ZS_VCFIO.py
# Contains the VCF class which encapsulated a relatively
# simple dictionary data structure. Use this for parsing
# and filtering a VCF in a performant but high-memory
# use way.

import os, codecs, gzip, subprocess, platform
from contextlib import contextmanager
from ZS_SeqIO import FastASeq
import ZS_Utility

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            return "utf-16"
        except UnicodeDecodeError:
            print(f"VCF class can't tell what codec '{fileName}' is!!")

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

class StandardProgramRunners:
    '''
    A class of static methods to help with running standard programs
    involved in VCF creation or editing.
    '''
    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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

    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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
    
    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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

    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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

    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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

    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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

    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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

    @staticmethod
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
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
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

class SNPStatics:
    '''
    A class of static methods to help with using SNP data as represented
    in VCF files.
    '''
    @staticmethod
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
            featureType -- a string indicating the type of feature we're localising to
                        e.g., a child feature of the mRNA like "CDS" or "exon".
        '''
        
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
    
    @staticmethod
    def edit_reference_to_haplotype_sequence(referenceSeq, haplotypeEditList,
                                         strand, VALIDATE_STRICT=False):
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
                            this is useful during testing and development but shouldn't be
                            necessary if confident that the code functions correctly
        Returns:
            editedSeq -- an edited version of the input referenceSeq with all variations made
        '''
        
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
        editedSeq = referenceSeq if strand == "+" else FastASeq.get_reverse_complement(self=None, staticSeq=referenceSeq)
        return editedSeq.upper() # we'd like all sequences to be upper cased

class SimpleGenotypeIterator:
    '''
    A class to iterate through a VCF file and yield the genotype calls.
    The first value yielded is the sample IDs, with subsequent yields
    providing the following parameters:
        - chrom -- a string indicating the chromosome
        - pos -- an int indicating the position
        - ref -- a string indicating the reference allele
        - alt -- a list of strings indicating the alternate allele(s)
        - posGenotypeDict -- a dictionary with sample IDs as keys and
                             genotype calls as lists of integers (0 for ref, 1 for alt,
                             and so on...)
    '''
    def __init__(self, file_location, impute_missing=False):
        assert os.path.exists(file_location), \
            f"SimpleIterator class can't find file existing at '{file_location}'"
        
        self.fileLocation = file_location
        self.imputeMissing = impute_missing
    
    def parse(self):
        with open_vcf_file(self.fileLocation) as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n").replace('"', '').split("\t") # replace quotations may help with Excel files
                
                # Handle header lines
                if line.startswith("#CHROM"):
                    samples = l[9:] # This gives us the ordered sample IDs
                    yield samples
                    continue
                if line.startswith("#"):
                    continue
                
                # Extract relevant details of the SNP
                chrom = l[0]
                pos = int(l[1])
                ref = l[3]
                alt = l[4].split(",")
                
                # Determine which field position we're extracting to get our GT value
                fieldsDescription = l[8]
                if ":" not in fieldsDescription:
                    gtIndex = -1
                else:
                    gtIndex = fieldsDescription.split(":").index("GT")
                
                # Format a dictionary to store sample genotypes for this position
                posGenotypeDict = {}
                ongoingCount = 0 # This gives us the index for our samples header list 
                for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                    # Grab our genotype
                    if gtIndex != -1:
                        genotype = sampleResult.split(":")[gtIndex]
                    else:
                        genotype = sampleResult
                    
                    # Edit genotype to have a consistently predictable separator
                    "We don't care if the VCF is phased or not for this function"
                    genotype = genotype.replace("/", "|")
                    
                    # Impute empty genotypes (if applicable) or skip otherwise
                    if self.imputeMissing == True:
                        genotype = genotype.replace(".", "0")
                    else:
                        if "." in genotype:
                            ongoingCount += 1
                            continue
                    
                    # Parse and store genotype
                    samplePopulation = samples[ongoingCount]
                    posGenotypeDict[samplePopulation] = list(map(int, genotype.split("|")))
                    
                    ongoingCount += 1
                
                # Check to see if this genotype, after imputation, still has a variant allele
                alleles = set([allele for key, allelePair in posGenotypeDict.items() if key != "ref_alt" for allele in allelePair])
                if alleles == {0}: # skip if it doesn't
                    continue
                
                # Yield result
                yield chrom, pos, ref, alt, posGenotypeDict
    
    def __iter__(self):
        for x in self.parse():
            yield x

class VCF:
    def __init__(self, file_location):
        assert os.path.exists(file_location), \
            f"VCF class can't find file existing at '{file_location}'"
        
        self.fileLocation = file_location
        
        self.comments = {
            "header": [], # contains everything prior to the ##contig= comments
            "contigs": [], # contains the ##contig= comments
            "footer": [] # contains everything after the ##contig= comments
        }
        self.samples = [] # contains the sample IDs as defined by the #CHROM line from splitLine[9:]
        self.variants = {} # contains the variant calls indexed by chrom->pos
        
        self.isVCF = True
        self.parse_vcf() 
    
    def parse_vcf(self):
        '''
        This function will parse a VCF file and set a dictionary which can be used
        for various purposes
        
        Parameters:
            self.fileLocation -- a string indicating the file location of the VCF file
        Sets:
            variants -- a dictionary with structure like:
                        {
                            'chrom1': {
                                pos1: {"REF": ___, "ALT": ___, "GT": ___, "AD: ___, ...},
                                pos2: { ... },
                                ...
                            },
                            'chrom2': ...,
                            ...
                        }
        '''
        beforeContigComments = True
        with open_vcf_file(self.fileLocation) as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n ")
                sl = l.split("\t")
                
                # Store comment lines prior to contig= lines
                if beforeContigComments == True and line.startswith("#") and not "contig=" in line:
                    self.comments["header"].append(l)
                    continue
                # Store comment lines after contig= lines
                elif beforeContigComments == False and line.startswith("#") and not line.startswith("#CHROM") and not "contig=" in line:
                    self.comments["footer"].append(l)
                    continue
                # Store contig= comment lines
                elif line.startswith("#") and "contig=" in line:
                    self.comments["contigs"].append(l)
                    beforeContigComments = False
                    continue
                # Store #CHROM header line
                elif line.startswith("#CHROM"):
                    self.samples = sl[9:]
                    continue
                # Raise error if we didn't hit any comment conditions
                elif line.startswith("#"):
                    raise Exception("Found a comment line I am not equipped to handle --> {0}".format(line))
                
                # Handle content lines
                else:
                    chrom, pos, id, ref, alt, \
                        qual, filter, info, format = sl[0:9]
                    self.variants.setdefault(chrom, {})
                    self.variants[chrom][pos] = {
                        "ID": id,
                        "REF": ref,
                        "ALT": alt,
                        "QUAL": qual,
                        "FILTER": filter,
                        "INFO": info.split(";"),
                        "FORMAT": format.split(":")
                    }
                    
                    # Store sample details in compact format
                    for i in range(len(self.samples)):
                        sampleDetails = sl[9+i].split(":")
                        self.variants[chrom][pos][self.samples[i]] = sampleDetails # sampleIDs[i] == sampleName
                    
                    # If DP is listed in info but not per-sample, impute it now
                    infoDP = VCF.parse_info(info.split(";"), "DP")
                    if infoDP != None and "DP" not in self.variants[chrom][pos]["FORMAT"] and "AD" in self.variants[chrom][pos]["FORMAT"]:
                        self.variants[chrom][pos]["FORMAT"].append("DP")
                        self._impute_dp(chrom, pos, infoDP)
    
    @staticmethod
    def parse_info(info, tag):
        '''
        Helper static function to look through the INFO field for a specific tag's value,
        returning said value.
        
        Parameters:
            info -- a list indexed by self.variants[chrom][pos]["INFO"]
            tag -- the specific {tag}= to find a value for. If it doesn't exist, instead returns None
        '''
        for tagValue in info:
            if tagValue.startswith(f"{tag}="):
                return tagValue.split("=")[1]
        return None
    
    def _impute_dp(self, chrom, pos, infoDP):
        '''
        Simply put, imputes DP for each sample if it's not found for any samples.
        '''
        adIndex = self.variants[chrom][pos]["FORMAT"].index("AD")
        
        sampleADs = [
            self.variants[chrom][pos][id][adIndex]
                for id in self.samples
        ]
        
        adsRatio = [
            sum(map(int, ad.split(","))) / int(infoDP)
                for ad in sampleADs
        ]
        
        dpsImpute = [
            int(adR*int(infoDP))
                for adR in adsRatio
        ]
        
        for i in range(len(self.samples)):
            self.variants[chrom][pos][self.samples[i]].append(str(dpsImpute[i]))
    
    def del_variant(self, contig, pos):
        '''
        Parameters:
            contig -- a string indicating a contigID that is indexed in self.variants
            pos -- a string or int indicating a position that is indexed in
                   self.variants[contig]
        '''
        assert contig in self.variants, \
            f"Cannot delete {contig}-{pos} since {contig} doesn't exist in the VCF"
        assert str(pos) in self.variants[contig], \
            f"Cannot delete {pos} from {contig} since {pos} isn't indexed under that contig"
        
        # Delete the variant
        del self.variants[contig][pos]
        
        # Delete the contig if no more variants exist
        if len(self.variants[contig]) == 0:
            self.del_contig(contig)
    
    def del_sample(self, sampleID):
        '''
        Simply put, eliminates a sample from the dictionary structure.
        
        Parameters:
            sampleID -- a string indicating a sample that is listed in self.samples
        '''
        assert sampleID in self.samples, \
            f"'{sampleID}' doesn't exist in the VCF"
        
        for contigID, contigDict in self.variants.items():
            for pos, posDict in contigDict.items():
                del posDict[sampleID]
        
        del self.samples[self.samples.index(sampleID)]
    
    def del_contig(self, contig):
        '''
        Removes the contig both from self.variants, as well as self.comments["contigs"]
        
        Parameters:
            contig -- a string indicating a contigID that is indexed in self.variants
        '''
        del self.variants[contig]
        
        contigComment = [contigComment for contigComment in self.comments["contigs"] if f"contig=<ID={contig}" in contigComment]
        assert len(contigComment) == 1, \
            f"Somehow {contig} matches against more than one VCF contig comment? Hits = {contigComment}"
        
        commentIndex = self.comments["contigs"].index(contigComment[0])
        del self.comments["contigs"][commentIndex]
    
    def write_vcf(self, outputFileName):
        '''
        Takes the VCF represented by this object and write it to a new VCF file.
        Implicitly makes sure the file is ordered by contig ID sorting.
        '''
        assert not os.path.isfile(outputFileName), \
            f"VCF won't write to '{outputFileName}' since a file already exists here!"
        
        with open(outputFileName, "w") as fileOut:
            # Get our contigs comment line
            contigComments = [cComment for cComment in self.comments["contigs"] if cComment.split(",")[0].split("ID=")[1] in self.variants]
            
            # Write comments
            fileOut.write("\n".join(self.comments["header"]))
            if len(self.comments["header"]) != 0: # spacer
                fileOut.write("\n")
                
            fileOut.write("\n".join(sorted(
                contigComments, key=lambda x: x.split(",")[0].split("ID=")[1]
            )))
            if len(contigComments) != 0: # spacer
                fileOut.write("\n")
            
            fileOut.write("\n".join(self.comments["footer"]))
            if len(self.comments["footer"]) != 0: # spacer
                fileOut.write("\n")
            
            fileOut.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}\n".format(
                "\t".join(self.samples)
            ))
            
            # Write contents
            for contigID in sorted(self.variants.keys()):
                for pos, posDict in self.variants[contigID].items():
                    lineOut = "{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}".format(
                        CHROM=contigID,
                        POS=pos,
                        ID=posDict["ID"],
                        REF=posDict["REF"],
                        ALT=posDict["ALT"],
                        QUAL=posDict["QUAL"],
                        FILTER=posDict["FILTER"]
                    )
                    lineOut += "\t{0}\t{1}\t{2}\n".format(
                        ";".join(posDict["INFO"]),
                        ":".join(posDict["FORMAT"]),
                        "\t".join(
                            [":".join(posDict[sID]) for sID in self.samples]
                        )
                    )
                    fileOut.write(lineOut)
    
    def write_geno(self, outputFileName):
        '''
        Takes the VCF represented by this object and write it to a new .geno file.
        Implicitly makes sure the file is ordered by contig ID sorting.
        '''
        assert not os.path.isfile(outputFileName), \
            f"VCF won't write to '{outputFileName}' since a file already exists here!"
        
        with open(outputFileName, "w") as fileOut:
            # Write comments
            fileOut.write("#CHROM\tPOS\t{0}\n".format(
                "\t".join(self.samples)
            ))
            
            # Write contents
            for contigID in sorted(self.variants.keys()):
                for pos, posDict in self.variants[contigID].items():
                    # Extract relevant info
                    gtOptions = [posDict["REF"], *posDict["ALT"].split(",")] # This enables us to deal with multi-allelic calling
                    gtIndex = posDict["FORMAT"].index("GT")
                    
                    # Parse genotype per sample
                    genotypes = []
                    for sampleID in self.samples:
                        if posDict[sampleID] == ["."]: # handle freebayes non-called
                            genotypes.append("N/N")
                            continue
                        
                        # Multi-allelic compatible genotype finding
                        sampleGT = posDict[sampleID][gtIndex]
                        for i in range(0, len(gtOptions)):
                            sampleGT = sampleGT.replace(str(i), gtOptions[i]).replace(".", "N") # need to switch to N to satisfy the .geno requirements
                        
                        # Store results
                        genotypes.append(sampleGT)
                    
                    # Format output line
                    lineOut = "{CHROM}\t{POS}\t{GENOTYPES}\n".format(
                        CHROM=contigID,
                        POS=pos,
                        GENOTYPES="\t".join(genotypes)
                    )
                    fileOut.write(lineOut)
    
    def __getitem__(self, key):
        return self.variants[key]
    
    def __setitem__(self, key, item):
        self.variants[key] = item
    
    def __len__(self):
        return len(self.variants)
    
    def __iter__(self):
        return iter(self.variants)
    
    def __contains__(self, item):
        return item in self.variants
    
    def has_key(self, key):
        return key in self.variants
    
    def keys(self):
        return self.variants.keys()
    
    def values(self):
        return self.variants.values()
    
    def items(self):
        return self.variants.items()
    
    def __repr__(self):
        return "<VCF object;file='{0}';num_contigs={1};num_variants={2}>".format(
            self.fileLocation,
            len(self.variants),
            sum([len(v) for v in self.variants.values()])
        )

class PhasedVCF:
    def __init__(self, file_location):
        assert os.path.exists(file_location), \
            f"PhasedVCF class can't find file existing at '{file_location}'"
        
        self.fileLocation = file_location
        
        self.samples = [] # contains the sample IDs as defined by the #CHROM line from splitLine[9:]
        self.variants = {} # contains the variant calls indexed by chrom->pos->sampleID
        
        self.isPhasedVCF = True
    
    def parse_whatshap_vcf(self):
        '''
        This function will parse a VCF file produced by WhatsHap.
        
        Parameters:
            self.fileLocation -- a string indicating the file location of the VCF file
        Sets:
            variants -- a dictionary with structure like:
                        ## TBD?
                        {
                            'chrom1': {
                                pos1: {"REF": ___, "ALT": ___, "GT": ___, "AD: ___, ...},
                                pos2: { ... },
                                ...
                            },
                            'chrom2': ...,
                            ...
                        }
        '''
        with open_vcf_file(self.fileLocation) as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n ")
                sl = l.split("\t")
                
                # Store #CHROM header line
                if line.startswith("#CHROM"):
                    self.samples = sl[9:]
                    continue
                # Skip over any other comment lines
                elif line.startswith("#"):
                    continue
                
                # Handle content lines
                else:
                    chrom, pos, id, ref, alt, \
                        qual, filter, info, format = sl[0:9]
                    pos = int(pos)
                    formatList = format.split(":")
                    ref_alt = [ref, *alt.split(",")]
                    
                    # Establish dictionary structure
                    self.variants.setdefault(chrom, {})
                    self.variants[chrom].setdefault(pos, {x: None for x in self.samples})
                    
                    # Store sample details if appropriate
                    for i in range(len(self.samples)):
                        sampleDetailsDict = { _format:_value for _format,_value in zip(formatList, sl[9+i].split(":")) }
                        sampleDetailsDict["phased"] = True if \
                            ("|" in sampleDetailsDict["GT"] or len(set(sampleDetailsDict["GT"].split("/"))) == 1) \
                            else False
                        sampleDetailsDict["ref_alt"] = ref_alt
                        sampleDetailsDict["GT"] = [ ref_alt[int(x)] for x in sampleDetailsDict["GT"].replace("|", "/").split("/") ]
                        
                        # Store the first occurrence of this position
                        if self.variants[chrom][pos][self.samples[i]] == None:
                            self.variants[chrom][pos][self.samples[i]] = sampleDetailsDict
                        # Handle duplicated positions
                        else:
                            if self.variants[chrom][pos][self.samples[i]]["phased"] == False:
                                # Overwrite the unphased call with the phased call
                                "We assume WhatsHap's ability to phase this call means it's the correct one"
                                if sampleDetailsDict["phased"] == True:
                                    self.variants[chrom][pos][self.samples[i]] = sampleDetailsDict
                            else:
                                if sampleDetailsDict["phased"] == True:
                                    raise Exception((f"Found a duplicate phased call at {chrom}:{pos} for sample {self.samples[i]}, " + 
                                                    "both of which were phased; I don't know which call to use!"))
                                
                                # Overwrite the call with the one with the higher combined AD
                                elif sum(map(int, sampleDetailsDict["AD"].split(","))) > \
                                    sum(map(int, self.variants[chrom][pos][self.samples[i]]["AD"].split(","))):
                                        self.variants[chrom][pos][self.samples[i]] = sampleDetailsDict
    
    def __getitem__(self, key):
        return self.variants[key]
    
    def __setitem__(self, key, item):
        self.variants[key] = item
    
    def __len__(self):
        return len(self.variants)
    
    def __iter__(self):
        return iter(self.variants)
    
    def __contains__(self, item):
        return item in self.variants
    
    def has_key(self, key):
        return key in self.variants
    
    def keys(self):
        return self.variants.keys()
    
    def values(self):
        return self.variants.values()
    
    def items(self):
        return self.variants.items()
    
    def __repr__(self):
        return "<PhasedVCF object;file='{0}';num_contigs={1};num_variants={2}>".format(
            self.fileLocation,
            len(self.variants),
            sum([len(v) for v in self.variants.values()])
        )

if __name__ == "__main__":
    pass
