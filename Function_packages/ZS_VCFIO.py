#! python3
# ZS_VCFIO.py
# Contains the VCF class which encapsulated a relatively
# simple dictionary data structure. Use this for parsing
# and filtering a VCF in a performant but high-memory
# use way.

import os, codecs, gzip
from contextlib import contextmanager

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

if __name__ == "__main__":
    pass
