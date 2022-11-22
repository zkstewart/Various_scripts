#! python3
# ZS_VCFIO.py
# Contains the VCF class which encapsulated a relatively
# simple dictionary data structure. Use this for parsing
# and filtering a VCF in a performant but high-memory
# use way.

import re, os, codecs

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            pass
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                pass
            return "utf-16"
        except UnicodeDecodeError:
            print(f"VCF class can't tell what codec '{fileName}' is!!")

class VCF:
    def __init__(self, file_location):
        assert os.path.isfile(file_location), \
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
        with open(self.fileLocation, "r", encoding=get_codec(self.fileLocation)) as fileIn:
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
                        "INFO": info,
                        "FORMAT": format
                    }
                    
                    # Parse sample details according to FORMAT 
                    for i in range(len(self.samples)):
                        sampleFormat = format.split(":")
                        sampleDetails = sl[9+i].split(":")
                        
                        sampleDict = {}
                        for x in range(len(sampleFormat)):
                            sampleDict[sampleFormat[x]] = sampleDetails[x]
                        
                        # Store result
                        self.variants[chrom][pos][self.samples[i]] = sampleDict # sampleIDs[i] == sampleName
                    
                    # If DP is listed in info but not per-sample, impute it now
                    if "DP=" in info:
                        infoDP = [x.split("=")[1] for x in info.split(";") if x.startswith("DP=")][0]
                    
                        if "DP" not in sampleDict and "AD" in sampleDict: # most recently used sampleDict is fine
                            self._impute_dp(chrom, pos, infoDP)
    
    def _impute_dp(self, chrom, pos, infoDP):
        '''
        Simply put, imputes DP for each sample if it's not found for any samples.
        '''
        sampleADs = [self.variants[chrom][pos][id]["AD"] for id in self.samples]
        adsRatio = [sum(map(int, ad.split(","))) / int(infoDP) for ad in sampleADs]
        dpsImpute = [int(adR*int(infoDP)) for adR in adsRatio]
        for i in range(len(self.samples)):
            self.variants[chrom][pos][self.samples[i]]["DP"] = str(dpsImpute[i])
    
    def __getitem__(self, key):
        return self.variants[key]
    
    def __setitem__(self, key, item):
        self.variants[key] = item
    
    def __len__(self):
        return len(self.variants)
    
    def __delitem__(self, key):
        if key not in self.variants:
            raise ValueError(f"'{key}' not found in this VCF object")
        else:
            del self.variants[key]
    
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
        # return "<VCF object;file='{0}';num_variants={1};{2}>".format(
        #     self.fileLocation,
        #     len(self.variants),
        #     ";".join(["num_{0}={1}".format(key, len(self.types[key])) for key in self.types.keys()])
        # )
        return "<VCF object;file='{0}';num_contigs={1};num_variants=TBD>".format(
            self.fileLocation,
            len(self.variants),
            len(self.variants.values())
        )

if __name__ == "__main__":
    pass
