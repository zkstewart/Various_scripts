#! python3
# ZS_BAMIO.py
# Implemented as a bamnostic.AlignmentFile Class with extra
# functions added for my own purposes.

import os, pickle, math
import bamnostic as bs
import numpy as np
import pandas as pd

def get_chunking_points(numberToChunk, chunkSize):
    '''
    This is a general purpose function to take in a number of "things"
    that you want to chunk, and find out how to chunk them evenly.
    
    Params:
        numberToChunk -- an integer value, possibly derived from a list length as example.
        chunkSize -- an integer value for the desired number of things per chunk.
    '''
    assert isinstance(numberToChunk, int)
    assert isinstance(chunkSize, int)
    if numberToChunk <= chunkSize:
        raise Exception("Chunking only valid if chunkSize is smaller than numberToChunk")
    
    numChunks = int(numberToChunk / chunkSize)
    rawNum = numberToChunk / numChunks # This line is more relevant in the multithreading code I took this from, but it's okay to just leave it.
    numRoundedUp = round((rawNum % 1) * numChunks, 0) # By taking the decimal place and multiplying it by the num of chunks, we can figure out how many chunks need to be rounded up
    
    chunkPoints = []
    ongoingCount = 0
    for i in range(numChunks):
        if i+1 <= numRoundedUp: # ngl I don't remember why this is needed; I'm borrowing this code from something I wrote a while back
            chunkPoints.append(math.ceil(rawNum) + ongoingCount) # Round up the rawNum, and also add our ongoingCount which corresponds to the number of things already put into a chunk
            ongoingCount += math.ceil(rawNum)
        else:
            chunkPoints.append(math.floor(rawNum) + ongoingCount)
            ongoingCount += math.floor(rawNum)
        if ongoingCount >= numberToChunk: # Without this check, if we have more chunks than things to chunk, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
            break  # This doesn't actually affect program function, but for aesthetic reasons and for clarity of how this function works, I prevent this from occurring.
    
    return chunkPoints

class BAM(bs.AlignmentFile):
    def __init__(self, file_location):
        super().__init__(file_location, "rb")
        
        self.fileLocation = file_location
        self.coverage = None
        self.coverage_histogram = None
        self.gbc = None
    
    def set_blank_coverage(self):
        '''
        This function sets the .coverage field of this instance to be a
        dictionary with keys equal to contig names, and values being
        numpy arrays of 0's equal in length to the length of contigs
        specified within the BAM header.
        '''
        assert "SQ" in self.header.keys(), \
            "SQ not specified in BAM header; I can't tell what contigs it references!"
        
        self.coverage = {}
        raiseWarning = False
        for sqValue in self.header["SQ"]:
            try:
                self.coverage[sqValue["SN"]] = np.zeros(sqValue["LN"])
            except:
                raiseWarning = True
        if raiseWarning:
            print("Warning: One or more header SQ values lacked a SN and/or LN value") 
    
    def compute_coverage(self, gff3Obj=None):
        '''
        Computes coverage for contigs at each position. It will populate the
        .coverage field with a dictionary with structure like:
        {
            contigID_1: np.array([pos1CovInt, ...]),
            contigID_2: ...,
            ...
        }
        
        When this function has completed running, the result will be pickled in the
        file with name ${self.fileLocation}.cov.pkl. If a file exists with this name,
        it will load it in assuming it is an already pickled result, since this function
        otherwise takes a very long time.
        
        If you provide a ZS_GFF3IO.GFF3 object, this function will operate in a behaviour
        assuming the BAM file has been mapped to chromosomal contigs, but we actually
        want coverage values for genic regions. It will use the GFF3 object to do this,
        obtaining counts that cover CDS regions specifically.
        
        Parameters:
            gff3Obj -- a ZS_GFF3IO object containing an indexed GFF3 which corresponds
                       to the BAM file.
        '''
        FLAGS_TO_SKIP = [
            (2048, "supplementary alignment"), (1024, "read is PCR or optical duplicate"),
            (512, "read fails platform/vendor quality checks"), (256, "not primary alignment"),
            (4, "read unmapped")
        ]
        
        # Figure out if this has already been computed and picked
        pickleFileName = "{0}.cov.pkl".format(self.fileLocation)
        pickleExists = os.path.isfile(pickleFileName)
        
        # If it does not, compute it now
        if not pickleExists:
            self.set_blank_coverage()
            for read in self:
                # Check read flags to see if we want to count it
                flags = bs.utils.flag_decode(read.flag)
                shouldSkip = any([skipFlag in flags for skipFlag in FLAGS_TO_SKIP])
                if shouldSkip:
                    continue
                
                # Get the aligned base coordinates from interpretation of the CIGAR
                prevPos = None
                coords = []
                for base, refPos in bs.utils.cigar_alignment(read.seq, read.cigarstring):
                    # Handle first iteration
                    if prevPos == None and base != ".":
                        coords.append([refPos, None])
                        prevPos = refPos
                    else:
                        # Handle read deletions
                        if coords[-1][1] == None and base == ".":
                            coords[-1][1] = prevPos # go back to last found position
                        # Handle reference deletions
                        elif coords[-1][1] == None and refPos != prevPos + 1:
                            coords[-1][1] = prevPos # go back to last found position
                        # Start a new coord range if we filled in the last one
                        if coords[-1][1] != None and base != ".":
                            coords.append([refPos, None])
                        prevPos = refPos
                if coords[-1][1] == None:
                    coords[-1][1] = prevPos
                
                # Convert alignment coordinates to reference coordinates
                coords = [[coord[0] + read.reference_start, coord[1] + read.reference_start] for coord in coords]
                
                # Modify coverages value
                for start, end in coords:
                    self.coverage[read.reference_name][start:end + 1] += 1 # offset range not considering the last position
            
            # If a GFF3 was provided, get counts to gene features
            if gff3Obj != None:
                geneCoverage = {}
                for geneFeature in gff3Obj.types["gene"]:
                    contig = geneFeature.contig
                    cdsFeatures = [feature for feature in geneFeature.retrieve_all_children() if feature.type == "CDS"]
                    assert contig in self.coverage, \
                        "'{0}' contig is missing from original BAM alignment; GFF3 doesn't match?".format(contig)
                    
                    # Figure out the exon coordinates for this gene
                    "This will work on an open-shut algorithm, where when all opens are shut, the exon has ended"
                    coords = [["{0}-start".format(feature.start), "{0}-end".format(feature.end)] for feature in cdsFeatures]
                    coords = [value for coord in coords for value in coord] # flatten list
                    coords.sort(key = lambda x: int(x.split("-")[0]))
                    
                    mergedCoords = []
                    opens = 0
                    for coord in coords:
                        if "start" in coord:
                            if opens == 0:
                                mergedCoords.append([int(coord.split("-")[0]), None])
                            opens += 1
                        if "end" in coord:
                            opens -= 1
                            if opens == 0:
                                mergedCoords[-1][1] = int(coord.split("-")[0])
                    
                    # Extract counts for the exon regions
                    geneCovArray = np.concatenate(
                        [self.coverage[contig][coord[0]-1:coord[1]] for coord in mergedCoords] # convert 1-based GFF3 coords to 0-based
                    )
                    
                    # Get the strand for this gene & flip if needed
                    if geneFeature.strand == "-":
                        geneCovArray = geneCovArray[::-1]
                    
                    # Store it in the dictionary
                    geneCoverage[geneFeature.ID] = geneCovArray
                
                # Overwrite chromosome coverage with gene coverage
                self.coverage = geneCoverage
            
            # Store as pickle
            pickle.dump(self.coverage, open(pickleFileName, "wb"))
        # If it does exist, load it and store it in this object
        else:
            self.coverage = pickle.load(open(pickleFileName, "rb"))
    
    def summarise_coverage_into_histogram(self, binPctSize=10, gff3Obj=None, mappingDict=None):
        '''
        After running the compute_coverage() method, this function will
        summarise coverage values into per-contig histograms. Each bin
        will be composed of a set percentage of the contig's length.
        
        Providing the GFF3 object during compute_coverage() handles cases
        where you have mapped to chromosomes and you want to extract gene
        counts. Providing it HERE instead handles cases where you've mapped
        to mRNA transcripts where the CDS may be in the +ve or -ve strand.
        As such, providng it here will flip the counts where appropriate.
        
        If the IDs between your GFF don't match the BAM, you will need to
        provide the mappingDict as well. This should relate the keys in
        self.coverage_histogram to the gene/mRNA IDs in the GFF3 object.
        When multiple mappings exist e.g., multiple mRNAs point to the same
        locus ID, we will combine the counts in this function.
        
        There's a lot going on here admittedly, so think before you act.
        
        Parameters:
            binPctSize -- an integer or float value that dictates the size of histogram
                          bins as a percentage from 1 -> 50. This value must provide a 
                          whole number when 
            gff3Obj -- a ZS_GFF3IO object containing an indexed GFF3 which corresponds
                       to the BAM file.
        
        Sets:
            .coverage_histogram -- a dictionary with structure like:
                                   {
                                       contigID_1: np.array([pos1Cov, pos2Cov, ...]),
                                       ...
                                   }
        '''
        assert self.coverage is not None, \
            "Run compute_coverage() before calling this method!"
        try:
            assert binPctSize < 100 and 100 % binPctSize == 0, \
                "binPctSize must be a value less than 100 that cleanly divides 100"
        except:
            "binPctSize must be a float or integer value"
        
        self.coverage_histogram = {}
        for contig, coverage in self.coverage.items():
            # Figure out the histogram bin boundaries
            contigLength = len(coverage)
            binSize = int(contigLength / binPctSize)
            binBoundaries = get_chunking_points(contigLength, binSize)
            
            # Summarise coverage per bin
            self.coverage_histogram[contig] = np.zeros(len(binBoundaries))
            for i in range(len(binBoundaries)):
                boundaryStart = 0 if i == 0 else binBoundaries[i-1] + 1
                boundaryEnd = binBoundaries[i]
                self.coverage_histogram[contig][i] = np.sum(self.coverage[contig][boundaryStart:boundaryEnd + 1]) # +1 to make range inclusive
        
        # If a mappingDict is provided, modify keys and merge counts where needed
        if mappingDict != None:
            new_coverage_histogram = {}
            for contig, coverage in self.coverage_histogram.items():
                if contig not in mappingDict:
                    print("'{0}' contig from BAM is not in the mapping file; summarise_coverage_into_histogram() is dropping this".format(contig))
                    continue
                
                new_coverage_histogram.setdefault(mappingDict[contig], np.zeros(coverage.size))
                new_coverage_histogram[mappingDict[contig]] += coverage
            self.coverage_histogram = new_coverage_histogram
        
        # If GFF3 object is provided, flip any counts where needed
        if gff3Obj != None:
            for contig, coverage in self.coverage_histogram.items():
                assert contig in gff3Obj, \
                    "'{0}' contig not found in GFF3; summarise_coverage_into_histogram() failed!".format(contig)
                
                try:
                    strand = gff3Obj[contig].strand
                except:
                    "'{0}' contig in GFF3 does not have a strand field; summarise_coverage_into_histogram() failed!".format(contig)
                
                if strand == "-":
                    self.coverage_histogram[contig] = coverage[::-1]
    
    def qc_genebody_coverage(self):
        '''
        QC FUNCTION: This method will compute the genebody coverage statistic.
        By default, contigs are treated as if they are gene model CDS, and the
        coverage histogram will be compared as-is.
        
        Sets:
            .gbc -- the GeneBody Coverage for this instance as a Pandas DataFrame where
                    row indices == gene names, and column values = min-max normalised (by
                    row) genebody coverage ratios.
        '''
        assert self.coverage_histogram != None, \
            "Run summarise_coverage_into_histogram() before calling this method!"
        
        # Establish pandas dataframe
        df = pd.DataFrame(self.coverage_histogram)
        df = df.T # transpose so contigs are rows and columns are bins
        
        # Drop rows where all values are 0
        df = df.loc[(df != 0).any(axis=1)]
        
        # Min-max normalise each row
        df = df.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)
        
        # Store as instance field
        self.gbc = df
    
    def __repr__(self):
        return "<BAM object;file='{0}'".format(
            self.fileLocation
        )

if __name__ == "__main__":
    pass
