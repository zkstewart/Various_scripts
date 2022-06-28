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
    
    def compute_coverage(self):
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
                
                # Extract relevant alignment details
                contig = read.reference_name
                start = read.reference_start # 0-based
                end = read.reference_end # 0-based
                
                # Modify coverages value
                self.coverage[contig][start:end + 1] += 1 # offset range not considering the last position
            # Store as pickle
            pickle.dump(self.coverage, open(pickleFileName, "wb"))
        # If it does exist, load it and store it in this object
        else:
            self.coverage = pickle.load(open(pickleFileName, "rb"))
    
    def summarise_coverage_into_histogram(self, binPctSize=10):
        '''
        After running the compute_coverage() method, this function will
        summarise coverage values into per-contig histograms. Each bin
        will be composed of a set percentage of the contig's length.
        
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
    
    def qc_genebody_coverage(self, gff3Obj=None, mappingDict=None):
        '''
        QC FUNCTION: This method will compute the genebody coverage statistic.
        By default, contigs are treated as if they are gene model CDS, and the
        coverage histogram will be compared as-is.
        
        If the BAM file contains alignments to genomic contigs, providing a ZS_GFF3IO
        object will compute coverage statistics over the genic regions (exons) found
        within the chromosomal contigs. It will also ensure that coverage statistics
        are computed considering the gene's strand i.e., all genes will be considered
        as if they are +ve stranded going from 5' to 3'.
        
        Parameters:
            gff3Obj -- a ZS_GFF3IO object containing an indexed GFF3 which corresponds
                       to the BAM file either directly or through a mappingDict intermediary.
            mappingDict -- a dictionary containing mappings between keys in the BAM file
                           to keys in the GFF3 object e.g.,
                           {
                               bam_ID1: GFF3_ID1,
                               ...
                           }
            outputPlotName -- a string indicating the location to write the plot to. If
                              this file exists or the parameter is not specified, it will
                              be automatically created with the format ""
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
        
        # If a mappingDict was provided, see if we should drop any rows
        "Mappings that don't exist are interpreted to be mappings that SHOULDN'T exist"
        rowsToDrop = []
        if mappingDict != None:
            for index in df.index:
                if index not in mappingDict:
                    print("'{0}' cannot be found within the mappingDict; this row will be dropped".format(index))
                    rowsToDrop.append(index)
        if rowsToDrop != []:
            df = df.drop(rowsToDrop)
        
        # If a GFF3 was provided, get counts to gene features
        ### TBD!!
        rowsToDrop = []
        if gff3Obj != None:
            for index, row in df.iterrows():
                # If mappingDict is provided, map the BAM ID to the GFF3 ID
                if mappingDict != None:
                    featureID = mappingDict[index] # this is guaranteed to exist because of the above rowsToDrop
                else:
                    featureID = index
                
                # Get the strand for this gene
                try:
                    strand = gff3Obj[featureID].strand
                except:
                    print("'{0}' cannot be found within the GFF3; qc_genebody_coverage failed".format(featureID))
                    return
                
                # Flip if needed
                if strand == "-":
                    reversedRow = list(row.array[::-1])
                    for i in range(len(row)):
                        df.at[index,i] = reversedRow[i]
        
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
