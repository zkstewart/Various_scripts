#! python3
# ZS_GFF3IO.py
# Contains the GFF3 class and its associated Feature
# Class. Can be used for (memory-heavy) GFF3 file parsing

import re, sys, os
import pandas as pd
from collections import OrderedDict
from ncls import NCLS

sys.path.append(os.path.dirname(__file__))
from ZS_SeqIO import FastASeq

class Feature:
    '''
    A Feature object is intended to be quite flexible. Set its attributes to
    whatever you want. The purpose of it is to track children for which it
    provides a useful method to chase down all the kids.
    
    However, it is also designed first-and-foremost to be a GFF3 Feature.
    That means it will behave in a way that facilitates that purpose.
    
    An important thing to consider is how this class works with children.
    There's no observer pattern here, so the state of the Feature when
    added as a child will dictate how it's loaded here. Changing the
    child feature afterwards will NOT update any field setting in the
    parent. So, fully specify your child Feature before adding it!
    '''
    def __init__(self):
        self.children = []
        self.types = {}
        self.isFeature = True
    
    def add_attributes(self, dict):
        '''
        This method will associate a provided dictionary
        to be fields within this Feature instance
        
        Parameters:
            dict -- any dictionary with key: value pairs
        '''
        for key, value in dict.items():
            self.__dict__[key] = value
    
    def add_child(self, childFeature):
        '''
        This method is intended to add child Features to this Feature.
        In the process of adding the child feature, it will create
        fields as necessary for this Feature to house the child.
        The field key will be determined by the "type" field of the
        added Feature.
        
        So for example, if _this_ is a gene Feature, adding an mRNA
        feature to it will create a self.mRNA key, which points to
        a list containing the added mRNA child.
        
        Lastly, it also adds things to the self.types dict akin to
        how the GFF3 class .types dict works.
        
        Parameters:
            childFeature -- should be a Feature object, but any object
                            will be accepted.
        '''
        self.children.append(childFeature)
        try:
            self.__dict__.setdefault(childFeature.type, [])
            self.__dict__[childFeature.type].append(childFeature)
            
            self.types.setdefault(childFeature.type, [])
            self.types[childFeature.type].append(childFeature)
        except:
            pass
    
    def retrieve_child(self, childID):
        '''
        Parameters:
            childID -- a string indicating the ID field of the child to return
        Returns:
            child -- the child object/Feature that was found, OR a None value
                     when the child does not exist.
        '''
        childList = self.retrieve_all_children()
        for child in childList:
            try:
                if child.ID == childID:
                    return child
            except:
                pass
        return None
    
    def retrieve_all_children(self):
        childList = []
        for child in self.children:
            childList.append(child)
            childList += child.retrieve_all_children()
        return childList
    
    def __getitem__(self, key):
        return self.retrieve_child(key)
    
    def __repr__(self):
        reprPairs = []
        attrsToShow = ["ID", "type", "Parent", "coords"]
        
        for attr in attrsToShow:
            try:
                reprPairs.append("{0}={1}".format(attr, self.__dict__[attr]))
            except:
                pass
        
        return "<{0}>".format(";".join(reprPairs))

class GFF3:
    def __init__(self, file_location, strict_parse=True):
        self.fileLocation = file_location
        self.features = OrderedDict()
        self.types = {}
        self.contigs = set()
        self.parentTypes = set() # tells us which feature types to expect as being parents
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        
        self.isGFF3 = True
        self.parse_gff3(strictParse=strict_parse)
    
    @staticmethod
    def make_feature_case_appropriate(featureType):
        if featureType.lower() == "gene":
            return "gene"
        elif featureType.lower() == "mrna":
            return "mRNA"
        elif featureType.lower() == "exon":
            return "exon"
        elif featureType.lower() == "cds":
            return "CDS"
        else:
            return featureType
    
    @staticmethod
    def longest_isoform(geneFeature):
        '''
        We pick out the representative gene based on length. If length is identical,
        we'll end up picking the entry listed first in the gff3 file since our > condition
        won't be met. I doubt this will happen much or at all though.
        '''
        assert hasattr(geneFeature, "mRNA"), \
            "Longest isoform finding can only occur on features that have .mRNA children"
        
        longestMrna = [None, 0]
        for mrnaFeature in geneFeature.mRNA:
            if hasattr(mrnaFeature, "CDS"):
                featType = "CDS"
            else:
                featType = "exon"
            
            mrnaLen = 0
            for subFeature in mrnaFeature.__dict__[featType]:
                mrnaLen += (subFeature.end - subFeature.start + 1)
                
            if mrnaLen > longestMrna[1]:
                longestMrna = [mrnaFeature, mrnaLen]
        return longestMrna[0]
    
    def parse_gff3(self, strictParse=True):
        # Gene object loop
        lineCount = 0
        with open(self.fileLocation, 'r') as fileIn:
            for line in fileIn:
                lineCount += 1
                line = line.replace('\r', '')
                
                # Skip filler and comment lines
                if line == "\n" or line.startswith("#"):
                    continue
                
                # Extract information from this line
                try:
                    contig, source, featureType, start, end, \
                        score, strand, frame, attributes \
                        = line.rstrip('\t\n').split('\t')
                except:
                    if strictParse == True:
                        raise ValueError(
                            f"Line #{lineCount} (\"{line}\") does not meet GFF3 standards; parsing failed"
                        )
                    else:
                        print(
                            f"Line #{lineCount} (\"{line}\") does not meet GFF3 standards;" +
                            " strict parsing is disabled, so we'll just continue and hope for the best"
                        )
                        continue
                
                splitAttributes = []
                for a in attributes.split("="):
                    if ";" in a:
                        splitAttributes += a.rsplit(";", maxsplit=1)
                    else:
                        splitAttributes.append(a)
                attributesDict = {splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)}
                
                self.contigs.add(contig)
                
                # Ensure case conformity
                featureType = GFF3.make_feature_case_appropriate(featureType)
                
                # Skip un-indexable features
                if 'ID' not in attributesDict and featureType.lower() != "cds": # see the human genome GFF3 biological_region values for why this is necessary
                    continue
                
                # Handle parent-level features
                if 'Parent' not in attributesDict: # If no Parent field this should BE the parent
                    featureID = attributesDict["ID"]
                    
                    # End parsing if duplicate ID is found (strictParse is True)
                    if strictParse == True:
                        assert featureID not in self.features, \
                            "'{0}' feature occurs twice indicating poorly formatted file; \
                            for debugging, line #{1} == {2}; parsing will stop now".format(featureID, lineCount, line)
                    # Skip parsing if duplicate ID is found (strictParse is False)
                    else:
                        if featureID in self.features:
                            print(
                                "'{0}' feature occurs more than once indicating poorly formatted file; \
                                strict parsing is disabled, so we will just continue and hope for the best".format(featureID)
                            )
                            continue
                    
                    # Create feature and index it
                    feature = Feature()
                    self.features[featureID] = feature
                    self.types.setdefault(featureType, [])
                    self.types[featureType].append(feature)
                    self.parentTypes.add(featureType)
                    
                    # Populate feature with details
                    feature.add_attributes(attributesDict)
                    feature.add_attributes({
                        "contig": contig, "source": source, "type": featureType,
                        "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                        "score": score, "strand": strand, "frame": frame
                    })
                
                # Handle subfeatures
                else:
                    parents = attributesDict["Parent"].split(',')
                    
                    # Loop through parents and associate feature to them
                    for parentID in parents:
                        # Flexibly obtain a feature ID for normal features and for CDS that may lack IDs
                        try:
                            featureID = attributesDict["ID"]
                        except:
                            if featureType.lower() != "cds":
                                raise AttributeError(
                                    "'{0}' feature lacks an ID attribute and hence cannot be indexed; \
                                    for debugging, line #{1} == {2}; parsing will stop now".format(featureType, lineCount, line)
                                )
                            else:
                                featureID = "{0}.cds".format(parentID)
                        
                        # End parsing if parent doesn't exist (strictParse is True)
                        if strictParse == True:
                            assert parentID in self.features, \
                                "'{0}' feature points to a non-existing parent '{1}' at the time of parsing; \
                                for debugging, line #{2} == {3}; parsing will stop now".format(featureID, parentID, lineCount, line)
                        else:
                            if parentID not in self.features:
                                print("'{0}' feature points to a non-existing parent '{1}' at the time of parsing; ".format(featureID, parentID) + 
                                "strict parsing is disabled, so we will just continue and hope for the best")
                                continue
                        
                        # End parsing if duplicate ID is found
                        "strictParse won't negate this problem since it's simply unacceptable!"
                        if featureType.lower() != "cds": # CDS features are the ONLY feature allowed to have non-unique IDs
                            assert self.features[parentID].retrieve_child(featureID) == None, \
                                "'{0}' feature is associated to a parent '{1}' more than once; \
                                for debugging, line #{2} == {3}; parsing will stop now".format(featureID, parentID, lineCount, line)
                        
                        # Create feature and index it
                        feature = Feature()
                        if featureType.lower() != "cds": # since CDS features aren't guaranteed to have unique IDs, there's no point indexing them
                            self.features[featureID] = feature
                        self.types.setdefault(featureType, [])
                        self.types[featureType].append(feature)
                        
                        # Populate feature with details
                        feature.add_attributes(attributesDict)
                        feature.add_attributes({
                            "contig": contig, "source": source, "type": featureType,
                            "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                            "score": score, "strand": strand, "frame": frame
                        })
                        
                        # Add it as a child of the parent
                        "It's important to fully-specify the child before running add_child()"
                        self.features[parentID].add_child(feature)
            
            # Generate shortcut fields
            self.gene_values = self.types["gene"]
            try:
                self.mrna_values = self.types["mRNA"]
            except:
                self.mrna_values = None
            
            # Sort contigs
            self.contigs = list(self.contigs)
            try:
                self.contigs.sort(key = lambda x: list(map(int, re.findall(r'\d+', x)))) # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
            except:
                self.contigs.sort()
    
    def sort_CDS(self):
        '''
        This method will take any parent feature that indexes CDS features, and ensures
        that the CDS children are sorted in the way we would normally expect them to be.
        For example, a +ve stranded gene should have CDS sorted so that the first entry
        in a list is the "leftmost" on the chromosome (lesser coordinate). Alternatively,
        a -ve stranded gene should have CDS sorted so that the first entry in a list is
        the "rightmost" on the chromosome (greatest coordinate). The logic behind this
        is to have CDS features sorted in 5'->3' order, rather than focusing on the
        absolute chromosomal coordinates.
        '''
        assert "CDS" in self.types, \
            "There are no CDS features in this GFF3; sorting is irrelevant"
        
        for parentID in set([x.Parent for x in self.types["CDS"]]):
            parentFeature = self[parentID]
            if parentFeature.strand == "+":
                parentFeature.CDS.sort(key = lambda x: x.start)
            elif parentFeature.strand == "-":
                parentFeature.CDS.sort(key = lambda x: -x.end)
            else:
                raise ValueError(f"'{parentID} feature has CDS children but strand '{parentFeature.strand}' is unrecognised.")
    
    def sort_exon(self):
        '''
        Akin to sort_CDS(), this function does the same but for exon-indexing features.
        The differences comes in how we approach features without a proper strand
        being noted (+ve or -ve). In these cases, we'll just skip over the feature since,
        without strand info, we can't know how it's "supposed" to be sorted.
        '''
        assert "exon" in self.types, \
            "There are no exon features in this GFF3; sorting is irrelevant"
        
        for parentID in set([x.Parent for x in self.types["exon"]]):
            parentFeature = self[parentID]
            if parentFeature.strand == "+":
                parentFeature.exon.sort(key = lambda x: x.start)
            elif parentFeature.strand == "-":
                parentFeature.exon.sort(key = lambda x: -x.end)
            else:
                continue
    
    def infer_UTRs(self):
        '''
        This method will associate UTR Features to the GFF3 object when they can be
        inferred to exist. Specifically, where there's exons uncovered by CDS annotations,
        we can surmise that a UTR exists in that region.
        '''
        raise NotImplementedError()
    
    def create_ncls_index(self, typeToIndex="mRNA"):
        '''
        Creates an indexed NCLS structure that can be used to find range overlaps
        for the feature types of interest.
        
        Associates the created index to the .ncls field of this object instance.
        A hidden ._nclsIndex dictionary links the ncls indices to feature objects.
        
        Parameters:
            typeToIndex -- a string (case-sensitive) indicating the entry type
                           to index.
        '''
        assert typeToIndex in self.types, \
            "'{0}' not found as a feature type within the parsed GFF3 ('{1}')".format(typeToIndex, self.fileLocation)
        
        nclsIndex = {}
        starts, ends, ids = [], [], []
        
        # Add features of the specified type to our NCLS structure
        ongoingCount = 0
        for feature in self.types[typeToIndex]:
            starts.append(feature.start)
            ends.append(feature.end + 1) # NCLS indexes 0-based like a range so +1 to make this more logically compliant with gff3 1-based system
            ids.append(ongoingCount)
            nclsIndex[ongoingCount] = feature
            ongoingCount += 1
        
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        
        # Associate it to this instance
        self.ncls = ncls
        self._nclsType = typeToIndex
        self._nclsIndex = nclsIndex
    
    def ncls_finder(self, start, stop, field, value):
        '''
        Queries the NCLS structure to find Features that exist within the given
        start->stop range. Specifying the field and value will narrow results
        to only those that have a Feature .field with an equal (==) value.
        
        Parameters:
            start -- an integer indicating the start position of the feature to check
                     for overlaps
            end -- an integer indicating the end positon of the feature to check for
                   overlaps; this should be 1-based in GFF3 style e.g., a first
                   position of a feature would be start=1, end=1.
            field -- a string (case-sensitive) indicating the field of the Feature
                     object that we want to check. For example, if you want to find
                     features that overlap the given start->stop range on the contig
                     "X", you'd provide "contig" as the field so this function knows
                     to check the Feature.contig field for the value of "X".
            value -- a string (case-sensitive) indicating the value of the Feature
                     field we want to find. As in the above example, if you want to
                     find the value "X" within the .contig field, you'd provide "X" as
                     the value here.
        Returns:
            features -- a list containing Features that overlap the specified range.
                        These Features are NOT deepcopied, so handle them carefully.
        '''
        assert self.ncls != None and self._nclsIndex != None, \
            "Run create_ncls_index before you call this method!"
        
        overlaps = self.ncls.find_overlap(start, stop+1) # Although our ncls is already 1-based, find_overlap acts as a range. We need to +1 to keep everything logically 1-based.
        
        features = []
        for result in overlaps: # result == [start, end, index]
            feature = self._nclsIndex[result[2]]
            if feature.__dict__[field] == value:
                features.append(feature)
        
        # Return list
        return features
    
    @staticmethod
    def _get_contig_as_FastASeq(contigSequences, contigID):
        '''
        Hidden helper method of retrieve_sequence_from_FASTA(). Intended to work with
        many different object types i.e., a FASTA, FastASeq, or Biopython sequence.
        '''
        # Handle Biopython dict type
        if type(contigSequences).__name__ == "dict":
            assert contigID in contigSequences, \
                "'{0}' does not exist within the dictionary contig sequences object".format(contigID)
            assert type(contigSequences[contigID]).__name__ == "SeqRecord", \
                "contigSequences is a dictionary, but values aren't Biopython SeqRecords? Can't handle."
            
            sequence = str(contigSequences[contigID].seq)
        
        # Handle ZS_SeqIO.FASTA type
        elif hasattr(contigSequences, "isFASTA") and contigSequences.isFASTA is True:
            try:
                sequence = contigSequences[contigID].seq
            except:
                raise TypeError("'{0}' does not exist within the FASTA object".format(contigID))
        
        # Handle ZS_SeqIO.FastASeq type
        elif hasattr(contigSequences, "isFastASeq") and contigSequences.isFastASeq is True:
            assert contigID == contigSequences.id or contigID == contigSequences.alt, \
                "'{0}' contig does not match the provided FastASeq object".format(contigID)
            sequence = contigSequences.seq
        
        # Throw error for unhandled types
        else:
            raise TypeError("Unrecognised data type '{0}' for contig string retrieval".format(type(contigSequences).__name__))
        
        # Create a FastASeq object to return
        outputSeq = FastASeq(id=contigID, seq=sequence)
        return outputSeq
    
    @staticmethod
    def _get_contig_as_string(contigSequences, contigID):
        '''
        Hidden helper method of retrieve_sequence_from_FASTA(). Intended to work with
        many different object types i.e., a FASTA, FastASeq, or Biopython sequence.
        '''
        # Handle Biopython dict type
        if type(contigSequences).__name__ == "dict":
            assert contigID in contigSequences, \
                "'{0}' does not exist within the dictionary contig sequences object".format(contigID)
            assert type(contigSequences[contigID]).__name__ == "SeqRecord", \
                "contigSequences is a dictionary, but values aren't Biopython SeqRecords? Can't handle."
            
            sequence = str(contigSequences[contigID].seq)
        
        # Handle ZS_SeqIO.FASTA type
        elif hasattr(contigSequences, "isFASTA") and contigSequences.isFASTA is True:
            try:
                sequence = contigSequences[contigID].seq
            except:
                raise TypeError("'{0}' does not exist within the FASTA object".format(contigID))
        
        # Handle ZS_SeqIO.FastASeq type
        elif hasattr(contigSequences, "isFastASeq") and contigSequences.isFastASeq is True:
            assert contigID == contigSequences.id or contigID == contigSequences.alt, \
                "'{0}' contig does not match the provided FastASeq object".format(contigID)
            sequence = contigSequences.seq
        
        # Throw error for unhandled types
        else:
            raise TypeError("Unrecognised data type '{0}' for contig string retrieval".format(type(contigSequences).__name__))
        
        # Return the sequence as a string
        return sequence
    
    @staticmethod
    def _get_feature_coords(feature, exonOrCDS):
        '''
        Hidden function of retrieve_sequence_from_FASTA() to retrieve sorted
        coordinates lists for the exon or CDS features associated with what
        should be an mRNA feature. If it's not it'll probably crash, hence why
        this is a private method since I know exactly how it'll be used.
        '''
        coords = [f.coords for f in feature.__dict__[exonOrCDS]]
        frames = [f.frame for f in feature.__dict__[exonOrCDS]]
        forSorting = list(zip(coords, frames))
        
        if feature.strand == '+':
            forSorting.sort(key = lambda x: (int(x[0][0]), int(x[0][1])))
        else:
            forSorting.sort(key = lambda x: (-int(x[0][0]), -int(x[0][1])))
        
        coords = [c for c, f in forSorting]
        frames = [f for c, f in forSorting]
        
        return coords, frames
    
    def _get_artifacts_as_vcfDict(self):
        '''
        Hidden function which enables quick retrieval of recognised insertion, deletions,
        and substitution artifact entries from this GFF3. It will return an object
        VCF-like dictionary format which can be used for editing a sequence object.
        
        Returns:
            vcfDict -- a dictionary with structure like:
                       {
                           contigID1: [
                               [[start1, end1], artifact_type, residue],
                               ...
                           ],
                           ...
                       }
        '''
        ACCEPTED_INDEL_SUB_TYPES = ["deletion_artifact", "substitution_artifact", "insertion_artifact"]
        
        vcfDict = {}
        for artifactType in ACCEPTED_INDEL_SUB_TYPES:
            if artifactType in self.types:
                for artifactFeature in self.types[artifactType]:
                    vcfDict.setdefault(artifactFeature.contig, [])
                    residue = "." if not hasattr(artifactFeature, "residues") else artifactFeature.residues
                    vcfDict[artifactFeature.contig].append([
                        artifactFeature.coords, artifactFeature.type, residue
                    ])
        
        for value in vcfDict.values():
            value.sort(key = lambda x: ([-x[0][0], -x[0][1]]))
        
        return vcfDict
    
    def retrieve_coords(self, feature_or_featureID, sequenceType):
        '''
        Using this GFF3 instance, retrieve the exon or CDS coordinates for the corresponding
        sequenceID. Note that this function expects the provided sequenceID to directly
        point to a feature that contains CDS and/or exon values. It used to be able to accept
        gene features and return all the subfeature values, but this proved to be unwieldy.
        
        Parameters:
            feature_or_featureID -- a string corresponding to a feature within this GFF3 object.
                                    or just a feature object in general
            sequenceType -- a string corresponding to the type of sequence to retrieve
                            i.e., in the list ["CDS", "exon"]
        Returns:
            featureCoords -- a list of lists with format like:
                             [
                                 [start_1, end_1],
                                 [start_2, end_2],
                                 ...
                             ]
            startingFrames -- a list of integers indicating what the starting frame should
                              be in any translations (if applicable)
            featureTypes -- a list of strings indicating what type of feature's details
                            have been returned e.g., "mRNA" or "lnc_RNA".
        '''
        VALID_TYPES = ["cds", "exon"]
        assert sequenceType.lower() in VALID_TYPES, \
            "'{0}' is not recognised as a valid sequenceType; should be in list {1}".format(sequenceType.lower(), VALID_TYPES)
        
        # Figure out if we're handling a feature or featureID
        if isinstance(feature_or_featureID, str):
            featureID = feature_or_featureID
            assert featureID in self.features, \
                "'{0}' is not recognised as a feature within this GFF3".format(featureID)
            feature = self.features[featureID]
        elif type(feature_or_featureID).__name__ == "Feature" \
            or type(feature_or_featureID).__name__ == "ZS_GFF3IO.Feature" \
            or (hasattr(feature_or_featureID, "isFeature") and feature_or_featureID.isFeature is True):
                feature = feature_or_featureID
                featureID = feature.ID
        
        # Validate that the feature contains the relevant fields
        if sequenceType.lower() == "cds":
            assert hasattr(feature, "CDS"), \
                "CDS feature type is requested of feature '{0}' which lacks CDS".format(featureID)
        elif sequenceType.lower() == "exon":
            assert hasattr(feature, "exon"), \
                "exon feature type is requested of feature '{0}' which lacks exon".format(featureID)
        
        # Get the coordinates required for sequenceType retrieval
        if sequenceType.lower() == "exon":
            featureCoord, featureFrame = GFF3._get_feature_coords(feature, "exon")
            featureType = feature.type
            featureID = feature.ID
        elif sequenceType.lower() == "cds":
            featureCoord, featureFrame = GFF3._get_feature_coords(feature, "CDS")
            featureType = feature.type
            featureID = feature.ID
        
        # Reverse the coord lists if we're looking at a '-' model so we start at the 3' end of the gene model
        '''
        I truly have no idea why I do this. It was in my legacy code I wrote years back,
        and I have to assume I did it for a reason. I THINK it's because some GFF3s have
        truncated sequences, and when it's truncated at the 3' end of a -ve stranded gene,
        it won't accommodate that. So when we flip it, we need to take the final frame as
        our starting frame unlike the usual, non-truncated case where it'll always start
        in the 0-frame. Sorting beforehand is hence just a way of sorting the frames, not
        the coords.
        '''
        if feature.strand == '-':
            featureCoord.reverse()
        
        # Get the starting frame for the sequence
        startingFrame = featureFrame[0]
        return featureCoord, startingFrame, featureType
    
    def retrieve_sequence_from_FASTA(self, contigSequences, feature_or_featureID, sequenceType):
        '''
        Using this GFF3 instance, retrieve the exon or CDS sequence for the corresponding
        sequenceID.
        
        Parameters:
            contigSequences -- a ZS_SeqIO.FASTA, a ZS_SeqIO.FastASeq, or a Biopython
                               SeqIO.to_dict() object that contains the contig sequence
                               that sequenceID is located on/within.
            feature_or_featureID -- a string corresponding to a feature within this GFF3 object.
                                    or just a feature object in general
            sequenceType -- a string corresponding to the type of sequence to retrieve
                            i.e., in the list ["CDS", "exon"]
        Returns:
            FastASeq_objs -- a list containing ZS_SeqIO.FastASeq objects. Why? Because if
                             you provide a gene ID here, there may be multiple mRNAs that
                             are its children.
            featureTypes -- a list containing strings indicating the feature types that
                            have been returned e.g., mRNAs or lnc_RNAs.
            startingFrames -- a list containing the starting frame for any potential
                              translations to occur (if applicable) as a string.
        '''
        # Figure out if we're handling a feature or featureID
        if isinstance(feature_or_featureID, str):
            featureID = feature_or_featureID
            assert featureID in self.features, \
                "'{0}' is not recognised as a feature within this GFF3".format(featureID)
            feature = self.features[featureID]
        elif type(feature_or_featureID).__name__ == "Feature" \
            or type(feature_or_featureID).__name__ == "ZS_GFF3IO.Feature" \
            or (hasattr(feature_or_featureID, "isFeature") and feature_or_featureID.isFeature is True):
                feature = feature_or_featureID
                featureID = feature.ID
        
        # Get the coordinates required for sequenceType retrieval
        "Relevant validations are performed by retrieve_coords()"
        featureCoord, startingFrame, featureType = self.retrieve_coords(feature, sequenceType)
        
        # Get the contig sequence as a string, regardless of input type
        contigString = GFF3._get_contig_as_string(contigSequences, feature.contig)
        
        # Get VCF-dict object
        vcfDict = self._get_artifacts_as_vcfDict()
        
        # Create sequence by piecing together exon / CDS bits
        sequence = ""
        for start, end in featureCoord:
            sequenceBit = contigString[start-1:end] # 1-based correction to start to make it 0-based
            
            # Edit the sequence bit with any suggestions from vcfDict
            if feature.contig in vcfDict:
                for editCoords, editType, editResidues in vcfDict[feature.contig]: # gives us a list of lists
                    editStart, editEnd = editCoords # deconstruct the [start, end] coordinate value
                    bitStart = editStart - start # gives 0-based position
                    bitEnd = editEnd - start # also 0-based
                    if start <= editEnd and editStart <= end: # if it overlaps
                        if editType == 'deletion_artifact':
                            sequenceBit = sequenceBit[:bitStart] + sequenceBit[bitEnd+1:]
                        elif editType == 'substitution_artifact':
                            sequenceBit = sequenceBit[:bitStart] + editResidues + sequenceBit[bitEnd+1:]
                        elif editType == 'insertion_artifact':
                            sequenceBit = sequenceBit[:bitStart] + editResidues + sequenceBit[bitEnd:]
            
            # Store it
            sequence += sequenceBit
        
        # Reverse complement if necessary
        if feature.strand == "-":
            sequence = FastASeq.get_reverse_complement(self=None, staticSeq=sequence)
        
        # Create FastASeq object to represent this sequence
        FastASeq_obj = FastASeq(id=featureID, seq=sequence)
        return FastASeq_obj, featureType, startingFrame
    
    def __getitem__(self, key):
        return self.features[key]
    
    def __setitem__(self, key, item):
        self.features[key] = item
    
    def __len__(self):
        return len(self.features)
    
    def __delitem__(self, key):
        del self.features[key]
    
    def __iter__(self):
        return iter(self.features)
    
    def __contains__(self, item):
        return item in self.features
    
    def has_key(self, key):
        return key in self.features
    
    def keys(self):
        return self.features.keys()
    
    def values(self):
        return self.features.values()
    
    def items(self):
        return self.features.items()
    
    def __repr__(self):
        return "<GFF3 object;file='{0}';num_contigs={1};{2}>".format(
            self.fileLocation,
            len(self.contigs),
            ";".join(["num_{0}={1}".format(key, len(self.types[key])) for key in self.types.keys()])
        )

class LinesGFF3(GFF3):
    '''
    This subclass of GFF3 is intended to add various features that I don't want
    to clutter the base class with. Specifically, it allows for behaviours seen
    in my older GFF3 parsing class that is (as of writing this comment) still
    strewn throughout my Genome_analysis_scripts git repo.
    '''
    def __init__(self, file_location, strict_parse=True):
        super().__init__(file_location, strict_parse)
        self.isLinesGFF3 = True
    
    def add_comments(self): # This function is just add_lines but with the gene lines section gutted
        '''
        This function provides special functionalities when handling GFF3s of a format I (zkstewart)
        have created. Specifically, it allows us to hold onto any comments which can contain information
        as to the sequence. This proved important in the past (not sure as of now) when handling files
        output by PASA since it would sometimes annotate an incorrect frame which would give a different
        translation. We needed to get the specific translation that PASA output in its own comment.
        
        As said, it's unsure how useful this is today, but it's just legacy code so eh.
        '''
        KNOWN_HEAD_COMMENTS = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND', '# GEMOMA ANNOTATION', '# APOLLO ANNOTATION')
        KNOWN_FOOT_COMMENTS = ('#PROT')
        assert self.fileLocation != None
        
        with open(self.fileLocation, 'r') as file_in:
            for line in file_in:
                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                
                # Skip filler lines
                if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}: # If this is true, it's a blank line or a comment line with no information in it
                    continue
                
                # Handle known header comment lines
                if line.startswith(KNOWN_HEAD_COMMENTS):
                    # Extract gene ID
                    mrna_ID = line.split(': ')[1].split(' ')[0].rstrip(',') # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                    try:
                        gene_ID = self[mrna_ID].Parent
                    except:
                        gene_ID = self[mrna_ID].parent
                    # Add to Feature object
                    if not hasattr(self[gene_ID], "lines"):
                        self[gene_ID].add_attributes({'lines': {0: [line], 1: [], 2: []}})
                    else:
                        self[gene_ID].lines[0].append(line)
                # Handle known footer comment lines
                elif line.startswith(KNOWN_FOOT_COMMENTS):
                    # Extract gene ID
                    gene_ID = line.split()[2] # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                    
                    # Add to Feature object
                    if not hasattr(self[gene_ID], "lines"):
                        self[gene_ID].add_attributes({'lines': {0: [line], 1: [], 2: []}})
                    else:
                        self[gene_ID].lines[2].append(line)
                # Handle all other lines
                else:
                    pass
    
    def add_lines(self):
        '''
        This function provides special functionalities when handling GFF3s. Specifically, it
        allows us to write out new GFF3s using the raw lines obtained from an existing GFF3.
        This can be useful when trying to subset a GFF3 to only certain features.
        
        Honestly, there'd be much better ways to handle this. But, I want to maintain my
        legacy functions as much as possible without requiring whole-scale rewrites.
        '''
        KNOWN_HEAD_COMMENTS = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND', '# GEMOMA ANNOTATION', '# APOLLO ANNOTATION')
        KNOWN_FOOT_COMMENTS = ('#PROT')
        assert self.fileLocation != None
        
        with open(self.fileLocation, 'r') as file_in:
            for line in file_in:
                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                
                # Skip filler lines
                if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}: # If this is true, it's a blank line or a comment line with no information in it
                    continue
                sl = line.rstrip('\n').split('\t')
                
                # Handle known header comment lines
                if line.startswith(KNOWN_HEAD_COMMENTS):
                    # Extract gene ID
                    mrna_ID = line.split(': ')[1].split(' ')[0].rstrip(',') # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                    try:
                        gene_ID = self[mrna_ID].Parent
                    except:
                        gene_ID = self[mrna_ID].parent
                    # Add to Feature object
                    if not hasattr(self[gene_ID], "lines"):
                        self[gene_ID].add_attributes({'lines': {0: [line], 1: [], 2: []}})
                    else:
                        self[gene_ID].lines[0].append(line)
                
                # Handle known footer comment lines
                elif line.startswith(KNOWN_FOOT_COMMENTS):
                    # Extract gene ID
                    gene_ID = line.split()[2] # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                    
                    # Add to Feature object
                    if not hasattr(self[gene_ID], "lines"):
                        self[gene_ID].add_attributes({'lines': {0: [line], 1: [], 2: []}})
                    else:
                        self[gene_ID].lines[2].append(line)
                
                # Handle feature detail lines
                elif not line.startswith('#'):
                    # Immediately skip over any malformatted lines
                    if len(sl) < 9:
                        continue
                    
                    # Extract gene ID whilst accounting for multi-parent features
                    attributes = sl[8].split(';')
                    if sl[2] in self.parentTypes:
                        for attr in attributes:
                            if attr.startswith('ID='): # For parent-type lines, the ID= is our gene/feature ID
                                gene_ID = attr[3:].strip('\n') # This trims off the ID= bit and any new lines
                    else:
                        parentID = None
                        for attr in attributes:
                            if attr.startswith('Parent='): # For every other type of line, the Parent= field should tell us the geneID or mrnaID
                                parentID = attr[7:].strip('\n') # This trims off the Parent= bit and any new lines
                        if parentID == None: # This will handle biological_region and other values which lack ID= and Parent= fields; we don't index these since they are (currently) of no interest
                            continue
                        if "," in parentID:
                            gene_ID = parentID.split(",")
                        else:
                            gene_ID = self[parentID].ID if not hasattr(self[parentID], "Parent") else self[parentID].Parent
                    
                    # Add to lines dict
                    if type(gene_ID) != list:
                        if not hasattr(self[gene_ID], "lines"):
                            self[gene_ID].add_attributes({'lines': {0: [], 1: [line], 2: []}})
                        else:
                            self[gene_ID].lines[1].append(line)
                    else:
                        for parent in gene_ID:
                            parent_text = line.split('Parent=')[1].split(';')[0] # This will extract just the bit of the comment from Parent= to any potential ; after
                            new_line = line.replace(parent_text, parent)
                            
                            if not hasattr(self[parent], "lines"):
                                self[parent].add_attributes({'lines': {0: [], 1: [new_line], 2: []}})
                            else:
                                self[parent].lines[1].append(new_line)
                # All other lines are ignored
                else:
                    pass
    
    def pasaprots_extract(self):
        '''
        Using the comments obtained from add_comments() or [TBD()], extracts any protein sequences
        from the foot comments. Could be useful if the GFF3's frame values are incorrect with respect
        to the actual translation desired, which was a problem in the past with PASA/Augustus/some program
        in my method of gene annotation. Unsure as to its relevance today. 
        '''
        self.pasa_prots = {}
        
        for geneFeature in self.gene_values:
            # Skip any genes not indexed with lines
            if not hasattr(geneFeature, "lines"):
                continue
            
            # Parse each foot comment to extract the protein sequence
            foot_comments = geneFeature.lines[2]
            for comment in foot_comments:
                split_comment = comment.rstrip('\r\n').split('\t')
                
                # Extract the mRNA ID
                mrnaID = split_comment[0].split(' ')[1] # Format for PASA comments after ' ' split should be ['#PROT', mrnaID, geneID]
                
                # Extract the sequence
                sequence = split_comment[1]
                
                # Add into output dict
                assert mrnaID not in self.pasa_prots # If this assertion fails, GFF3 comment format is flawed - there is a duplicate mRNA ID
                self.pasa_prots[mrnaID] = sequence

if __name__ == "__main__":
    pass
