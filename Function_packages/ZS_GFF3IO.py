#! python3
# ZS_GFF3IO.py
# Contains the GFF3 class and its associated Feature
# Class. Can be used for (memory-heavy) GFF3 file parsing

import re
import pandas as pd
from collections import OrderedDict
from ncls import NCLS
from copy import deepcopy

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
    def __init__(self, file_location):
        self.fileLocation = file_location
        self.features = OrderedDict()
        self.types = {}
        self.contigs = set()
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        
        self.parse_gff3()
    
    def _make_feature_case_appropriate(self, featureType):
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
    
    def parse_gff3(self):
        # Gene object loop
        with open(self.fileLocation, 'r') as fileIn:
            for line in fileIn:
                line = line.replace('\r', '')
                
                # Skip filler and comment lines
                if line == "\n" or line.startswith("#"):
                    continue
                
                # Extract information from this line
                contig, source, featureType, start, end, \
                    score, strand, frame, attributes \
                    = line.rstrip('\n').split('\t')
                attributesDict = {a.split("=")[0]: a.split("=")[1] for a in attributes.split(";")}
                self.contigs.add(contig)
                
                # Ensure case conformity
                featureType = self._make_feature_case_appropriate(featureType)
                
                # Skip un-indexable features
                if 'ID' not in attributesDict: # see the human genome GFF3 biological_region values for why this is necessary
                    continue
                
                # Handle parent-level features
                if 'Parent' not in attributesDict: # If no Parent field this should BE the parent
                    featureID = attributesDict["ID"]
                    
                    # End parsing if duplicate ID is found
                    assert featureID not in self.features, \
                        "'{0}' feature occurs twice indicating poorly formatted file; \
                        for debugging, line == {1}; parsing will stop now".format(featureID, line)
                    
                    # Create feature and index it
                    feature = Feature()
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
                
                # Handle subfeatures
                else:
                    featureID = attributesDict["ID"]
                    parents = attributesDict["Parent"].split(',')
                    for parentID in parents:
                        # End parsing if parent doesn't exist
                        assert parentID in self.features, \
                            "'{0}' feature points to a non-existing parent '{1}' at the time of parsing; \
                            for debugging, line == {2}; parsing will stop now".format(featureID, parentID, line)
                        
                        # End parsing if duplicate ID is found
                        if featureType.lower() != "cds": # CDS features are the ONLY feature allowed to have non-unique IDs
                            assert self.features[parentID].retrieve_child(featureID) == None, \
                                "'{0}' feature is associated to a parent '{1}' more than once; \
                                for debugging, line == {1}; parsing will stop now".format(featureID, parentID, line)
                        
                        # Create feature and index it
                        feature = Feature()
                        if featureType.lower() != "cds": # since CDS features don't have unique IDs, there's no point indexing them
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
            self.mrna_values = self.types["mRNA"]
            
            # Sort contigs
            self.contigs = list(self.contigs)
            try:
                self.contigs.sort(key = lambda x: list(map(int, re.findall(r'\d+', x)))) # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
            except:
                self.contigs.sort()
    
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
                        These Features are deepcopied, so you can do whatever you want
                        to them without altering the GFF3 object.
        '''
        assert self.ncls != None and self._nclsIndex != None, \
            "Run create_ncls_index before you call this method!"
        
        overlaps = self.ncls.find_overlap(start, stop+1) # Although our ncls is already 1-based, find_overlap acts as a range. We need to +1 to keep everything logically 1-based.
        
        features = []
        for result in overlaps: # result == [start, end, index]
            features.append(self._nclsIndex[result[2]]) # this gets the Feature object from the ._nclsIndex dictionary
        features = deepcopy(features) # deepcopy so we don't alter the underlying NCLS structure
        
        # Narrow down our features to hits with a .field == value match
        features = self._ncls_feature_narrowing(features, field, value)
        
        # Return list
        return features
    
    def _ncls_feature_narrowing(self, nclsEntries, field, value):
        '''
        Hidden method for narrowing retrieved NCLS hits to the specified
        field==value match. Usually used to narrow down to contig hits,
        so field="contig" and value equal to whatever contig you want to find
        hits on.
        '''
        raiseWarning = False
        for k in range(len(nclsEntries)-1, -1, -1):
            try:
                if nclsEntries[k].__dict__[field] != value:
                    del nclsEntries[k]
            except:
                raiseWarning = True
        if raiseWarning:
            print("Warning: '{0}' wasn't found as a field in at least one of the found Features; your results might be incomplete!".format(field))
        return nclsEntries
    
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

if __name__ == "__main__":
    pass
