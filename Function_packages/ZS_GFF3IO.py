#! python3
# ZS_GFF3IO.py
# Contains the GFF3 class and its associated Feature
# Class. Can be used for (memory-heavy) GFF3 file parsing

import re, sys, os
import pandas as pd
from pyfaidx import Fasta
from collections import OrderedDict
from ncls import NCLS
from hashlib import md5

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
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
    
    def del_child(self, childFeatureID, childFeatureType):
        '''
        This method is intended to delete child Features from this Feature.
        In the process of deleting the child feature, it will clean up
        fields as necessary if they become empty.
        
        It will also make sure this feature's coordinate values make sense
        after deletion of the child by calling .update_coordinates().
        
        Parameters:
            childFeatureID -- the child features .ID attribute value
        '''
        # From .children
        parentChildIndices = [i for i in range(len(self.children)) if self.children[i].ID == childFeatureID]
        for index in reversed(parentChildIndices):
            del self.children[index]
        
        # From .types
        parentTypeIndices = [i for i in range(len(self.types[childFeatureType])) if self.types[childFeatureType][i].ID == childFeatureID]
        for index in reversed(parentTypeIndices):
            del self.types[childFeatureType][index]
        if len(self.types[childFeatureType]) == 0: # if this .type[key] is now empty
            del self.types[childFeatureType]
        
        # From attributes associated with .types
        parentAttributeIndices = [i for i in range(len(self.__dict__[childFeatureType])) if self.__dict__[childFeatureType][i].ID == childFeatureID]
        for index in reversed(parentAttributeIndices):
            del self.__dict__[childFeatureType][index]
        if len(self.__dict__[childFeatureType]) == 0: # if this .key is now empty
            delattr(self, childFeatureType)
        
        self.update_coordinates()
    
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
    
    def format_as_gff3(self):
        '''
        This method will attempt to render a GFF3-correct format of the
        data this object contains. Several assumptions are made which, if you
        haven't done anything truly weird, will hold true.
        '''
        # Separate object attributes according to 1) expected GFF3 columns, and 2) anything else
        gff3Attributes = ["contig", "source", "type", "start", "end", "score", "strand", "frame"]
        for attribute in gff3Attributes:
            assert attribute in self.__dict__, \
                f"Feature lacks expected GFF3 value '{attribute}'; we can't format it appropriately"
        
        objectAttributes = ["children", "types", "isFeature", "coords"] # these come from the Feature class itself
        detailAttributes = [attribute for attribute in self.__dict__ \
            if attribute not in objectAttributes and attribute not in gff3Attributes and attribute not in self.types]
        detailAttributes.sort(key = lambda x: (0 if x.lower() == "id" else 2, 1 if x.lower() == "parent" else 2)) # always bring ID and Parent to front
        
        # Format a string and return
        gff3Line = "\t".join([str(self.__dict__[attribute]) for attribute in gff3Attributes])
        gff3Line += "\t{0}".format(
            ";".join(
                ["{0}={1}".format(attribute, str(self.__dict__[attribute])) for attribute in detailAttributes]
            )
        )
        
        return gff3Line
    
    def update_coordinates(self):
        '''
        This method assumes that this Feature object is being used in a typical fashion
        for GFF3 feature storage. This means it will have .coord, .start, and .end attributes.
        
        If the children of this object have been changed, the coordinates of this parent might
        also need to change. For example, a gene feature's start and end values are contigent upon
        the earliest start and latest end of all its children values.
        
        Hence, this methods's goal is to make sure this feature's coordinates properly reflect
        the children that it contains, if any.
        '''
        earliestStart = None
        latestEnd = None
        
        for childFeature in self.children:
            if earliestStart == None or childFeature.start < earliestStart:
                earliestStart = childFeature.start
            if latestEnd == None or childFeature.end > latestEnd:
                latestEnd = childFeature.end
        
        if earliestStart != None:
            self.start = earliestStart
            self.coords[0] = earliestStart
        if latestEnd != None:
            self.end = latestEnd
            self.coords[1] = latestEnd
    
    def update_id(self, newID, GFF3_obj=None):
        '''
        This method assumes that this Feature object is being used in a typical fashion
        for GFF3 feature storage. This means it will have a .ID value, and its immediate
        children will have .Parent values of the same value.
        
        This method's goal is to update a Feature's .ID value such that its children
        reflect this change in their .Parent values (if applicable).
        
        Parameters:
            newID -- a string value to become the new .ID value for this Feature
            GFF3_obj -- OPTIONAL; if you provide this, the GFF3 class will handle
                        the change in ID for this Feature and all its children
                        with appropriate indexing changes managed as well.
        '''
        if not isinstance(newID, str):
            raise ValueError("newID should be a string")
        
        # Let GFF3 object handle the change
        if GFF3_obj != None:
            if self in GFF3_obj:
                GFF3_obj.update_id(self, newID)
            else:
                raise ValueError(f"Feature '{self.ID}' not found in GFF3 object")
        # Apply change locally to this Feature
        else:
            self.ID = newID
            for childFeature in self.children:
                childFeature.Parent = self.ID
    
    def reformat_id(self, newID, GFF3_obj=None):
        '''
        Unlike update_id which merely changes the current Feature's ID value and its
        immediate children's Parent values, this method will change the ID of this
        Feature, and also change the ID of all its children Features recursively 
        in a way that reflects the new ID.
        
        Parameters:
            newID -- a string value to become the new .ID value for this Feature, 
                     the .Parent ID for children, and provide the format for new
                     child .ID values.
            GFF3_obj -- OPTIONAL; if you provide this, the GFF3 class will handle
                        the change in ID for this Feature and all its children
                        with appropriate indexing changes managed as well.
        '''
        if not isinstance(newID, str):
            raise ValueError("newID should be a string")
        
        # Let GFF3 object handle the change
        if GFF3_obj != None:
            if self.ID in GFF3_obj:
                GFF3_obj.reformat_id(self, newID)
            else:
                raise ValueError(f"Feature '{self.ID}' not found in GFF3 object")
        # Apply change locally to this Feature
        else:
            self.ID = newID
            
            # Update the immediate child features
            typesCount = {}
            for childFeature in self.children:
                typesCount.setdefault(childFeature.type, 1)
                childFeature.Parent = self.ID
                
                # Format the new ID
                newID = f"{newID}.{childFeature.type}{typesCount[childFeature.type]}"
                typesCount[childFeature.type] += 1
                
                # Recursively update children
                childFeature.reformat_id(newID)
    
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
    def __init__(self, file_location, strict_parse=True, fix_duplicated_ids=False, slim_index=False):
        self.fileLocation = file_location
        self.features = OrderedDict()
        self.types = {}
        self.contigs = set()
        self.parentTypes = set() # tells us which feature types to expect as being parents
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        
        self.isGFF3 = True
        self.parse_gff3(strict_parse=strict_parse, fix_duplicated_ids=fix_duplicated_ids, slim_index=slim_index)
    
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
    
    @staticmethod
    def _recursively_write_feature_details(feature, fileHandle):
        '''
        Hidden helper to use recursion for writing GFF3 features to file
        '''
        # Write parent details
        fileHandle.write(feature.format_as_gff3() + "\n")
        
        # Iteratively write children details
        for childFeature in sorted(feature.children, key = lambda x: x.start):
            GFF3._recursively_write_feature_details(childFeature, fileHandle)
    
    @staticmethod
    def _get_contig_as_FastASeq(contigSequences, contigID):
        '''
        Hidden helper method of retrieve_sequence_from_FASTA(). Intended to work with
        many different object types i.e., a FASTA, FastASeq, or Biopython sequence.
        '''
        sequence = GFF3._get_contig_as_string(contigSequences, contigID)
        
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
            if contigID not in contigSequences:
                raise KeyError("'{0}' does not exist within the dictionary contig sequences object".format(contigID))
            if type(contigSequences[contigID]).__name__ != "SeqRecord":
                raise TypeError("contigSequences is a dictionary, but values aren't Biopython SeqRecords? Can't handle.")
            sequence = str(contigSequences[contigID].seq)
        
        # Handle ZS_SeqIO.FASTA type
        elif hasattr(contigSequences, "isFASTA") and contigSequences.isFASTA is True:
            try:
                sequence = contigSequences[contigID].seq
            except:
                raise KeyError("'{0}' does not exist within the FASTA object".format(contigID))
        
        # Handle ZS_SeqIO.FastASeq type
        elif hasattr(contigSequences, "isFastASeq") and contigSequences.isFastASeq is True:
            if contigID != contigSequences.id and contigID != contigSequences.alt:
                raise ValueError("'{0}' contig does not match the provided FastASeq object".format(contigID))
            sequence = contigSequences.seq
        
        # Handle pyfaidx.Fasta type
        elif type(contigSequences) == Fasta:
            if contigID not in contigSequences:
                raise KeyError("'{0}' does not exist within the pyfaidx.Fasta sequences object".format(contigID))
            sequence = str(contigSequences[contigID])
        
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
    
    def parse_gff3(self, strict_parse=True, fix_duplicated_ids=False, full_warning=False, slim_index=False):
        '''
        Parameters:
            strict_parse -- a boolean indicating whether this function should
                            die if the GFF3 doesn't meet strict format standards,
                            or if failing annotation values should simply be skipped
            full_warning -- a boolean indicating whether every single warning should
                            be printed, or just the first 10.
            slim_index -- a boolean indicating whether the GFF3 indexed should be fully
                          featured (slim_index == False) or if a slim index should be created
                          with minimal features detailed (slim_index == True)
        '''
        # Set up warning handling system
        warningContainer = { # this dict acts like a JSON for data storage
            "warningCount": 0,
            "warningLimit": 10,
            "hasHandledWarnings": False
        }
        def _handle_warning_message(warningContainer, message):
            if full_warning == True or warningContainer["warningCount"] < warningContainer["warningLimit"]:
                print(message)
                warningContainer["warningCount"] += 1
            if warningContainer["hasHandledWarnings"] == False and warningContainer["warningCount"] == warningContainer["warningLimit"]:
                print("Further warning messages will be suppressed.")
                warningContainer["hasHandledWarnings"] = True
        
        # Setup for slim parsing functionality
        def _format_attributes(attributes, slim_index):
            SLIM_ATTRIBUTES = ["id", "parent"]
            
            splitAttributes = []
            for a in attributes.split("="):
                if ";" in a:
                    splitAttributes += a.rsplit(";", maxsplit=1)
                else:
                    splitAttributes.append(a)
            
            if not slim_index:
                attributesDict = {splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)}
            else:
                attributesDict = {
                    splitAttributes[i]: splitAttributes[i+1]
                    for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)
                    if splitAttributes[i].lower() in SLIM_ATTRIBUTES
                }
            return attributesDict
        
        # Gene object loop
        lineCount = 0
        with open(self.fileLocation, 'r') as fileIn:
            for line in fileIn:
                lineCount += 1
                line = line.replace('\r', '').replace('"', '') # remove carriage returns and quotes
                
                # Skip filler and comment lines
                if line == "\n" or line.startswith("#"):
                    continue
                
                # Extract information from this line
                try:
                    contig, source, featureType, start, end, \
                        score, strand, frame, attributes \
                        = line.rstrip('\t\n').split('\t')
                except:
                    if strict_parse == True:
                        raise ValueError(
                            f"Error: Line #{lineCount} (\"{line}\") does not meet GFF3 standards; parsing failed"
                        )
                    else:
                        _handle_warning_message(warningContainer,
                            f"Warning: Line #{lineCount} (\"{line}\") does not meet GFF3 standards;" +
                            " strict parsing is disabled, so we'll just continue and hope for the best"
                        )
                        continue
                
                # Format attributes dictionary
                attributesDict = _format_attributes(attributes, slim_index)
                
                # Ensure case conformity
                featureType = GFF3.make_feature_case_appropriate(featureType)
                
                # Fix GFF3s which did not give exons an ID
                "I'm looking at you banana genome hub. You shouldn't do this."
                if 'ID' not in attributesDict and featureType.lower() == "exon":
                    parentID = attributesDict["Parent"].split(',')
                    assert len(parentID) == 1, \
                        ("I tried to fix missing exon IDs but found a sequence with >1 parent ", +
                         "i.e., {0}".format(attributesDict["Parent"]))
                    parentID = parentID[0]
                    
                    parentFeature = self.features[parentID]
                    try:
                        numExons = len(parentFeature.exon)
                    except:
                        numExons = 0
                    attributesDict["ID"] = f"{parentID}.exon{numExons+1}"
                
                # Skip un-indexable features
                if 'ID' not in attributesDict and featureType.lower() != "cds": # see the human genome GFF3 biological_region values for why this is necessary
                    continue
                self.contigs.add(contig) # we can index this contig now that we've skipped un-indexable features
                
                # Handle parent-level features
                if 'Parent' not in attributesDict: # If no Parent field this should BE the parent
                    featureID = attributesDict["ID"]
                    
                    # End parsing if duplicate ID is found (strict_parse is True)
                    if strict_parse == True:
                        assert featureID not in self.features, \
                            (f"Error: '{featureID}' feature occurs twice indicating poorly formatted file;" + 
                            f" for debugging, line #{lineCount} == {line}; parsing will stop now")
                    # Skip parsing if duplicate ID is found (strict_parse is False)
                    else:
                        if featureID in self.features:
                            _handle_warning_message(warningContainer,
                                f"Warning: '{featureID}' feature occurs more than once indicating poorly formatted file;"+
                                " strict parsing is disabled, so we will just continue and hope for the best"
                            )
                            continue
                    
                    # Create feature and populate it with details
                    feature = Feature()
                    feature.add_attributes(attributesDict)
                    if not slim_index:
                        feature.add_attributes({
                            "contig": contig, "source": source, "type": featureType,
                            "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                            "score": score, "strand": strand, "frame": frame
                        })
                    else:
                        feature.add_attributes({
                            "contig": contig, "type": featureType,
                            "start": int(start), "end": int(end),
                            "strand": strand
                        })
                    
                    # Index feature
                    self.features[featureID] = feature
                    self.types.setdefault(featureType, [])
                    self.types[featureType].append(feature)
                    self.parentTypes.add(featureType)
                
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
                                    f"Error: '{featureType}' feature lacks an ID attribute and hence cannot be indexed;" +
                                    f" for debugging, line #{lineCount} == {line}; parsing will stop now"
                                )
                            else:
                                featureID = f"{parentID}.cds"
                        
                        # End parsing if parent doesn't exist (strict_parse is True)
                        if strict_parse == True:
                            assert parentID in self.features, \
                                (f"Error: '{featureID}' feature points to a non-existing parent '{parentID}' at the time of parsing;" + 
                                f" for debugging, line #{lineCount} == {line}; parsing will stop now")
                        else:
                            if parentID not in self.features:
                                _handle_warning_message(warningContainer,
                                    f"Warning: '{featureID}' feature points to a non-existing parent '{parentID}' at the time of parsing;" + 
                                    " strict parsing is disabled, so we will just continue and hope for the best"
                                )
                                continue
                        
                        # End parsing if duplicate ID is found and we aren't fixing it
                        if featureType.lower() != "cds": # CDS features are the ONLY feature allowed to have non-unique IDs
                            if self.features[parentID].retrieve_child(featureID) != None:
                                if fix_duplicated_ids:
                                    for i in range(2, len(self.features[parentID].retrieve_all_children())+3):
                                        if self.features[parentID].retrieve_child(f"{featureID}.{i}") == None:
                                            featureID = f"{featureID}.{i}"
                                            attributesDict["ID"] = featureID
                                            break
                                else:
                                    raise ValueError(
                                        (f"Error: '{featureID}' feature is associated to a parent '{parentID}' more than once;" + 
                                        f" for debugging, line #{lineCount} == {line}; parsing will stop now")
                                    )
                        
                        # Create feature and populate it with details
                        feature = Feature()
                        feature.add_attributes(attributesDict)
                        if not slim_index:
                            feature.add_attributes({
                                "contig": contig, "source": source, "type": featureType,
                                "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                                "score": score, "strand": strand, "frame": frame
                            })
                        else:
                            feature.add_attributes({
                                "contig": contig, "type": featureType,
                                "start": int(start), "end": int(end),
                                "strand": strand
                            })
                        
                        # Index feature
                        if featureType.lower() != "cds": # since CDS features aren't guaranteed to have unique IDs, there's no point indexing them
                            self.features[featureID] = feature
                        self.types.setdefault(featureType, [])
                        self.types[featureType].append(feature)
                        
                        # Double-index for special child attributes
                        self._index_products(feature)
                        
                        # Add it as a child of the parent
                        "It's important to fully-specify the child before running add_child()"
                        self.features[parentID].add_child(feature)
            
            # Generate shortcut fields
            try:
                self.gene_values = self.types["gene"]
            except:
                print("WARNING: No gene features found in GFF3 file '{0}'".format(self.fileLocation))
                self.gene_values = None
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
                           to index OR an iterable of strings indicating multiple
                            types to index.
        '''
        if isinstance(typeToIndex, str):
            typeToIndex = [typeToIndex]
        
        for indexType in typeToIndex:
            assert indexType in self.types, \
                "'{0}' not found as a feature type within the parsed GFF3 ('{1}')".format(indexType, self.fileLocation)
        
        nclsIndex = {}
        starts, ends, ids = [], [], []
        
        # Add features of the specified type to our NCLS structure
        ongoingCount = 0
        for indexType in typeToIndex:
            for feature in self.types[indexType]:
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
    
    def retrieve_sequence_from_FASTA(self, contigSequences, feature_or_featureID, sequenceType,
                                     skipComplement=False):
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
            skipComplement -- OPTIONAL; a boolean indicating whether this function should
                              automatically reverse complement -ve stranded features.
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
        if feature.strand == "-" and skipComplement is False:
            sequence = FastASeq.get_reverse_complement(self=None, staticSeq=sequence)
        
        # Create FastASeq object to represent this sequence
        FastASeq_obj = FastASeq(id=featureID, seq=sequence)
        return FastASeq_obj, featureType, startingFrame
    
    def write(self, outputFileName, force=False, parentTypesToSkip=[]):
        '''
        Writes the GFF3 object out to file at the location provided. Formatting will be done
        in a strictly correct GFF3 style. However, it will be done in unconventional manner
        wherein parent types will be ordered first, before contig order and positional order.
        
        Params:
            outputFileName -- a string indicating the file location to write to. This file
                              must not exist or an error will be raised (unless Force==True)
            force -- a boolean indicating whether any existing file should be overwritten (True)
                     or if an error should be raised instead (False)
            parentTypesToSkip -- a list containing strings which correspond to parent types
                                 which should be omitted from the output file.
        '''
        # Validate parameters
        if not isinstance(outputFileName, str):
            raise TypeError("GFF3 can only be written to a string file path; provided outputFileName type unrecognised")
        if os.path.isfile(outputFileName) and force is False:
            raise FileExistsError(f"GFF3 can't be written to '{outputFileName}' since it already exists and force is False")
        for parentType in parentTypesToSkip:
            assert parentType in self.parentTypes, \
                f"'{parentType}' is not a recognised parent type in this GFF3 object"
        
        # Write to file
        with open(outputFileName, "w") as fileOut:
            for parentType in self.parentTypes:
                if parentType in parentTypesToSkip:
                    continue
                else:
                    for feature in sorted(self.types[parentType], key = lambda x: (x.contig, x.start)):
                        GFF3._recursively_write_feature_details(feature, fileOut)
    
    def add_feature(self, feature):
        '''
        This method will receive a fully-formed Feature, which includes children,
        and will integrate it into this GFF3 object.
        
        Parameters:
            feature -- a Feature object which is fully specified with/without children.
        '''
        
        if feature.type.upper() != "CDS":
            # Add the Feature type if it's new
            if not hasattr(feature, "Parent"):
                self.parentTypes.add(feature.type)
            
            # Integrate this Feature
            self.features[feature.ID] = feature
            self.types.setdefault(feature.type, [])
            self.types[feature.type].append(feature)
            
            # Integrate children
            for childFeature in feature.children:
                self.add_feature(childFeature)
        else:
            self.types.setdefault("CDS", [])
            self.types["CDS"].append(feature)
    
    def update_id(self, feature, newID):
        '''
        This method's goal is to update a Feature's .ID value alongside its children's
        .Parent values.
        
        Parameters:
            feature -- a Feature part of this GFF3 object
            newID -- a string value to be used for replacing the original ID
                     of the featue and any child .ID and .Parent values.
        '''
        if not isinstance(newID, str):
            raise ValueError("newID should be a string")
        if not feature in self:
            raise ValueError(f"Feature with ID '{feature.ID}' not found in this GFF3 object")
        
        # Re-index the feature
        self.reindex_feature(feature, newID)
        
        # Update child .Parent values
        for childFeature in feature.children:
            childFeature.Parent = feature.ID
    
    def reformat_id(self, feature, newID):
        '''
        This method's goal is to update a Feature's .ID value such that its children's
        .Parent values reflect this, as well as reformatting the child's ID value to
        reflect the new parent ID.
        
        Parameters:
            feature -- a Feature part of this GFF3 object
            newID -- a string value to be used for replacing the original ID
                     of the featue and any child .ID and .Parent values.
        '''
        if not isinstance(newID, str):
            raise ValueError("newID should be a string")
        if (feature.type != "CDS") and (not feature in self):
            raise ValueError(f"Feature with ID '{feature.ID}' not found in this GFF3 object")
        
        # Re-index the feature
        self.reindex_feature(feature, newID)
        
        # Update the immediate child features
        typesCount = {}
        for childFeature in feature.children:
            typesCount.setdefault(childFeature.type, 1)
            childFeature.Parent = feature.ID
            
            # Format the new child ID
            childID = f"{newID}.{childFeature.type}{typesCount[childFeature.type]}"
            typesCount[childFeature.type] += 1
            
            # Recursively update children
            self.reformat_id(childFeature, childID)
    
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
    
    def reindex_feature(self, feature, newID):
        '''
        Function to reindex a feature after a change in ID. It will specifically
        update the .features dictionary and the .types dictionary.
        
        Parameters:
            feature -- a Feature object to have its ID changes
            newID -- a string indicating what the new feature.ID value should be
        '''
        # Clear .features of the feature indexed by ID
        if feature.type != "CDS": # CDS features are not indexed in the .features dictionary
            self.features.pop(feature.ID)
        
        # Index the feature with its new ID
        feature.ID = newID
        if feature.type != "CDS":
            self.features[feature.ID] = feature
    
    def _index_products(self, feature):
        '''
        Hidden helper for indexing features with the .product attribute. These are found in
        some GFF3 files, and often any protein files generated from these will have the
        product ID associated to the translated features. This can prove problematic when
        trying to find their coordinates in the GFF3 since we'll find no matches.
        
        This method will index the product so it's discoverable within the GFF3 object.
        It won't be listed as a parentType (it should always be under a gene parent), and
        
        '''
        if hasattr(feature, "Product"):
            "We want the case of this to be predictable"
            feature.product = feature.Product
            delattr(feature, "Product")
        if hasattr(feature, "product"):
            self.features[feature.product] = feature
            self.types.setdefault("product", [])
            self.types["product"].append(feature)
    
    def eliminate_feature(self, thisFeature, isInRecursion=False):
        '''
        When we're eliminating a feature, there's a few places we need
        to look for references. At a whole class level, the GFF3 object will
        index anything other than CDS in the .feature dictionary. That's an
        easy target. It will also store everything in the .types dictionary,
        which itself indexes the various feature types alongside a list of
        every feature of that type. It's less efficient to do this since
        they do not have keys, but it's still doable.
        
        Parameters:
            thisFeature -- the feature to remove from this GFF3 object
            isInRecursion -- a boolean that you shouldn't touch. It will
                             let this method know if it's being called
                             from within a recursion or not. If it is
                             in a recursion, it won't try to touch Parents
                             lest things break. It is only the initial
                             method call that will escalate elimination
                             up to the Parent if applicable
        '''
        # Recurse into this feature's children first (depth-first deletion)
        for childFeature in thisFeature.children:
            self.eliminate_feature(childFeature, isInRecursion=True)
        
        # Eliminate from .features
        if thisFeature.ID in self.features:
            del self.features[thisFeature.ID]
        
        # Eliminate from .types
        indexToRemove = [
            index
            for index, f in enumerate(self.types[thisFeature.type])
            if f.ID == thisFeature.ID
        ]
        assert len(indexToRemove) == 1, \
            f"eliminate_feature ERROR: Found {len(indexToRemove)} instances of '{thisFeature.ID}' in '{thisFeature.type}' type"
        del self.types[thisFeature.type][indexToRemove[0]]
        
        # Eliminate the type if it's now empty
        if len(self.types[thisFeature.type]) == 0:
            del self.types[thisFeature.type]
        
        # Eliminate from the parent
        if isInRecursion is False and hasattr(thisFeature, "Parent"):
            parentFeature = self[thisFeature.Parent]
            parentFeature.del_child(thisFeature.ID, thisFeature.type)
            
            # Eliminate parent feature if it's now child free
            if len(parentFeature.children) == 0:
                self.eliminate_feature(self[thisFeature.Parent])
    
    def __getitem__(self, key):
        return self.features[key]
    
    def __setitem__(self, key, item):
        self.features[key] = item
    
    def __len__(self):
        return len(self.features)
    
    def __delitem__(self, key):
        if key not in self.features:
            raise ValueError(f"'{key}' not found in this GFF3 object")
        else:
            featureToBeDeleted = self.features[key]
            self.eliminate_feature(featureToBeDeleted)
    
    def __iter__(self):
        return iter(self.features)
    
    def __contains__(self, item):
        return item.ID in self.features
    
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
    
    def __hash__(self):
        '''
        This hashing function assumes the file name for the GFF3 is identical,
        which should be a reasonable assumption in most cases. It otherwise doesn't
        care where the GFF3 is located, so long as its contents remain identical
        which is checked by referring to its .types values.
        '''
        hashString = "<GFF3 object;file='{0}';num_contigs={1};{2}>".format(
            os.path.basename(self.fileLocation),
            len(self.contigs),
            ";".join(["num_{0}={1}".format(key, len(self.types[key])) for key in self.types.keys()])
        )
        hash = int(md5(hashString.encode("utf-8")).hexdigest(), 16)
        return hash

class LinesGFF3(GFF3):
    '''
    This subclass of GFF3 is intended to add various features that I don't want
    to clutter the base class with. Specifically, it allows for behaviours seen
    in my older GFF3 parsing class that is (as of writing this comment) still
    strewn throughout my Genome_analysis_scripts git repo.
    '''
    def __init__(self, file_location, strict_parse=True, fix_duplicated_ids=False, slim_index=False):
        super().__init__(file_location,
                         strict_parse=strict_parse,
                         fix_duplicated_ids=fix_duplicated_ids,
                         slim_index=slim_index)
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
