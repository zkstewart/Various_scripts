#! python3
# gff3_id_mapping.py
# Uses the ZS_GFF3IO.GFF3 Class to parse a GFF3 file and
# allow flexible logic for mapping features to one or more
# attribute values.

import os, argparse, re
from collections import OrderedDict

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
    def __init__(self, file_location, strict_parse=True):
        self.fileLocation = file_location
        self.features = OrderedDict()
        self.types = {}
        self.contigs = set()
        
        self.parse_gff3(strictParse=strict_parse)
    
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
    
    def parse_gff3(self, strictParse=True):
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
                if 'ID' not in attributesDict and featureType.lower() != "cds": # see the human genome GFF3 biological_region values for why this is necessary
                    continue
                
                # Handle parent-level features
                if 'Parent' not in attributesDict: # If no Parent field this should BE the parent
                    featureID = attributesDict["ID"]
                    
                    # End parsing if duplicate ID is found (strictParse is True)
                    if strictParse == True:
                        assert featureID not in self.features, \
                            "'{0}' feature occurs twice indicating poorly formatted file; \
                            for debugging, line == {1}; parsing will stop now".format(featureID, line)
                    # Skip parsing if duplicate ID is found (strictParse if False)
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
                                    for debugging, line == {2}; parsing will stop now".format(featureType, line)
                                )
                            else:
                                featureID = "{0}.cds".format(parentID)
                        
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
            self.mrna_values = self.types["mRNA"]
            
            # Sort contigs
            self.contigs = list(self.contigs)
            try:
                self.contigs.sort(key = lambda x: list(map(int, re.findall(r'\d+', x)))) # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
            except:
                self.contigs.sort()
    
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

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.gff3File):
        print('I am unable to locate the directory containing .stats files (' + args.statsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def gff3_id_mapper(forEach, map, to, gff3Obj):
    '''
    Provides a way to form a structured query of the GFF3 object to
    retrieve details in a list.
    
    Parameters:
        forEach -- a string indicating the type of GFF3 feature to check.
                   Whatever value is given here, must be a key within
                   the gff3Obj.types dictionary.
        map -- a string indicating the "primary key" that will be
               the left-most column of the mapping list.
        to -- a list containing strings of any valid keys that can be
              obtained from the GFF3 feature that is being iterated over
              with the forEach parameter.
        gff3Obj -- a ZS_GFF3IO.GFF3 instance.
    Returns:
        mapping -- a list containing sublists of mappings.
    '''
    # Conform and validate inputs
    forEach = gff3Obj._make_feature_case_appropriate(forEach)
    assert forEach in gff3Obj.types, \
        "'{0}' forEach value does not exist within the GFF3; \
        valid options include '{1}'".format(forEach, list(gff3Obj.types.keys()))
    
    # Perform the query operation
    mapping = []
    for feature in gff3Obj.types[forEach]:
        assert map in feature.__dict__, \
            "'{0}' map key does not exist within one or more features; \
            valid keys include '{1}'; check that letters have correct \
            upper/lower casing".format(map, list(feature.__dict__.keys()))
        mappingDetail = feature.__dict__[map]
        
        toDetails = []
        for toKey in to:
            assert toKey in feature.__dict__, \
                "'{0}' to key does not exist within one or more features; \
                valid keys include '{1}'; check that letters have correct \
                upper/lower casing".format(map, list(feature.__dict__.keys()))
            toDetails.append(feature.__dict__[toKey])
        
        mapping.append([mappingDetail] + toDetails)
    return mapping

if __name__ == "__main__":
    # User input
    usage = """%(prog)s accepts a GFF3 file and provides a simple way to format
    a statement that will extract relevant information mappings in a TSV output file.
    
    The statement consists of 3 parameters, namely being -forEach __ -map __ -to __.
    Logically, this should form a kind of sentence where we are saying e.g., "for each mRNA,
    map ID to Parent". In this example, the output file would give columns linking the mRNA 
    ID to its parent gene ID. Multiple -to __ keys are accepted, and each additional key will
    give an extra column in the output file. Separate -to __ keys with a space.
    
    Valid keys are any attribute key=value pairing in column 9 of the GFF3. Alternatively,
    the first 8 columns have their own key names which are: contig, source, ___ (specify in
    -forEach), start, end, score, strand, frame. Keys are case-sensitive!
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    ## Logical structure parameters
    p.add_argument("-forEach", dest="forEach", required=True,
        help="""
            The value here should relate to any feature type in column 3 of the GFF3;
            examples include gene or mRNA.
        """)
    p.add_argument("-map", dest="map", required=True,
        help="""
            The value here should relate to any valid key that will be the left-most column
            of the output and serve as the anchor for other details.
        """)
    p.add_argument("-to", dest="to", required=True, nargs="+",
        help="One or more valid keys separated with a space.")
    ## File I/O parameters
    p.add_argument("-g", dest="gff3File", required=False,
                help="Specify the location of the GFF3 file")
    p.add_argument("-o", dest="outputFileName", required=False,
                help="Specify the name and location for the output file")
    
    args = p.parse_args()
    validate_args(args)
    
    # Load GFF3
    gff3Obj = GFF3(args.gff3File, False) # non-strict parsing
    
    # Obtain data linkings
    mapping = gff3_id_mapper(args.forEach, args.map, args.to, gff3Obj)
    
    # Write output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("#{0}\t{1}\n".format(args.map, "\t".join(args.to))) # header
        for row in mapping:
            fileOut.write("{0}\n".format("\t".join(row)))
