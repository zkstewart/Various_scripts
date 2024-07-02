#! python3
# ZS_GO.py
# Various statics for working with gene ontologies

import re, requests
from goatools import obo_parser

class GO:
    '''
    Parameters:
        obos -- a string OR list of strings pointing to the location of GO files parseable with
                obo_parser.GODag
    '''
    GO_REGEX = re.compile(r"GO:\d{7}")
    
    def __init__(self, obos):
        self.obos = obos # implicitly sets self.dags
        self.dags = self._parse_obos()
        self.queriedGOs = {}
    
    @property
    def obos(self):
        return self._obos
    
    @obos.setter
    def obos(self, obos):
        if type(obos) == str:
            obos = [obos]
        assert type(obos) == list, "obos must be a string or list of strings"
        
        self._obos = obos
    
    def _parse_obos(self):
        dags = []
        for obo in self.obos:
            dags.append(obo_parser.GODag(obo))
        return dags
    
    def convert_codes_to_names(self, goTerms):
        '''
        Receives the input goTerms object and converts all GO codes to their respective names.
        The input object can be flexibly specified as a string, a list, or a list of lists.
        The output will be returns in the same format as the input, but with GO codes replaced
        by their names. Note that values which do not pattern match as a GO code will
        be returned as-is.
        
        Parameters:
            goTerms -- a flexible object that can be a string of a single GO code,
                       a list of GO codes, or list of lists containing GO codes
        Returns
            goNames -- an output whose format mimics the input, but with terms converted
                       to names where possible
        '''
        # Handle a string (i.e., GO code) value
        if type(goTerms) == str:
            if not self.GO_REGEX.match(goTerms):
                return goTerms
            else:
                # Fix obsoletion if necessary
                self.fix_obsoletion(goTerms)
                
                # Obtain the GO objects
                gos = self[goTerms]
                
                # Handle None returns
                if gos == None:
                    "This happens when the GO term is not found in the obo file and API query fails"
                    return None
                
                # Handle list returns
                elif type(gos) == list:
                    "This happens when the GO term is obsolete and we found replacement(s)"
                    return goNames
                
                # Handle single returns
                else:
                    "This happens when the GO term is found in the obo file"
                    return gos.name
        
        # Recursively handle lists
        else:
            goNames = []
            for goTerm in goTerms:
                goNames.append(self.convert_codes_to_names(goTerm))
            
            return goNames
    
    def fix_obsoletion(self, goTerm):
        '''
        Function to check if a GO term is obsolete and, if so, attempt to fix it.
        Non-obsolete GO terms will be a no-op.
        
        Parameters:
            goTerm -- a string representing a GO term to check for obsoletion and fix if necessary
        '''
        if goTerm in self:
            return None
        elif goTerm in self.queriedGOs:
            return None
        else:
            self.queriedGOs[goTerm] = query_go_api(goTerm)
            return None
    
    def fix_obsoletions(self, goTermList):
        '''
        Assistant function to pipe a list of GO terms through fix_obsoletion
        
        Parameters:
            goTermList -- a list containing strings of GO terms
        '''
        for goTerm in goTermList:
            self.fix_obsoletion(goTerm)
    
    def __contains__(self, key):
        for dag in self.dags:
            if key in dag:
                return True
        if key in self.queriedGOs:
            return True
        
        return False
    
    def __getitem__(self, key):
        assert type(key) == str, "GO key must be a string"
        
        for dag in self.dags:
            if key in dag:
                return dag[key]
        if key in self.queriedGOs:
            return self.queriedGOs[key]
        
        # return self.fix_obsoletion(key)
        raise KeyError(f"GO term '{key}' not found in any of the parsed obo files")
    
    def __repr__(self):
        return "<GO object;obos='{0}';num_terms={1}>".format(
            ", ".join(self.obos),
            ", ".join(str(len(dag)) for dag in self.dags)
        )

def fix_obsoletions(goList, goObo, queriedGOs):
    '''
    Helper function to receive a list of GO terms, alongside the parsed go.obo
    and a dictionary potentially containing previous API queries. Using this, it
    will check if the GOs in the GO list are obsolete and, if so, attempt to fix
    them via API query.
    
    Parameters:
        goList -- a list of string GO terms e.g., ['GO:0033644', 'GO:0016020'].
        goObo -- a obo_parser.GODag object which has parsed a go.obo file.
        queriedGOs -- a dictionary with previous API query hits (if any) that
                      will be modified by this function; has structure like:
                      {
                          'goID1': ['replacementID1'],
                          'goID2': ['replacementID1', 'replacementID2'],
                          'goID3': None, # if no replacement is available
                          ...
                      }
    '''
    fixedGOs = []
    for goTerm in goList:
        
        # If we have a hit, everything's good!
        if goTerm in goObo:
            fixedGOs.append(goTerm)
        
        # Otherwise, if we found an obsoletion, handle it
        else:
            # Perform an API query if we haven't seen this GO yet
            if not goTerm in queriedGOs:
                replacements = query_go_api(goTerm)
                queriedGOs[goTerm] = replacements
            
            # Get the API query result
            apiGOs = queriedGOs[goTerm]
            
            # Store result if we were successful
            if apiGOs != None:
                fixedGOs += apiGOs
    return fixedGOs

def query_go_api(goTerm):
    '''
    Parameters:
        goTerm -- a GO term to query for replacement to its obsoletion
    Returns:
        replacementList -- a list containing GO terms which replace this one, OR
                           None if no replacements were found
    '''
    # Query the API
    apiUrl = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A{0}/history".format(
        goTerm.split(":")[1]
    )
    response = requests.get(apiUrl)
    json = response.json()
    
    # Handle errors
    if "messages" in json:
        message = json["messages"][0]
        print(f"# WARNING: Querying the GO term '{goTerm}' resulted in an API error")
        print(f"# > The return message was: '{message}'")
        print("# > No result will be available for this GO term")
        return None
    
    # Handle null hits
    if json["results"] == []:
        print(f"# WARNING: Querying the GO term '{goTerm}' returned no results")
        print("# > No result will be available for this GO term")
        return None
    
    # Parse results out of the JSON response
    resultsDict = json["results"][0]
    
    # Handle situations where the GO isn't actually obsolete
    if resultsDict["isObsolete"] != True:
        print(f"# WARNING: The GO term '{goTerm}' does not appear to be obsolete")
        print("# > However, it is not found in the go.obo file?")
        print("# > Ancestor terms will be unavailable for this GO term")
        return [goTerm]
    
    # Process results and try to find replacements or considers
    else:
        replacedBys = set()
        considers = set()
        for historyDict in resultsDict["history"]:
            historyText = historyDict["text"]
            if "consider GO:" in historyText:
                considerMatch = re.search(r"consider (GO:\d{7})", historyText)
                if considerMatch != None:
                    considers.add(considerMatch.group(1))
            if "replaced_by GO:" in historyText:
                replacedByMatch = re.search(r"replaced_by (GO:\d{7})", historyText)
                if replacedByMatch != None:
                    replacedBys.add(replacedByMatch.group(1))
    
    # Return the best value, or None if no replacements were possible
    if len(replacedBys) > 0:
        return list(replacedBys)
    elif len(considers) > 0:
        return list(considers)
    else:
        return None

if __name__ == "__main__":
    pass
