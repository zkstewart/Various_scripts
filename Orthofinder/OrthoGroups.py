#! python3

import os
from Bio import SeqIO

class OrthoGroups:
    '''
    Attributes:
        groups -- a dictionary with structure like:
                  {
                      orthogroupID: {
                          species1: [seq1, seq2, ...],
                          species2: [seq3, seq4, ...],
                          ...
                      }
                  }
        sequences -- a dictionary with structure like:
                     {
                         species1: SeqIO.to_dict(file1),
                         species2: SeqIO.to_dict(file2),
                         ...
                     }
        species -- a dictionary with structure like:
                   {
                       species1: {
                           orthogroup1: [seq1, seq2, ...],
                           ...
                       }
                       species2: {
                           orthogroup2: [seq3, seq4, ...],
                           ...
                       },
                       ...
                   }
    '''
    def __init__(self, orthogroupsTSV, fastaDict):
        '''
        Parameters:
            orthoGroupsTSV -- the path to the Orthogroups.tsv file
            fastaDict -- a dictionary with structure like:
                         {
                             species1: path1,
                             species2: path2,
                             ...
                         }
        '''
        self.tsv = orthogroupsTSV
        self.fastaDict = fastaDict
        
        self.groups, self.species = {}, {}
        self._parse_groups()
        
        self.sequences = self._load_sequences()
    
    @property
    def tsv(self):
        return self._tsv
    
    @tsv.setter
    def tsv(self, value):
        assert os.path.isfile(value), f"OrthoGroups file '{value}'  file not found"
        
        self._tsv = value
    
    @property
    def fastaDict(self):
        return self._fastaDict
    
    @fastaDict.setter
    def fastaDict(self, value):
        assert isinstance(value, dict), "fastaDict must be a dictionary"
        for k, v in value.items():
            assert os.path.isfile(v), f"FASTA file '{v}' for species '{k}' not found"
        
        self._fastaDict = value
    
    def _parse_groups(self):
        with open(self.tsv , "r") as orthoFile:
            for line in orthoFile:
                l = line.replace('"', '').strip() # remove any oddities from the line
                sl = l.split("\t")
                
                # Handle first line
                if self.species == {}:
                    header = sl[1:] # remove the "Orthogroup" column header
                    self.species = {header[i]: {} for i in range(len(header)) }
                
                # Handle body lines
                else:
                    orthogroupID = sl[0]
                    self.groups[orthogroupID] = {}
                    
                    for i in range(1, len(sl)): # start at 1 to skip the orthogroup ID
                        thisSpecies = header[i-1] # i-1 since header has the orthogroupID removed
                        
                        # Split sequence IDs by comma
                        thisSequences = sl[i].split(", ")
                        if thisSequences == [""]:
                            thisSequences = []
                        thisSequences = set(thisSequences)
                        
                        # Index value in groups subdict
                        self.groups[orthogroupID].setdefault(thisSpecies, set())
                        self.groups[orthogroupID][thisSpecies].update(thisSequences)
                        
                        # Index value in species subdict
                        if thisSequences != set():
                            self.species[thisSpecies][orthogroupID] = thisSequences
    
    def _load_sequences(self):
        sequences = {}
        for species, fasta in self.fastaDict.items():
            sequences[species] = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        
        return sequences
    
    def __repr__(self):
        totalSeqs = sum([len(v) for v in self.species.values()])
        
        reprStr = f"<OrthoGroups object;tsv='{self.tsv}';" + \
                  f"num_species={len(self.species)};" + \
                  f"num_orthogroups={len(self.groups)};" + \
                  f"num_sequences={totalSeqs}>"
        return reprStr

def fastaDict_formatter(fastaDirs, header, suffixes):
    '''
    A helper function to generate the fastaDirs dictionary required
    for the OrthoGroups class.
    
    Parameters:
        fastaDirs -- a list of strings indicating directories to check for
                     FASTA files.
        header -- the header of the Orthogroups.tsv file, and used for
                  checking through the indicated fastaDirs for the relevant files.
        suffixes -- a list of strings indicating the file suffixes to look for.
                    You should use this list to prevent issues with duplicates
                    being found e.g., if you have a species "species1" in your
                    analysis (and, hence, in the header) and you have multiple
                    files like "species1.fasta", "species1.aa" etc., you should
                    specify the suffixes as ["fasta"] if you want to ignore the ".aa"
                    files.
    
    Returns:
        fastaDict -- a dictionary with validated FASTA file locations, and formatted
                     for use with the OrthoGroups class.
    '''
    fastaDict = {}
    for fastaDir in fastaDirs:
        for file in os.listdir(fastaDir):
            filePrefix, fileSuffix = file.rsplit(".", maxsplit=1)
            
            if filePrefix in header and fileSuffix in suffixes:
                assert filePrefix not in fastaDict, f"Duplicate species found: {filePrefix}"
                fastaDict[filePrefix] = os.path.join(fastaDir, file)
    return fastaDict
