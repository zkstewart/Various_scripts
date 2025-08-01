#! python3

import os, argparse, requests, re, untangle, time

CONVERSION_FILE = 'dbto_dbfrom_conversion.tsv'
FASTA_FILE = 'sequences.fasta'

class NCBI_API:
    def __init__(self, email, apiKey):
        self.email = email
        self.apiKey = apiKey
        
        self.base = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.link = "elink.fcgi?"
        self.fetch = "efetch.fcgi?"
    
    def url_key(self):
        return f"&api_key={self.apiKey}"
    
    def link_dbs(self, dbfrom, dbto, identifier):
        '''
        Parameters:
            dbfrom -- a string indicating the NCBI database that the identifier
                      belongs to
            dbto -- a string indicating the NCBI database that you want the identifier
                    to be converted to
            identifier -- the string identifier from dbfrom that you want to convert
                          to its corresponding dbto identifier
        Returns:
            linkedID -- the converted identifier string for the dbto database
        '''
        url = f"{self.base}{self.link}dbfrom={dbfrom}&db={dbto}&id={identifier}"#{self.url_key}"
        response = requests.get(url)
        
        xmlText = re.sub(r"<!DOCTYPE.+?>", "", response.text)
        document = untangle.parse(xmlText)
        try:
            linkedID = document.eLinkResult.LinkSet.LinkSetDb.Link.Id.cdata
        except:
            linkedID = None
        
        return linkedID
    
    def fetch_fasta(self, dbfrom, identifier):
        '''
        Parameters:
            dbfrom -- a string indicating the NCBI database that the identifier
                      belongs to
            identifier -- the string identifier that you want to grab a FASTA for
        Returns:
            fastaString -- a multiline string for the FASTA
        '''
        url = f"{self.base}{self.fetch}db={dbfrom}&id={identifier}&rettype=fasta&retmode=text"#{self.url_key}"
        response = requests.get(url)
        
        fastaString = response.text.rstrip("\n")
        
        return fastaString
    
    def batch_fetch(self, dbfrom, dbto, identifiers):
        '''
        Parameters:
            dbfrom -- a string indicating the NCBI database that the identifier
                      belongs to
            dbto -- a string indicating the NCBI database that you want the identifier
                    to be converted to
            identifiers -- a list of the string identifiers from dbfrom that you want to
                           grab FASTA sequences from dbto for
        Returns:
            linkedIDs -- a list of converted identifier strings for the dbto database
            fastaStrings -- a list of multiline strings for the FASTA
        '''
        linkedIDs, fastaStrings = [], []
        for identifier in identifiers:
            if dbto != None and dbfrom != dbto:
                linkedID = self.link_dbs(dbfrom, dbto, identifier)
            else:
                linkedID = identifier
            
            if linkedID == None:
                print(f"Failed to link '{identifier}'; skipping...")
                linkedIDs.append(f"{identifier}\t.")
                continue
            linkedIDs.append(f"{identifier}\t{linkedID}")
            
            fastaString = self.fetch_fasta(dbto, linkedID)
            fastaStrings.append(fastaString)
            
            time.sleep(1) # no more than 3 API queries per second without API key
        return linkedIDs, fastaStrings

def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.identifiersFile):
        raise FileNotFoundError(f"-i '{args.identifiersFile} does not exist!")
    # Validate output file locations
    args.conversionFile = os.path.join(os.path.abspath(args.outputDirectory), CONVERSION_FILE)
    args.fastaFile = os.path.join(os.path.abspath(args.outputDirectory), FASTA_FILE)
    
    if os.path.exists(args.outputDirectory):
        if os.path.exists(args.conversionFile):
            raise FileExistsError(f"'{args.conversionFile}' file already exists; set a different -o location")
        
        if os.path.exists(args.fastaFile):
            raise FileExistsError(f"'{args.fastaFile}' file already exists; set a different -o location")
    else:
        os.makedirs(args.outputDirectory)
        print(f"Created output directory '{args.outputDirectory}' as " + 
              "part of argument validation")

def parse_identifiers(identifiersFile):
    identifiers = []
    with open(identifiersFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n\t ")
            if l != "":
                identifiers.append(l)
    return identifiers

def main():
    usage = """%(prog)s will receive a list of Entrez sequence identifiers,
    optionally convert them from one database type to another, and fetch
    them as FASTA sequences. Outputs will be written to a directory,
    which will be populated with the 'dbto_dbfrom_conversion.tsv' and
    'sequences.fasta' files
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="identifiersFile",
                   required=True,
                   help="Specify the input file listing sequence identifiers")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="""Output directory where FASTA and conversion list files
                   are to be written""")
    # Optional
    p.add_argument("--dbfrom", dest="dbfrom",
                   required=False,
                   help="""Specify the database your identifiers belong to;
                   default == 'nuccore'""",
                   default="nuccore")
    p.add_argument("--dbto", dest="dbto",
                   required=False,
                   help="""Specify the database to convert identifiers to prior
                   to FASTA fetching; if unspecified, no conversion will occur""",
                   default=None)
    
    args = p.parse_args() # sets .conversionFile; .fastaFile
    validate_args(args)
    
    # Parse identifiers file
    identifiers = parse_identifiers(args.identifiersFile)
    
    # Create API object
    api = NCBI_API(None, None) # email, apiKey; not implemented
    
    # Run batch query
    linkedIDs, fastaStrings = api.batch_fetch(args.dbfrom, args.dbto, identifiers)
    
    # Write outputs to file
    with open(args.conversionFile, "w") as fileOut:
        for value in linkedIDs:
            fileOut.write(f"{value}\n")
    
    with open(args.fastaFile, "w") as fileOut:
        for value in fastaStrings:
            fileOut.write(f"{value}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
