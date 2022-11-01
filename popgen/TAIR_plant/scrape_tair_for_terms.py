#! python3
# scrape_tair_for_terms.py
# Script to take a TSV file containing columns which have TAIR
# URLs, and scrape each page's contents for certain terms

import os, argparse, requests, re, time
from bs4 import BeautifulSoup

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.tsvFileName):
        print('I am unable to locate the TSV file (' + args.tsvFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def parse_urls_from_tsv(tsvFile, urlColumnName):
    tairURLs = []
    with open(tsvFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header line
            if firstLine is True:
                if urlColumnName not in sl:
                    raise ValueError(f"{urlColumnName} not found in first line of TSV file")
                
                urlIndex = sl.index(urlColumnName)
                firstLine = False
                continue
            # Handle content lines
            else:
                url = sl[urlIndex]
                tairURLs.append(url)
    return tairURLs

def retry_scrape(url, savedHTML):
    # Get web page HTML
    #headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36'}
    headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36'}
    contents = None
    while contents == None:
        try:
            page = requests.get(url, headers=headers)
            contents = page.text
            with open(savedHTML, "w") as fileOut:
                fileOut.write(contents)
            
            # Courtesy sleep to reduce load on server
            time.sleep(10)
        except:
            print("Connection refused; sleeping for 1 min")
            time.sleep(60)
    
    return contents

def clean_mention(mention):
    return mention.strip(" \t").replace("\r", "").replace("\n", "")

def scrape_tair_for_term(tairURL, terms, htmlsDir):
    termMentions = {term: [] for term in terms}
    
    # Get web page HTML
    savedHTML = os.path.join(htmlsDir, tairURL.split("name=")[1])
    if os.path.isfile(savedHTML):
        contents = open(savedHTML, "r").read()
    else:
        contents = retry_scrape(tairURL, savedHTML)
    
    # Scrape this page for terms
    uniqueMentions = set() # handles times where more than 1 term shows up in the same section
    for term in terms:
        if term in contents:
            mentions = re.findall(r">([^<>]*?" + term + r".*?)<", contents.lower(), re.DOTALL)
            for m in mentions:
                if m not in uniqueMentions:
                    termMentions[term].append(clean_mention(m))
                    uniqueMentions.add(m)
    
    # Locate publication links
    publicationUrlRegex = re.compile(r"<a href=(/servlets/TairObject\?type=publication.+?id=.+?)>", re.DOTALL)
    publicationURLs = [f"https://www.arabidopsis.org{urlSuffix}" for urlSuffix in publicationUrlRegex.findall(contents)]
    
    # Scrape each publication link for terms
    for pubUrl in publicationURLs:
        # Get web page HTML
        savedPubHTML = os.path.join(htmlsDir, pubUrl.split("id=")[1])
        
        if os.path.isfile(savedPubHTML):
            pubContents = open(savedPubHTML, "r").read()
        else:
            pubContents = retry_scrape(pubUrl, savedPubHTML)
        
        # Parse for terms
        for term in terms:
            if term in pubContents:
                mentions = re.findall(r">([^<>]*?" + term + r".*?)<", pubContents.lower(), re.DOTALL)
                for m in mentions:
                    if m not in uniqueMentions:
                        termMentions[term].append(clean_mention(m))
                        uniqueMentions.add(m)
    
    return termMentions

def get_scraping_results(tairURLs, terms, htmlsDir):
    scrapeResults = []
    for tairURL in tairURLs:
        # Skip null rows
        if tairURL == ".":
            scrapeResults.append(["."]*len(terms))
            continue
        
        # Get TAIR mentions
        termMentions = scrape_tair_for_term(tairURL, terms, htmlsDir)
        
        # Skip pages with no mentions
        if [m for mentions in list(termMentions.values()) for m in mentions] == []:
            scrapeResults.append(["."]*len(terms))
            continue
        # Format results for pages with mentions
        else:
            thisResult = []
            for term in terms:
                if termMentions[term] == []:
                    thisResult.append(".")
                else:
                    thisResult.append(" ... ".join(termMentions[term]))
            
            scrapeResults.append(thisResult)
    
    return scrapeResults

def main():
    # User input
    usage = """%(prog)s receives a TSV file containing annotations wherein one column
    lists a TAIR URL (or "." for null). That URL will be scraped looking for one or more
    terms that are user-specifiable. Note that the TSV file must have a header line (first
    line, with or without a #) and you need to specify which column contains the URLs.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="tsvFileName",
                   required=True,
                   help="Input TSV file produced by QTL-seq")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    p.add_argument("-c", dest="columnName",
                   required=True,
                   help="Specify the header for the column containing TAIR URLs")
    p.add_argument("--terms", dest="terms", nargs="+",
                   required=True,
                   help="Specify one or more terms to look for in the TAIR webpages")
    args = p.parse_args()
    validate_args(args)
    
    # Set up location for HTML saving to limit scraping impact
    htmlsDir = os.path.join(os.path.dirname(args.outputFileName), "scraped_htmls")
    os.makedirs(htmlsDir, exist_ok=True)
    
    # Parse TSV file
    tairURLs = parse_urls_from_tsv(args.tsvFileName, args.columnName)
    
    # Perform scraping of TAIR website
    scrapeResults = get_scraping_results(tairURLs, args.terms, htmlsDir)
    
    # Format new output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("#{0}\n".format("\t".join([t + "_scrape" for t in args.terms])))
        for result in scrapeResults:
            fileOut.write("{0}\n".format("\t".join(result)))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
