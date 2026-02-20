#!/usr/bin/env python3

import os, re, pyperclip

# Paste document contents
wait = input("Copy the document text then press ENTER")
document = pyperclip.paste()

# Paste citations from reference manager
wait = input("Copy the in-text citations (e.g., from Zotero) then press ENTER")
refmanager = pyperclip.paste()

# Identify in-text citations from the document
parenthesesRegex = re.compile(r"\(.+?\)")
yearRegex = re.compile(r",\s(1|2)\d{3}$")
multiyearRegex = re.compile(r",\s(1|2)\d{3}\/(1|2)\d{3}$")

intextSet = set()
for parentheses in parenthesesRegex.findall(document):
    components = [ x.strip("() ") for x in parentheses.split(";") ]
    for component in components:
        if yearRegex.search(component) or multiyearRegex.search(component):
            intextSet.add(component.replace("’", "'")) # get rid of word-formatted characters

# Format in-text citations from reference manager as an equivalent set
multiplesRegex = re.compile(r",\s\d{4}")

tmpcitations = [ x.strip("() ").replace("’", "'") for x in refmanager.split(";") ] # get rid of word-formatted characters
citationsSet = set()
for c in tmpcitations:
    years = multiplesRegex.findall(c)
    if len(years) == 1:
        citationsSet.add(c)
    else:
        author = c.split(years[0])[0]
        for year in years:
            citationsSet.add(f"{author}{year}")

# Find overlaps
common = citationsSet.intersection(intextSet)
intextOnly = intextSet.difference(common)
citationsOnly = citationsSet.difference(common)

# Find explainable differences
def found_within(set1, set2):
    found = set()
    for s1 in set1:
        for s2 in set2:
            s1mod = s1.replace("&", "and")
            s2mod = s2.replace("&", "and")
            
            if ((s1 in s2) or (s2 in s1)) or ((s1mod in s2mod) or (s2mod in s1mod)):
                found.add(s1)
                found.add(s2)
    return set1.difference(found), set2.difference(found)

citationsOnly, intextOnly = found_within(citationsOnly, intextOnly)

# Print out results message
if len(citationsOnly) == 0 and len(intextOnly) == 0:
    print("# Success! No discrepancies noticed between document and reference manager")
else:
    print("# Discrepancies noticed between document and reference manager\n")
    if len(intextOnly) != 0:
        print("## Document-only references include:")
        for value in sorted(intextOnly):
            print(value)
        print("")
    if len(citationsOnly) != 0:
        print("## Reference manager-only references include:")
        for value in sorted(citationsOnly):
            print(value)
