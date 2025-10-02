#! python3

import os, argparse
from goatools import obo_parser

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputTextFile):
                print('I am unable to locate the input text file (' + args.inputTextFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        for oboFile in args.oboFile:
                if not os.path.isfile(oboFile):
                        print('I am unable to locate the go-basic.obo file (' + oboFile + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
                quit()

def text_go_parse(textFile, blankCharacter):
        # Setup
        codeDict = {}
        # Main function
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        line = line.rstrip('\r\n')
                        if "\t" in line:
                                line = line.split("\t")[1]
                        if line == blankCharacter:
                                continue
                        else:
                                codes = line.split('; ')
                                for code in codes:
                                        codeDict[code] = ''
        return codeDict

#### USER INPUT SECTION
usage = """This program will read in a text file listing GO codes separated with '; '
characters alongside the go-basic.obo file and produces an output text file listing the
names for these ontologies. '.' or '0' can be used to mark lines with no GO results, and this
will be carried over into the output file. It will also tolerate a TSV format where the second
column contains the GO codes (the first column is ignored here).
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="inputTextFile",
               required=True,
               help="Input text file name listing GO codes")
p.add_argument("-g", "-goobo", dest="oboFile",
               required=True,
               nargs="+",
               help="Input go-basic.obo file.")
p.add_argument("-o", "-output", dest="outputFileName",
               required=True,
               help="Output text file name.")
# Optional
p.add_argument("--blank", dest="blankCharacter",
               choices=[".", "0"],
               help="Specify which character marks a blank (default=='.').",
               default=".")

args = p.parse_args()
validate_args(args)

# Parse the input text file to identify which GO codes are listed
codeDict = text_go_parse(args.inputTextFile, args.blankCharacter)

# Parse .obo file
gos = [ obo_parser.GODag(oboFile) for oboFile in args.oboFile ]

# Hardcoded GO replacements for problem entries
obsoletedGOs = ['GO:0097034', 'GO:0097033']   # I can't find any suitable replacements for these terms
replacedGOs = {
    'GO:0140603': 'GO:0016887', 'GO:0036425': 'GO:0036424','GO:0005671': 'GO:0140671',
    'GO:2000574': 'GO:0140659', 'GO:0102132': 'GO:0004316', 'GO:0102131': 'GO:0004316',
    'GO:0009405': 'GO:0052031', 'GO:0015002': 'GO:0016491', 'GO:0052331': 'GO:0044179',
    'GO:0000186': 'GO:0000165', 'GO:2000575': 'GO:0140661', 'GO:0015491': 'GO:0008324',
    'GO:0006557': 'GO:0004014', 'GO:0000187': 'GO:0000165', 'GO:0018298': 'GO:0043687',
    'GO:0033577': 'GO:0006486', 'GO:0102552': 'GO:0016992', 'GO:2000582': 'GO:0140660',
    'GO:0031936': 'GO:0031047', 'GO:0102553': 'GO:0016992', 'GO:0070827': 'GO:0006325',
    'GO:0060968': 'GO:0040029', 'GO:0045857': 'GO:0044092', 'GO:0042766': 'GO:0140658',
    'GO:0009305': 'GO:0036211', 'GO:0070122': 'GO:0008233', 'GO:0010847': 'GO:0006325',
    'GO:0016584': 'GO:0140658', 'GO:0015301': 'GO:0008509', 'GO:0031935': 'GO:1902275',
    'GO:0031938': 'GO:0031509', 'GO:0005639': 'GO:0005637', 'GO:0005779': 'GO:0005778',
    'GO:0006471': 'GO:1990404', 'GO:0030173': 'GO:0000139', 'GO:0030176': 'GO:0005789',
    'GO:0031225': 'GO:0016020', 'GO:0031227': 'GO:0005789', 'GO:0031305': 'GO:0005743',
    'GO:0031307': 'GO:0005741', 'GO:0031362': 'GO:0009897', 'GO:0032592': 'GO:0031966',
    'GO:0043004': 'GO:0051220', 'GO:0043486': 'GO:0006338', 'GO:0046658': 'GO:0005886',
    'GO:0031224': 'GO:0016020', 'GO:0031226': 'GO:0005886', 'GO:0102488': 'GO:0017111',
    'GO:0102486': 'GO:0017111', 'GO:0102491': 'GO:0017111', 'GO:0102487': 'GO:0017111',
    'GO:0102490': 'GO:0017111', 'GO:0102489': 'GO:0017111', 'GO:0102485': 'GO:0017111',
    'GO:0031301': 'GO:0031090', 'GO:0031300': 'GO:0031090', 'GO:0018196': 'GO:0018193',
    'GO:0031228': 'GO:0000139', 'GO:0098573': 'GO:0031966', 'GO:0031358': 'GO:0009707',
    'GO:0031359': 'GO:0009707', 'GO:0031350': 'GO:0042170', 'GO:0031355': 'GO:0009527',
    'GO:0031354': 'GO:0009527', 'GO:0031351': 'GO:0042170', 'GO:0071556': 'GO:0098553',
    'GO:0071458': 'GO:0098554', 'GO:0000453': 'GO:0006364', 'GO:0015299': 'GO:0015078',
    'GO:0030285': 'GO:0030672', 'GO:0044214': 'GO:0005886', 'GO:0097056': 'GO:0001717',
    'GO:0016277': 'GO:0016274', 'GO:0031357': 'GO:0009706', 'GO:0031361': 'GO:0042651',
    'GO:1900049': 'GO:0140713', 'GO:1904576': 'GO:0034976', 'GO:0140323': 'GO:0008509',
    'GO:0055065': 'GO:0030003', 'GO:0072507': 'GO:0055080', 'GO:0072503': 'GO:0030003',
    'GO:0006875': 'GO:0030003', 'GO:0046916': 'GO:0030003', 'GO:0055076': 'GO:0030003',
    'GO:0055067': 'GO:0055080', 'GO:0030004': 'GO:0030003', 'GO:0072506': 'GO:0055081',
    'GO:0015298': 'GO:0008324', 'GO:0005451': 'GO:0008324', 'GO:0035511': 'GO:0035516',
    'GO:0044728': 'GO:0006304', 'GO:0033683': 'GO:0006289', 'GO:0008022': 'GO:0005515',
    'GO:0043631': 'GO:0016070', 'GO:0016307': 'GO:0052742', 'GO:0102756': 'GO:0009922',
    'GO:0035004': 'GO:0016303', 'GO:0000904': 'GO:0000902', 'GO:0017182': 'GO:0018202',
    'GO:0048017': 'GO:0035556', 'GO:1905114': 'GO:0007166', 'GO:0010769': 'GO:0022604',
    'GO:1900363': 'GO:0006378', 
} # Modifications to this point made 17-01-23

# Produce output text file of GO descriptions
with open(args.inputTextFile, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                line = line.rstrip('\r\n')
                isTsv = False
                if "\t" in line:
                        geneID, line = line.split("\t")
                        isTsv = True
                # Handle blank lines
                if line == args.blankCharacter:
                        if isTsv:
                                fileOut.write(f'{geneID}\t{args.blankCharacter}\n')
                        else:
                                fileOut.write(f'{args.blankCharacter}\n')
                # Handle GO code lines
                else:
                        # Get the names of GO entries
                        codes = line.split('; ')
                        names = []
                        for code in codes:
                                if code in replacedGOs:
                                        code = replacedGOs[code]
                                foundGo = False
                                for go in gos:
                                        try:
                                                names.append(go[code].name)
                                                foundGo = True
                                                break
                                        except:
                                                pass
                                assert foundGo == True, \
                                    f'GO code {code} not found in any of the obo files!'
                        # Format and write to file
                        if isTsv:
                                fileOut.write(geneID + "\t" + '; '.join(names) + '\n')
                        else:
                                fileOut.write('; '.join(names) + '\n')

# Done!
print('Program completed successfully!')
