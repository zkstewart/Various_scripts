#! python3

import re

def parse_change_tab(tabFileName):
    '''
    Function to parse the *_change.tab file from CAFE to obtain the changes occurring
    at each node in the tree.
    
    Parameters:
        tabFileName -- a string representing the path to the *_change.tab file
    Returns:
        changeDict -- a dictionary with structure like:
                      {
                          "1": {
                              "Family1": changeNum,
                              "Family2": changeNum,
                              ...
                          },
                          "2": { ... },
                          ...
                      }
    '''
    nodeRegex = re.compile(r"<(\d+)>")
    
    changeDict = {}
    with open(tabFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            # Handle header line
            if sl[0] == "FamilyID":
                nodeNames = sl[1:]
                assert all([ nodeRegex.search(x) for x in nodeNames ]), \
                    "All node names must contain '<number>'"
                
                nodeNums = [ nodeRegex.search(x).groups()[0] for x in nodeNames ]
                for num in nodeNums:
                    changeDict[num] = {}
            # Handle other lines
            else:
                familyName, counts = sl[0], list(map(int, sl[1:]))
                for i, num in enumerate(nodeNums):
                    changeDict[num][familyName] = counts[i]
    return changeDict

def parse_probabilities_tab(tabFileName, significance=0.05):
    '''
    Function to parse the *_branch_probabilities.tab file from CAFE to identify which
    families change in a statistically significant way at each node in the tree.
    
    Parameters:
        tabFileName -- a string representing the path to the *_branch_probabilities.tab file
        significance -- OPTIONAL; a float indicating the p-value threshold for significance;
                        default == 0.05
    Returns:
        significantDict -- a dictionary with structure like:
                           {
                               "1": set(["Family1", "Family2", ...]),
                               "2": { ... },
                               ...
                           }
    '''
    nodeRegex = re.compile(r"<(\d+)>")
    
    significantDict = {}
    with open(tabFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            # Handle header line
            if sl[0] == "FamilyID":
                nodeNames = sl[1:]
                assert all([ nodeRegex.search(x) for x in nodeNames ]), \
                    "All node names must contain '<number>'"
                
                nodeNums = [ nodeRegex.search(x).groups()[0] for x in nodeNames ]
                for num in nodeNums:
                    significantDict[num] = set()
            # Handle other lines
            else:
                familyName, pvalues = sl[0], sl[1:]
                for i, num in enumerate(nodeNums):
                    pvalue = pvalues[i]
                    if pvalue == "N/A":
                        continue
                    
                    if float(pvalue) <= significance:
                        significantDict[num].add(familyName)
    return significantDict

if __name__ == "__main__":
    pass
