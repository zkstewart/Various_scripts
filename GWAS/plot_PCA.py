import os, argparse
import numpy as np
import plotly.express as px

def validate_args(args):
    # Validate input file locations
    args.sscoreFile = os.path.abspath(args.sscoreFile)
    if not os.path.isfile(args.sscoreFile):
        raise FileNotFoundError(f"-s value '{args.sscoreFile}' is not a file or could not be located.")
    
    args.eigenvalFile = os.path.abspath(args.eigenvalFile)
    if not os.path.isfile(args.eigenvalFile): 
        raise FileNotFoundError(f"-e value '{args.eigenvalFile}' is not a file or could not be located.")
    
    # Validate output file
    if not args.outputFileName.endswith(".html"):
        raise ValueError(f"-o value '{args.outputFileName}' must end in .html")
    
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o value '{args.outputFileName}' already exists and will not be overwritten.")
    args.outputFileName = os.path.abspath(args.outputFileName)

def parse_sscore(fileName, iidCol="IID", fidCol="#FID"):
    '''
    Parse the PLINK2 .sscore file into a dictionary with inference of corresponding phenotype measurement types.
    
    Returns:
        sscoreDict -- a dictionary where keys are sample IDs and values give the family ID (fid)
                      along with any phenotypes (PHENO*) and principal components (PC*)
        phenoTypes -- a list containing booleans indicating whether each phenotype (in ascending
                      numerical order) is a categorical value (True) or is a continuous
                      measurement (False)
    '''
    sscoreDict = {}
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if firstLine:
                # Identify user-specifiable column indices
                try:
                    iidIndex = sl.index(iidCol)
                except ValueError:
                    raise ValueError(f"--sampleid '{iidCol}' does not appear in the sscore header '{sl}'")
                try:
                    fidIndex = sl.index(fidCol)
                except ValueError:
                    raise ValueError(f"--familyid '{fidCol}' does not appear in the sscore header '{sl}'")
                
                # Identify any other columns which may exist in variable number
                phenoIndices = [ i for i, value in enumerate(sl) if value.startswith("PHENO") ]
                pcIndices = [ i for i, value in enumerate(sl) if value.startswith("PC") and value.endswith("_AVG") ]
                
                firstLine = False
            else:
                sl = line.rstrip().split("\t")
                sampleID, familyID = sl[iidIndex], sl[fidIndex]
                phenos = [ sl[x] for x in phenoIndices ]
                pcs = [ float(sl[x]) for x in pcIndices ]
                sscoreDict[sampleID] = { "fid": familyID, "phenos": phenos, "pcs": pcs }
    
    # Figure out which PHENO attributes are categorical versus quantitative
    phenoTypes = []
    for i in range(len(phenoIndices)):
        phenoValues = [ v["phenos"][i] for v in sscoreDict.values() ]
        isCategorical = all([ pValue.isdigit() for pValue in phenoValues ])
        phenoTypes.append(isCategorical)
    
    return sscoreDict, phenoTypes

def parse_eigenval(fileName):
    '''
    Parses the .eigenval file produced by PLINK2 to identify the amount of variance explained
    by each component.
    
    Returns:
        explained -- a list of floats indicating the proportion (fraction of 1) of variance
                     explained by each component in ascending numerical order.
    '''
    eigenvals = []
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            val = line.rstrip()
            if val != "":
                eigenvals.append(float(val))
    total = sum(eigenvals)
    explained = [ x/total for x in eigenvals ]
    return explained

def main():
    usage = """%(prog)s receives a .sscore file produced by PLINK2 and visualises the principal components
    to allow identification of how many PCs are needed to explain any variance relevant to a e.g., subsequent
    GWAS analysis. Note that the --sampleid and --familyid values can be configured if the .sscore format has
    slight differences across versions, but its default is set to the PLINK2 behaviour as of this script's writing.
    """
    # Required arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-s", dest="sscoreFile",
                   required=True,
                   help="Location of .sscore file")
    p.add_argument("-e", dest="eigenvalFile",
                   required=True,
                   help="Location of .eigenval file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write output plot; file must end in .html")
    # Optional arguments
    p.add_argument("--sampleid", dest="sampleIDCol",
                   required=False,
                   help="""Optionally, specify the column header that denotes the sample IDs;
                   default=='IID'""",
                   default="IID")
    p.add_argument("--familyid", dest="familyIDCol",
                   required=False,
                   help="""Optionally, specify the column header that denotes the family IDs;
                   default=='#FID'""",
                   default="FID")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse input files
    sscoreDict, phenoTypes = parse_sscore(args.sscoreFile, iidCol="IID", fidCol="#FID")
    explained = parse_eigenval(args.eigenvalFile)
    sampleOrder = sorted(sscoreDict.keys())
    
    # Convert sscoreDict into a numpy array amenable to plotly handling
    components = np.array([
        sscoreDict[sampleID]["pcs"]
        for sampleID in sampleOrder
    ])
    
    # Extract and format aesthetics for the plot
    pclabels = {
        str(i): f"PC {i+1} ({var*100:.1f}%)"
        for i, var in enumerate(explained)
    }
    families = [
        sscoreDict[sampleID]["fid"]
        for sampleID in sampleOrder
    ]
    pcdimensions = len(explained)
    
    # Generate the plot
    fig = px.scatter_matrix(
        components,
        labels=pclabels,
        dimensions=range(pcdimensions),
        color=families,
        hover_name=sampleOrder
    )
    fig.update_traces(diagonal_visible=False)
    fig.write_html(args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
