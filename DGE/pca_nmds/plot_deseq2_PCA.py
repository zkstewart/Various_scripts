import os, argparse
import numpy as np
import pandas as pd
import plotly.express as px

def validate_args(args):
    # Validate input file locations
    args.pcaFile = os.path.abspath(args.pcaFile)
    if not os.path.isfile(args.pcaFile):
        raise FileNotFoundError(f"-i value '{args.pcaFile}' is not a file or could not be located.")
    
    # Validate output file
    if not args.outputFileName.endswith(".html"):
        raise ValueError(f"-o value '{args.outputFileName}' must end in .html")
    
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o value '{args.outputFileName}' already exists and will not be overwritten.")
    args.outputFileName = os.path.abspath(args.outputFileName)

def parse_pca_table(fileName, nameColumn="name"):
    '''
    Parse the DESeq2 produced table file into a dictionary, while ignoring the 'group' and
    the 'nameColumn' columns.
    
    Returns:
        pcaDict -- a dictionary where keys are sample IDs and values give the phenotypes/treatment
                   metadata along with the principal components (PC*)
    '''
    pcaDict = {}
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if firstLine:
                # Identify user-specifiable column indices
                try:
                    nameIndex = sl.index(nameColumn)
                except ValueError:
                    raise ValueError(f"--name '{nameColumn}' does not appear in the table header '{sl}'")
                
                # Identify the PC columns
                pcHeads = [ value for value in sl if value.startswith("PC") and value[2:].isdigit() ]
                pcHeads.sort(key = lambda x: int(x[2:]))
                pcIndices = [ sl.index(x) for x in pcHeads ]
                
                # Identify any phenotype/metadata columns
                "metadata can be found by exclusion without knowing what the column header values will be"
                phenotypes = [
                    (i, value)
                    for i, value in enumerate(sl)
                    if (not i == 0) and (not i == nameIndex) and (not i in pcIndices) and (not value == "group")
                ]
                phenoIndices, phenoLabels = zip(*phenotypes)
                
                firstLine = False
            else:
                sl = line.rstrip().split("\t")
                sampleID = sl[nameIndex]
                phenos = [ sl[x] for x in phenoIndices ]
                pcs = [ float(sl[x]) for x in pcIndices ]
                pcaDict[sampleID] = { "phenos": phenos, "pcs": pcs }
    
    return pcaDict, phenoLabels

def main():
    usage = """%(prog)s receives a .tsv file produced from DESeq2's plotPCA() function and visualises
    the principal components along with the sample metadata to assist in identifying any relevant trends. 
    """
    # Required arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="pcaFile",
                   required=True,
                   help="Location of PCA .tsv file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write output plot; file must end in .html")
    # Optional arguments
    p.add_argument("--name", dest="nameColumn",
                   required=False,
                   help="""Optionally, specify the column header that denotes the sample IDs;
                   default=='name'""",
                   default="name")
    p.add_argument("--colour", dest="colourColumn",
                   required=False,
                   help="""Optionally, specify the column header that denotes phenotype to colour;
                   if unspecified, no colour will be used""",
                   default=None)
    p.add_argument("--pcs", dest="pcsToPlot",
                   required=False,
                   nargs="+",
                   type=int,
                   help="""Optionally, specify the PCs by their number (e.g., '--pcs 1 2 3') to limit
                   plotting to just those components. The first PC is '1'. Default is to plot all PCs found in the
                   PCA file""",
                   default=[])
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse input file
    pcaDict, phenoLabels = parse_pca_table(args.pcaFile, nameColumn=args.nameColumn)
    sampleOrder = sorted(pcaDict.keys())
    
    # Convert pcaDict into a pandas DataFrame amenable to plotly handling
    components = np.array([
        pcaDict[sampleID]["pcs"]
        for sampleID in sampleOrder
    ])
    df = pd.DataFrame(components)
    
    # Add extra columns for metadata
    for i, pheno in enumerate(phenoLabels):
        df[pheno] = [
            pcaDict[sampleID]["phenos"][i]
            for sampleID in sampleOrder
        ]
    
    # Extract and format aesthetics for the plot
    exampleSampleID = list(pcaDict.keys())[0]
    pclabels = {
        str(i): f"PC {i+1}"
        for i, var in enumerate(pcaDict[exampleSampleID]["pcs"])
    }
    if args.colourColumn != None:
        if not args.colourColumn in phenoLabels:
            raise ValueError(f"--colour value '{args.colourColumn}' not found in PCA metadata columns ({phenoLabels})")
        
        pclabels["color"] = args.colourColumn
        
        colourGroups = [
            pcaDict[sampleID]["phenos"][phenoLabels.index(args.colourColumn)]
            for sampleID in sampleOrder
        ]
    
    if args.pcsToPlot == []:
        pcdimensions = range(len(pcaDict[exampleSampleID]["pcs"]))
    else:
        pcdimensions = [ x - 1 for x in args.pcsToPlot ]
        if any([ x < 0 for x in pcdimensions ]):
            raise ValueError("The lowest value given to --pcs should be 1")
    
    # Generate the plot
    if len(pcdimensions) > 2:
        fig = px.scatter_matrix(
            df,
            labels=pclabels,
            dimensions=pcdimensions,
            color=colourGroups,
            hover_name=sampleOrder,
            hover_data=phenoLabels
        )
        fig.update_traces(diagonal_visible=False)
    else:
        fig = px.scatter(
            df,
            labels=pclabels,
            x=pcdimensions[0],
            y=pcdimensions[1],
            color=colourGroups,
            hover_name=sampleOrder
        )
    
    fig.write_html(args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
