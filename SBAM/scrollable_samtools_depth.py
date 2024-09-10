#! python3
# scrollable_samtools_depth.py
# Script to create visualisations of the depth of coverage across a genome
# as calculated by samtools depth. Notably, this script creates line plots
# using the plotly library for scrollable inspection of relative depth
# across each chromosome.

import os, argparse, sys, pickle, math
import plotly.graph_objects as go
import numpy as np

from scipy.signal import decimate
from plot_samtools_depth import get_depth_for_dotting, WMA

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.depthFile):
        eprint(f'I am unable to locate the samtools depth TSV file ({args.depthFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        eprint(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        eprint("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        eprint(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")
    # Handle numeric arguments
    if args.wmaSize < 1:
        eprint("wmaSize must be an integer >= 1")
        quit()
    if args.decimateTo < 1:
        eprint("decimateTo must be an integer >= 1")
        quit()

def plotly_per_contig(dotsX, dotsY, wmaSize, 
                      outputDirectory,
                      statisticLabel="Smoothed Q13 depth",
                      decimateTo=100000):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the depth value for each dot
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        outputDirectory -- a string indicating the directory where output files will be written
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
        decimateTo -- OPTIONAL; an integer indicating the number of data points to approximately
                      decimate down to (default=10000)
    '''
    for contigID in dotsX.keys():
        # Derive our output file name and skip if already existing
        fileSuffix = "html"
        fileOut = os.path.join(outputDirectory, f"{contigID}.line.{fileSuffix}")
        if os.path.isfile(fileOut):
            print(f"WARNING: Plotly plot for '{contigID}' already found in output directory; skipping...")
            continue
        
        # Get plotting values
        x = np.array(dotsX[contigID])
        y = np.array(dotsY[contigID])
        smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contigID}' has too few data points to smooth; skipping...")
            continue
        
        # Resample data for lighter weight plotting
        if len(x) > decimateTo:
            xy = np.array([x, smoothedY])
            
            # Decimate data in increments of 10 until we reach our target
            "As recommended by scipy docs"
            decimateFactor = len(x) // decimateTo
            while decimateFactor > 1:
                thisDecimate = min(decimateFactor, 10)
                xy = decimate(xy, thisDecimate, axis=1)
                decimateFactor = len(xy[0]) // decimateTo
            x, smoothedY = xy
        
        # Min-max normalise y values for comparison
        #smoothedY = (smoothedY - np.min(smoothedY)) / (np.max(smoothedY) - np.min(smoothedY))
        
        # Configure plot
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x, y=smoothedY))
        
        # Add range slider
        fig.update_layout(
            ## Set title
            title_text = f"{contigID} {statisticLabel} plot",
            ## Set x-axis slider
            xaxis=dict(
                rangeslider=dict(
                    visible=True
                ),
                type="linear"
            ),
            ## Set y-axis range flexibility
            yaxis=dict(
                autorange=True,
                fixedrange=False
            )
        )
        
        # Save output file
        fig.write_html(fileOut)

def main():
    usage = """%(prog)s receives a samtools depth output file and creates
    an interative plotly line plot of the depth of coverage across each
    chromosome. This can be useful for identifying regions of low coverage
    or high coverage in a genome assembly.
    
    Coverage should be smoothed with --wmaSize to reduce noise, after
    which values will be min-max normalised to enable comparison of
    relative coverage across chromosomes.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-d", dest="depthFile",
                    required=True,
                    help="Specify the location of the input depth TSV file")
    p.add_argument("-o", dest="outputDirectory",
                    required=True,
                    help="Output directory where plot files will be written")
    ## Opts (behavioural)
    p.add_argument("--wmaSize", dest="wmaSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the number of previous values to consider
                   during weighted moving average calculation (default=20)""",
                   default=20)
    p.add_argument("--decimateTo", dest="decimateTo",
                   type=int,
                   required=False,
                   help="""Optionally, specify the maximum number of data points to allow
                   in an output plot (default=100000)""",
                   default=100000)
    
    args = p.parse_args()
    validate_args(args)
    
    # Figure out what our pickle file should be called
    pickleFile = os.path.join(
        args.outputDirectory,
        f"{os.path.basename(args.depthFile)}.pkl"
    )
    
    # Load pickle if it exists to skip computation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            dotsX, dotsY = pickle.load(fileIn)
    
    # Otherwise ...
    else:
        # Parse samtools depth data
        dotsX, dotsY = get_depth_for_dotting(args.depthFile)
        
        # Save data
        with open(pickleFile, "wb") as fileOut:
            pickle.dump([dotsX, dotsY], fileOut)
    
    # Create plots
    plotly_per_contig(dotsX, dotsY, args.wmaSize,
                      args.outputDirectory)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
