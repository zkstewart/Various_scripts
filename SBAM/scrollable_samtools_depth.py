#! python3
# scrollable_samtools_depth.py
# Script to create visualisations of the depth of coverage across a genome
# as calculated by samtools depth. Notably, this script creates line plots
# using the plotly library for scrollable inspection of relative depth
# across each chromosome.

import os, argparse, sys, pickle, math
import plotly.graph_objects as go
import numpy as np

from rdp import rdp

from plot_samtools_depth import get_depth_for_dotting, WMA

##
from typing import Union
##

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

def lang_simplification(input_points: np.ndarray, tolerance: Union[int, float], look_ahead: int) -> np.ndarray:
    """
    Code adapted from https://github.com/adrzystek/lang-simplification-algorithm
    
    Apply the Lang simplification algorithm on a list of points.
    
    :param input_points: points that the algorithm is to be applied on
    :param tolerance: the maximum distance that an intermediate point can be located from the segment
    :param look_ahead: how many points to look ahead while defining a segment
    :return: a subset of `input_points` that is the result of the algorithm run
    """
    output_points = input_points[0]
    key_index = 0
    while True:
        # Look within search region for points within tolerance
        segment = input_points[key_index : key_index + look_ahead + 1]
        points_from_segment_within_tolerance = get_points_within_tolerance(segment, tolerance)
        
        # Add the last point from the segment (which is within tolerance) to the output
        new_key = points_from_segment_within_tolerance[-1]
        output_points = np.column_stack((output_points, new_key))
        
        # Find the index of the new key in the input points list
        key_index = np.where(input_points == new_key)[0][0]
        
        # Exit condition if we reach the end of the input points list
        if key_index == len(input_points) - 1:
            break
    return output_points

def get_points_within_tolerance(points: np.ndarray, tolerance: Union[int, float]) -> np.ndarray:
    """
    This recurrent function returns a subset of `points` that fall within the tolerance limit.

    :param points: points forming a segment
    :param tolerance: the maximum distance between a point and the line between the first and the last point forming a
    segment
    :return: points for which the calculated distance falls below the tolerance level
    """
    points = points.copy()

    start_point = points[0]
    end_point = points[-1]
    intermediate_points = points[1:-1]

    if intermediate_points.shape[0] < 1:
        return points

    intermediate_points_distance = calculate_perpendicular_distance(start_point, end_point, points[1:-1])
    if np.any(intermediate_points_distance > tolerance):
        points = points[:-1]
        return get_points_within_tolerance(points, tolerance)

    return points

def calculate_perpendicular_distance(
    start_point: np.ndarray, end_point: np.ndarray, intermediate_points: np.ndarray,
) -> Union[np.ndarray, float]:
    """
    Calculate the shortest (ie. perpendicular) distance from point(s) to a line drawn between `start_point` and
    `end_point`.

    :param start_point: one of the two points that determine the line location
    :param end_point: another of the two points that determine the line location
    :param intermediate_points: a point or an array of points that the distance is to be calculated for
    :return: either a single distance from an intermediate point to the line or - if the `intermediate_points` is an
    array of more than 1 point - an array of distances
    """
    return np.abs(
        np.cross(end_point - start_point, start_point - intermediate_points) / np.linalg.norm(end_point - start_point)
    )

def plotly_per_contig(dotsX, dotsY,
                      simplificationAlg, wmaSize, 
                      outputDirectory,
                      statisticLabel="Smoothed Q13 depth"):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the depth value for each dot
        simplificationAlg -- a string indicating the simplification algorithm to use;
                             currently recognises "rdp" for Ramer-Douglas-Peucker and
                             "lang" for Lang simplification
        wmaSize -- an integer value indicating the number of previous values to consider
                   during weighted moving average calculation
        outputDirectory -- a string indicating the directory where output files will be written
        statisticLabel -- OPTIONAL; a string indicating what the statistical value represents
                          (default="Q13 depth")
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
        smoothedY = y
        #smoothedY = WMA(y, wmaSize)
        
        # Skip plotting if smoothing fails
        "This probably means there are not enough data points to smooth"
        if smoothedY is None:
            print(f"WARNING: '{contigID}' has too few data points to smooth; skipping...")
            continue
        
        # Simplify data for lighter weight plotting
        xy = np.column_stack((x, smoothedY))
        if simplificationAlg == "rdp":
            xy = rdp(xy)
        elif simplificationAlg == "lang":
            xy = lang_simplification(xy, 5, 20) # tolerance=5, look_ahead=20
        else:
            raise NotImplementedError(f"Unrecognised simplification algorithm '{simplificationAlg}'")
        x, smoothedY = xy[:, 0], xy[:, 1]
        
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
    p.add_argument("--alg", dest="simplificationAlg",
                   required=False,
                   choices=["rdp", "lang"],
                   help="""Optionally, specify the simplification algorithm to use;
                   currently recognises 'rdp' for Ramer-Douglas-Peucker and
                   'lang' for Lang simplification; default is 'lang'""",
                   default="lang")
    
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
    plotly_per_contig(dotsX, dotsY, args.simplificationAlg,
                      args.wmaSize, args.outputDirectory)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
