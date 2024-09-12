#! python3
# scrollable_samtools_depth.py
# Script to create visualisations of the depth of coverage across a genome
# as calculated by samtools depth. Notably, this script creates line plots
# using the plotly library for scrollable inspection of relative depth
# across each chromosome.

import os, argparse, sys, pickle
import plotly.graph_objects as go
import numpy as np

from rdp import rdp
from typing import Union

from plot_samtools_depth import get_depth_for_dotting

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

def simplify_depth(xy: np.ndarray, tolerance: Union[int, float]) -> np.ndarray:
    '''
    Simplify depth data for plotting by removing data points
    that are redundant or minimally informative.
    
    Algorithm works by obtaining the first and last points, in addition to
    any points that mark a departure from the previously stored point by more
    than the specified tolerance.
    
    Parameters:
        xy -- a 2D numpy array where the first column is the x-axis position
              and the second column is the y-axis value
        tolerance -- an integer or float value indicating the maximum increase
                     or decrease in y-axis value before a change is detected
                     and the point is retained
    Returns:
        simplifiedXY -- a 2D numpy array of the input data with some points removed
                        where deemed to be redundant or minimally informative
    '''
    def process_plateau(simplifiedXY, xPlateau, yPlateau):
        # Process plateaus with only one point
        if len(xPlateau) == 1:
            simplifiedXY.append([xPlateau[0], yPlateau[0]])
        
        # Process plateaus with multiple points
        else:
            # Obtain and store the median value of the plateau
            medianY = np.median(yPlateau)
            
            # Store the plateau's median value
            simplifiedXY.append([xPlateau[0], medianY])
            simplifiedXY.append([x-1, medianY])
    
    # Skip simplification if there are too few points
    if len(xy) < 3:
        return xy
    
    # Note first data points
    simplifiedXY = []
    xPlateau = [xy[0][0]]
    yPlateau = [xy[0][1]]
    
    # Begin iteration through remaining data points
    for i in range(1, len(xy)):
        x, y = xy[i]
        
        # Check for departure from this plateau's tolerance threshold
        if y > yPlateau[0] + tolerance or y < yPlateau[0] - tolerance:
            process_plateau(simplifiedXY, xPlateau, yPlateau)
            
            # Begin a new plateau at the point of departure
            xPlateau = [x]
            yPlateau = [y]
        
        # Otherwise, continue to build the plateau
        else:
            xPlateau.append(x)
            yPlateau.append(y)
    
    # Process the last plateau
    process_plateau(simplifiedXY, xPlateau, yPlateau)
    
    # Return the simplified data as a numpy array
    return np.array(simplifiedXY, dtype=float)

def plotly_per_contig(dotsX, dotsY, simplificationAlg, tolerance,
                      outputDirectory, statisticLabel="Smoothed Q13 depth"):
    '''
    Parameters:
        dotsX -- a dictionary linking chromosome IDs (keys) to lists of integers
                 indicating the position where a dot is located
        dotsY -- a dictionary linking chromosome IDs (keys) to lists of floats
                 indicating the depth value for each dot
        simplificationAlg -- a string indicating the simplification algorithm to use;
                             currently recognises "rdp" for Ramer-Douglas-Peucker, 
                             "lang" for Lang simplification, and "zks" for my simplification
        tolerance -- an integer or float value used by th simplification algorithms to indicate
                     the maximum increase or decrease in y-axis value before a change is detected
                     and the point is retained
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
        xy = np.column_stack((x, y))
        
        # Simplify data for lighter weight plotting
        if simplificationAlg == "rdp":
            xy = rdp(xy)
        elif simplificationAlg == "lang":
            xy = lang_simplification(xy, tolerance, 20) # tolerance=5, look_ahead=20
        elif simplificationAlg == "zks":
            xy = simplify_depth(xy, tolerance)
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
    chromosome. This can be useful for identifying regions of low
    or high coverage in a genome assembly.
    
    Coverage data is decimated/simplified for plotting to ensure that the
    output files are not too large. The user can specify the algorithm used
    for simplification, but note that for large genomes, the 'zks' algorithm
    is recommended for speed and efficiency.
    
    This 'zks' algorithm essentially looks for plateaus in the data where values
    do not change by more than a specified tolerance value. The median value of
    the plateau is then used as the starting and ending point for the plateau,
    hence simplifying all data in between.
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
    p.add_argument("--alg", dest="simplificationAlg",
                   required=False,
                   choices=["rdp", "lang", "zks"],
                   help="""Optionally, specify the simplification algorithm to use;
                   currently recognises 'rdp' for Ramer-Douglas-Peucker,
                   'lang' for Lang simplification, and "zks" for Zachary Kenneth
                   Stewart simplification; default is 'zks'""",
                   default="zks")
    p.add_argument("--tolerance", dest="tolerance",
                   required=False,
                   type=float,
                   help="""Optionally, specify the tolerance value for simplification;
                   generally speaking, this is the maximum increase or decrease in
                   depth value before a change is detected and the point is retained;
                   default is 10""",
                   default=10)
    
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
                      args.tolerance, args.outputDirectory)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
