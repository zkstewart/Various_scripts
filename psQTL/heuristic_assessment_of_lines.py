#! python3
# concat_pngs.py

import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

outFile = "qtl_test.tsv"

def func(x, a, b, c):
    return a*x**2 + b*x + c # quadratic function

####

# Locate files containing plot data
fileGroups = {}
for root, dirs, files in os.walk(os.getcwd()):
    if len(files) == 1 and files[0] == "chr1.0-10000000.call_line.tsv":
        callLineFile = os.path.join(root, files[0])
        params = os.path.dirname(root).split("/", maxsplit=7)[-1].replace("/", "_")
        fileGroups.setdefault(params, [])
        fileGroups[params].append(callLineFile)

# Assess each parameter set
results = [["pop_size", "pop_balance", "phenotype_error", "centre_exceeds_edges", "centre_is_peak",
           "centre_is_only_peak", "centre_exceeds_cutoff", "r_squared"]]
for params, files in fileGroups.items():
    thisResult = [[], [], [], [], []] # check1, check2, check3, check4, r_squared
    for file in files:
        # Parse data from call line file
        x, y = [], []
        with open(file, "r") as fileIn:
            for line in fileIn:
                if not line.startswith("contigID"):
                    contigID, position, ed, smoothedEd = line.strip().split("\t")
                    position = int(position)
                    smoothedEd = float(smoothedEd)
                    
                    x.append(position)
                    y.append(smoothedEd)
        
        # Run heuristic assessments
        # Check 1: Centre exceeds the edges
        CHECK1_THRESHOLD = 2
        centrePosition = len(x) // 2
        centreY = y[centrePosition]
        if centreY > (y[0] * CHECK1_THRESHOLD) and centreY > (y[-1] * CHECK1_THRESHOLD):
            check1 = True
        else:
            check1 = False
        
        # Check 2: Centre is a peak
        CHECK2_RADIUS = 50 # 50 Kbp
        peakY = max(y[centrePosition - CHECK2_RADIUS:centrePosition + CHECK2_RADIUS])
        if peakY == max(y):
            check2 = True
        else:
            check2 = False
        
        # Check 3: Centre is the only peak
        CHECK3_RADIUS = 500000 # 500 Kbp
        maxPositions = [i*1000 for i, j in enumerate(y) if j == peakY]
        if any([ abs(i - 5000000) > CHECK3_RADIUS for i in maxPositions ]):
            check3 = False
        else:
            check3 = True
        
        # Check 4: Centre ED exceeds a cutoff
        CHECK4_THRESHOLD = 0.5
        if centreY > CHECK4_THRESHOLD:
            check4 = True
        else:
            check4 = False
        
        # Score 1: Fit a curve to the data and check the R^2 value
        x = np.array(x)
        y = np.array(y)
        
        # Establish weights to force the curve to go down to 0 at the start and end
        y[[0, -1]] = 0 # set first and last points artificially to 0
        sigma = np.ones(len(x))
        sigma[[0, -1]] = 0.0001 # weight first and last points heavily
        
        # Fit the curve
        popt, pcov = curve_fit(func, x, y,
                               p0=(-0.01, 0, -0.01), # initial a, b, c values
                               sigma=sigma, # ensure line goes through y=0 at start and end
                               bounds=([-1, -1, -np.inf], [0, 1, np.inf])) # make sure a is negative
        
        # Calculate R^2 value
        residuals = y - func(x, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        # plt.plot(xdata, func(xdata, *popt), 'g--',
        #          label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
        # plt.scatter(x, y, label='data')
        # plt.xlabel('x')
        # plt.ylabel('y')
        # plt.title('Curve Fit')
        # plt.legend()
        # plt.show()
        
        # Store results
        thisResult[0].append(check1)
        thisResult[1].append(check2)
        thisResult[2].append(check3)
        thisResult[3].append(check4)
        thisResult[4].append(r_squared)
    
    # Summarise results for this parameter set
    thisResult[0] = np.mean(thisResult[0])
    thisResult[1] = np.mean(thisResult[1])
    thisResult[2] = np.mean(thisResult[2])
    thisResult[3] = np.mean(thisResult[3])
    thisResult[4] = np.median(thisResult[4])
    results.append(params.split("_") + thisResult) # TO-DO: Split params into individual columns

# Write results to file
with open(outFile, "w") as fileOut:
    for result in results:
        result = "\t".join([str(i) for i in result])
        fileOut.write(result + "\n")

print("Program completed successfully!")
