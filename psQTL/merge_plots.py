#! python3
# concat_pngs.py
# Join the replication psQTL analysis results into a single image for each parameter set

import os
import matplotlib.pyplot as plt

outDir = "merge_plots"

# Locate files for joined plotting
fileGroups = {}
for root, dirs, files in os.walk(os.getcwd()):
    if len(files) == 1 and files[0] == "chr1.0-10000000.call_line.tsv":
        callLineFile = os.path.join(root, files[0])
        params = os.path.dirname(root).split("/", maxsplit=7)[-1].replace("/", "_")
        fileGroups.setdefault(params, [])
        fileGroups[params].append(callLineFile)

# Plot each parameter set
for params, files in fileGroups.items():
    fig, ax = plt.subplots()
    for file in files:
        # Parse data from call line file
        x, y = [], []
        lastPosition = None
        lastEd = None
        with open(file, "r") as fileIn:
            for line in fileIn:
                if not line.startswith("contigID"):
                    contigID, position, ed, smoothedEd = line.strip().split("\t")
                    position = int(position)
                    smoothedEd = float(smoothedEd)
                    
                    if lastPosition == None:
                        x.append(position)
                        y.append(smoothedEd)
                    elif smoothedEd != lastEd:
                        x.append(lastPosition)
                        y.append(lastEd)
                        x.append(position)
                        y.append(smoothedEd)
                    lastPosition = position
                    lastEd = smoothedEd
        x.append(lastPosition)
        y.append(lastEd)
        
        # Plot data
        ax.plot(x, y)
    
    ax.set_xlabel("Position")
    ax.set_ylabel("$ED^4$")
    ax.set_title(params)
    fig.savefig(os.path.join(outDir, f"{params}.png"))
    plt.close(fig)

print("Program completed successfully!")
