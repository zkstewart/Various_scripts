#! python3
# concat_pngs.py
# Join the replication psQTL analysis results into a single image for each parameter set

import os
from PIL import Image

rawDir = "results"
outDir = "concat_results"

# Locate files that need concatenation
fileGroups = {}
for location in os.listdir():
    if os.path.isdir(location) and "results" in os.listdir(location):
        imageDir = os.path.join(location, "results")
        for imageFile in [ os.path.join(imageDir, f) for f in os.listdir(imageDir) if f.endswith(".png") ]:
            params = os.path.basename(imageFile).strip(".png")
            fileGroups.setdefault(params, [])
            fileGroups[params].append(imageFile)

# Concatenate files within each group
## See https://stackoverflow.com/questions/30227466/combine-several-images-horizontally-with-python
for params, files in fileGroups.items():
    images = [Image.open(file) for file in files]
    widths, heights = zip(*(i.size for i in images))
    
    total_width = max(widths) * 5
    max_height = max(heights) * 2
    
    new_im = Image.new('RGB', (total_width, max_height))
    
    x_offset = 0
    y_offset = 0
    ongoingCount = 0
    for im in images:
        new_im.paste(im, (x_offset, y_offset))
        x_offset += im.size[0]
        
        ongoingCount += 1
        if ongoingCount == 5:
            x_offset = 0
            y_offset += im.size[1]
    
    new_im.save(os.path.join(outDir, f"{params}.png"))

print("Program completed successfully!")
