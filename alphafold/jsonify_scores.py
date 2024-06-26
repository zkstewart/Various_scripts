"""
CODE OBTAINED FROM https://subtiwiki.uni-goettingen.de/v4/paeViewerDemo
MODIFIED BY ZACHARY K. STEWART 31-05-2024

ORIGINAL DOCSTRING FOLLOWS:
usage: jsonify_scores.py [-h] [-o OUTPUT_PATH] input_path

Unpickles the output of an AlphaFold-Multimer run and converts it into a
`.json` file which can be used by PAE Viewer.

Warning: The pickle module is not secure. Only unpickle data you trust. It is
possible to construct malicious pickle data which will execute arbitrary code
during unpickling. Never unpickle data that could have come from an untrusted
source, or that could have been tampered with.

positional arguments:
  input_path            Path to the `.pickle` file containing scores as
                        generated by AlphaFold.

options:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to the output JSON file. Will fail if file
                        already exists. Saves the output in the same directory
                        as the input file by default.

OBVIOUSLY, I HAVE CHANGED THIS A LOT. THE ONLY PART THAT REMAINS THE SAME IS THE
TRY:EXCEPT LOGIC AND JSON DUMPING AT THE END.
"""
import os, argparse, json, pickle, gzip
from contextlib import contextmanager

@contextmanager
def open_gz_pickle(filename):
    if str(filename).endswith(".gz"):
        with gzip.open(filename, "rb") as f:
            yield f
    else:
        with open(filename, "rb") as f:
            yield f

def validate_args(args):
    if not os.path.isdir(args.modelPath):
        print(f'I am unable to locate the model outputs directory ({args.modelPath})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()

if __name__ == '__main__':
    usage = """%(prog)s takes a directory as output by AlphaPulldown containing
    ranked_0.pdb, ranked_1.pdb, ... the ranking_debug.json file, and associated
    model pickle files (gzip'd or otherwise). It will produce json files containing
    PAE scores for each model, which can be used by the PAE Viewer at
    https://subtiwiki.uni-goettingen.de/v4/paeViewerDemo for example.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="modelPath",
                   required=True,
                   help="Specify the location of the multimer modelling outputs")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate input files
    pdbFiles = [
        os.path.join(args.modelPath, f)
        for f in os.listdir(args.modelPath)
        if f.startswith("ranked_") and f.endswith(".pdb")
    ]
    
    # Parse ranking json
    rankingFile = os.path.join(args.modelPath, "ranking_debug.json")
    if not os.path.isfile(rankingFile):
        raise FileNotFoundError(f"Unable to locate the ranking_debug.json file in '{args.modelPath}'!")
    rankings = json.load(open(rankingFile))
    
    # Locate model pickle files
    pklFiles = [
        os.path.join(args.modelPath, f)
        for f in os.listdir(args.modelPath)
        if f.startswith("result_model_") and (f.endswith(".pkl") or f.endswith(".pkl.gz"))
    ]
    
    # Validate that all expected files were located
    if not ( len(pdbFiles) == len(rankings["order"]) == len(pklFiles) ):
        raise ValueError((
            f"Expected to find an equal number of .pdb, pickle, and ranked models in the .json file, "
            f"but found {len(pdbFiles)}, {len(pklFiles)}, and {len(rankings['order'])} respectively."
        ))
    for ranking in rankings["order"]:
        if not any([ranking in f for f in pklFiles]):
            raise ValueError(f"Ranking '{ranking}' not found in the .pdb files.")
    
    # Jsonify the scores
    for pklFile in pklFiles:
        # Figure out which ranked model this is
        for modelIndex, modelID in enumerate(rankings["order"]):
            if modelID in pklFile:
                break
        
        # Format output file name based on model rank
        outputJsonName = os.path.join(args.modelPath, f"ranked_{modelIndex}.json")
        if os.path.exists(outputJsonName):
            print(f"Output file '{outputJsonName}' already exists. Skipping...")
            continue
        
        # Open the pickle file and write the json
        with open_gz_pickle(pklFile) as input_file, open(outputJsonName, "w") as output_file:
            scores = pickle.load(input_file)
            
            output_scores = {
                'pae': scores['predicted_aligned_error'].tolist()
            }
            
            try:
                output_scores['max_pae'] = scores[
                    'max_predicted_aligned_error'
                ].item()
            except KeyError:
                pass
            
            try:
                output_scores['plddt'] = scores['plddt'].tolist()
            except KeyError:
                pass
            
            try:
                output_scores['ptm'] = scores['ptm'].item()
            except KeyError:
                pass
            
            try:
                output_scores['iptm'] = scores['iptm'].item()
            except KeyError:
                pass
            
            json.dump(output_scores, output_file, ensure_ascii=False)
    
    print("Program completed successfully!")
