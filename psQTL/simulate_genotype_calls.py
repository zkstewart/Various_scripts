#! python3
# simulate_genotype_calls.py
# Simulate genotype calls for a population that segregates over a single QTL SNP.
# Used in evaluation of psQTL and the amount of samples needed to detect a QTL.

import os, random, argparse
import numpy as np
import pandas as pd

from itertools import product
from math import sqrt, ceil

from chromax import Simulator
from chromax.sample_data import genetic_map, genome

def random_sample(arr: np.array, size: int = 1) -> np.array:
    return arr[np.random.choice(len(arr), size=size, replace=False)]

def main():
    usage = """%(prog)s simulates genotype calls for a population that segregates over
    a single QTL SNP. It assumes a simple scenario with diploid heterozygous parents
    and a recessive trait. It will output a VCF and metadata file suitable for
    psQTL analysis.
    """
    
    # Parse command line arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-s", dest="seed",
                    required=True,
                    type=int,
                    help="Specify the seed for the simulation")
    p.add_argument("-o", dest="outputDirectory",
                    required=True,
                    help="Specify the location to write outputs")
    p.add_argument("--snpMB", dest="snpMB",
                    required=False,
                    type=int,
                    help="Specify the number of SNPs per megabase",
                    default=1000)
    p.add_argument("--genomeLength", dest="genomeLength",
                    required=False,
                    type=int,
                    help="Specify the genome length in base pairs",
                    default=10000000)
    p.add_argument("--centimorgans", dest="cmMB",
                    required=False,
                    type=float,
                    help="Specify the centiMorgan per megabase",
                    default=3.0)
    p.add_argument("--offspring", dest="numOffspring",
                    required=False,
                    type=int,
                    help="Specify the number of offspring to simulate",
                    default=10000)
    args = p.parse_args()
    
    ########################
    # SIMULATE POPULATIONS #
    ########################
    
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    qtlSNP = args.genomeLength//args.snpMB//2 # the QTL SNP is in the middle of the genome
    if qtlSNP % 1 != 0:
        raise ValueError(f"The QTL SNP index ({qtlSNP}) must be a whole number. Please check the genome length and SNP MB values " +
                        "to make sure they divide evenly.")
    
    # Generate diploid parent
    ploidy1 = np.ones(args.genomeLength//args.snpMB).reshape(1, args.genomeLength//args.snpMB) # ones indicate the positive trait; must be homozygous
    ploidy2 = np.zeros(args.genomeLength//args.snpMB).reshape(1, args.genomeLength//args.snpMB) # zeros indicate the negative trait; any presence means no good trait
    parentGenome = np.dstack((ploidy1, ploidy2))
    parentGenome = parentGenome.astype(bool)
    
    # Duplicate the parent to enable crossing
    parentGenomes = np.vstack((parentGenome, parentGenome))
    
    # Generate genetic map
    mapArray = [["CHR.PHYS", "cM", "Trait"]]
    for i in range(0, args.genomeLength//args.snpMB):
        physicalPosition = i*args.snpMB
        cMPosition = (physicalPosition / 1000000) * args.cmMB
        mapArray.append(["chr1", cMPosition, 0.01])
    geneticMap = pd.DataFrame(mapArray[1:], columns=mapArray[0])
    
    # Init the simulator
    simulator = Simulator(genetic_map=geneticMap, seed=args.seed)
    
    # Load the parent genome
    npyFile = os.path.join(args.outputDirectory, "parentGenomes.npy")
    np.save(npyFile, parentGenomes)
    f0 = simulator.load_population(npyFile)
    
    # Generate the F1 population
    f1, _ = simulator.random_crosses(f0, 1, n_offspring=args.numOffspring) # returns (n_crosses, n_individuals, n_loci, n_alleles)
    f1 = f1.reshape(args.numOffspring, args.genomeLength//args.snpMB, 2) # reshape to (n_individuals, n_loci, n_alleles)
    
    # Separate out progeny based on qtlSNP occurrence
    goodF1 = f1[(f1[:, qtlSNP, 0] == 1) & (f1[:, qtlSNP, 1] == 1)]
    badF1 = f1[(f1[:, qtlSNP, 0] == 0) | (f1[:, qtlSNP, 1] == 0)]
    
    # Generate a dummy genome file
    # MULTILINE_LEN = 80
    # with open("genome.fasta", "w") as fileOut:
    #     fileOut.write(">chr1\n")
    #     for i in range(0, args.genomeLength, MULTILINE_LEN):
    #         fileOut.write("A" * MULTILINE_LEN + "\n")
    
    # Symlink to the genome file instead
    os.symlink(os.path.abspath("genome.fasta"), os.path.join(args.outputDirectory, "genome.fasta"))
    
    ########################
    #    SIMULATE BULKS    #
    ########################
    
    # Set combination variables
    populationSize = range(10, 210, 10)
    populationBalance = [ x/100 for x in range(10, 60, 10) ]
    phenotypeErrorPct = [ x/100 for x in range(0, 55, 5) ]
    
    # Generate all combinations of variables
    productList = list(product(populationSize, populationBalance, phenotypeErrorPct))
    
    # Generate genotype calls for each QTL 'zone'
    numThings = 0
    for size, balance, phePct in productList:
        # Determine bulk1 (good) size
        numBulk1 = size * balance
        
        # Skip if the bulk size is not a whole number
        if numBulk1 % 1 != 0:
            continue
        
        # Convert to integer & derive the bulk2 (bad) size
        numBulk1 = int(numBulk1)
        numBulk2 = size - numBulk1
        
        # Skip if we can't simulate phenotype error with whole numbers
        if (numBulk1 * phePct % 1 != 0) or (numBulk2 * phePct % 1 != 0):
            continue
        
        # Randomly select individuals for each bulk
        bulk1 = random_sample(goodF1, numBulk1)
        bulk2 = random_sample(badF1, numBulk2)
        
        # Generate phenotype error
        numBulk1Error = int(numBulk1 * phePct)
        numBulk2Error = int(numBulk2 * phePct)
        
        # Introduce phenotype error
        if numBulk1Error > 0:
            index = numBulk1 - numBulk1Error
            bulk1, bulk1Error = bulk1[:index], bulk1[index:]
        if numBulk2Error > 0:
            index = numBulk2 - numBulk2Error
            bulk2, bulk2Error = bulk2[:index], bulk2[index:]
        
        if numBulk2Error > 0:
            bulk1 = np.concatenate((bulk1, bulk2Error))
        if numBulk1Error > 0:
            bulk2 = np.concatenate((bulk2, bulk1Error))
        
        # Format metadata
        metadata = []
        b1Samples, b2Samples = [], []
        ongoingCount = 1
        for _ in range(numBulk1):
            metadata.append(f"sample{ongoingCount}\tbulk1")
            b1Samples.append(f"sample{ongoingCount}")
            ongoingCount += 1
        for _ in range(numBulk2):
            metadata.append(f"sample{ongoingCount}\tbulk2")
            b2Samples.append(f"sample{ongoingCount}")
            ongoingCount += 1
        
        # Reformat into VCF
        vcfDict = {
            "#CHROM": ["chr1"] * (args.genomeLength//args.snpMB),
            "POS": [str(x*args.snpMB) for x in range(0, args.genomeLength//args.snpMB)], # physical position
            "ID": ["."] * (args.genomeLength//args.snpMB),
            "REF": ["A"] * (args.genomeLength//args.snpMB),
            "ALT": ["G"] * (args.genomeLength//args.snpMB),
            "QUAL": ["."] * (args.genomeLength//args.snpMB),
            "FILTER": ["PASS"] * (args.genomeLength//args.snpMB),
            "INFO": ["."] * (args.genomeLength//args.snpMB),
            "FORMAT": ["GT"] * (args.genomeLength//args.snpMB)
        }
        for sampleIDs, bulkArray in zip([b1Samples, b2Samples], [bulk1, bulk2]):
            for sampleID, sampleArray in zip(sampleIDs, bulkArray):
                vcfDict[sampleID] = []
                for genotype in sampleArray:
                    gt = list(genotype)
                    gt.sort()
                    vcfDict[sampleID].append(f"{gt[0]}/{gt[1]}")
        vcf = pd.DataFrame(vcfDict)
        
        # Generate files for psQTL analysis
        os.makedirs(args.outputDirectory, exist_ok=True)
        
        popOutDir = os.path.join(args.outputDirectory, str(size))
        os.makedirs(popOutDir, exist_ok=True)
        
        balanceOutDir = os.path.join(popOutDir, str(balance))
        os.makedirs(balanceOutDir, exist_ok=True)
        
        finalOutDir = os.path.join(balanceOutDir, str(phePct))
        os.makedirs(finalOutDir, exist_ok=True)
        
        metadataFile = os.path.join(finalOutDir, "metadata.txt")
        with open(metadataFile, "w") as fileOut:
            fileOut.write("\n".join(metadata))
        
        vcfFile = os.path.join(finalOutDir, "variants.vcf")
        vcf.to_csv(vcfFile, sep="\t", index=False)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
