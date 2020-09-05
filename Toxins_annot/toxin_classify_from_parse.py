#! python3
# toxin_classify_from_parse

# Intended to work after parse_domtblout_segregate.py to classify
# sequences as putative toxins on the basis of their domain matches.
# These domains are part of the toxin annotation project (more details
# forthcoming?).

import argparse, os
from Bio import SeqIO
from Bio.Seq import Seq

# Define functions for later use

## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputParse):
                print('I am unable to locate the parsed HMMER domtblout file (' + args.inputParse + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.inputFasta):
                print('I am unable to locate the FASTA file (' + args.inputFasta + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate that output directory is sensible
        args.outputDirectory = os.path.abspath(args.outputDirectory)
        if not os.path.isdir(args.outputDirectory):
                print('The output directory does not exist (' + args.outputDirectory + ')')
                print('Make sure you\'ve typed the location correctly and try again.')
                quit()
        return args

## Parse input
def read_parsed_domtblout(parsedFile):
        parseDict = {}
        with open(parsedFile, "r") as fileIn:
                for line in fileIn:
                        if line == "" or line =="\r\n" or line == "\n": continue
                        sl = line.split("\t")
                        while sl[-1] == "" or sl[-1] == "\r\n" or sl[-1] == "\n":
                                del sl[-1]
                        for i in range(1, len(sl)):
                                sl[i] = eval(sl[i])
                        parseDict[sl[0]] = sl[1:]
        return parseDict

## Hard-coded rules
def toxin_classifier(sequence, parsedDomainHits): # Receives format of [['Domain', startNum, stopNum, E-value], ...]
        # Acrorhagin
        ACRORHAGIN_DOM_NAME = "1_Acrorhagin_domain"
        ACRORHAGIN_FAMILY_NAME = "Acrorhagin"
        ACRORHAGIN_MIN_EVALUE = 0.01
        ACRORHAGIN_STARTS_WITHIN = 50
        ACRORHAGIN_ENDS_WITHIN = 50
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == ACRORHAGIN_DOM_NAME:
                if float(parsedDomainHits[0][3]) > ACRORHAGIN_MIN_EVALUE:
                        return None
                if int(parsedDomainHits[0][1]) > ACRORHAGIN_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + ACRORHAGIN_ENDS_WITHIN < len(sequence):
                        return None
                return ACRORHAGIN_FAMILY_NAME

        # PLA2
        PLA2_DOM_NAME = "10_PLA2_domain"
        PLA2_FAMILY_NAME = "PLA2"
        PLA2_L2D2_NAME = "50_Z13_Domain"
        PLA2_L2_MIN_EVALUE = 1e-9
        PLA2_L2_STARTS_WITHIN = 80
        PLA2_L2_ENDS_WITHIN = 50
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == PLA2_DOM_NAME and parsedDomainHits[1][0] == PLA2_L2D2_NAME:
                if float(parsedDomainHits[0][3]) > PLA2_L2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > PLA2_L2_MIN_EVALUE:
                        return None
                if int(parsedDomainHits[0][1]) > PLA2_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + PLA2_L2_ENDS_WITHIN < len(sequence):
                        return None
                return PLA2_FAMILY_NAME  + "_L2"
        PLA2_L1_MIN_EVALUE = 1e-4
        PLA2_L1_STARTS_WITHIN = 90
        PLA2_L1_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == PLA2_DOM_NAME:
                if float(parsedDomainHits[0][3]) > PLA2_L1_MIN_EVALUE:
                        return None
                if int(parsedDomainHits[0][1]) > PLA2_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + PLA2_L1_ENDS_WITHIN < len(sequence):
                        return None
                return PLA2_FAMILY_NAME  + "_L1"
        
        # SCRiP
        SCRIP_DOM_NAME = "11_SCRiP_domain"
        SCRIP_FAMILY_NAME = "SCRiP"
        SCRIP_MIN_EVALUE = 0.1
        SCRIP_STARTS_WITHIN = 110
        SCRIP_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == SCRIP_DOM_NAME:
                if float(parsedDomainHits[0][3]) > SCRIP_MIN_EVALUE:
                        return None
                if int(parsedDomainHits[0][1]) > SCRIP_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + SCRIP_ENDS_WITHIN < len(sequence):
                        return None
                return SCRIP_FAMILY_NAME

        # SA8
        SA8_DOM_NAME = "12_Sea_anemone_8_domain"
        SA8_FAMILY_NAME = "Sea_Anemone_8"
        SA8_MIN_EVALUE = 1e-4
        SA8_STARTS_WITHIN = 60
        SA8_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == SA8_DOM_NAME:
                if float(parsedDomainHits[0][3]) > SA8_MIN_EVALUE:
                        return None
                if int(parsedDomainHits[0][1]) > SA8_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + SA8_ENDS_WITHIN < len(sequence):
                        return None
                return SA8_FAMILY_NAME
        
        # PEPTIDASE S1
        S1_DOM_NAME = "14_Peptidase_S1_domain"
        S1_FAMILY_NAME = "Peptidase_S1"
        S1_L2D2_NAME = "13_ShK-like_domain"
        S1_L2D1_MIN_EVALUE = 1e-1
        S1_L2D2_MIN_EVALUE = 1e-60
        S1_L2_STARTS_WITHIN = 50
        S1_L2_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == S1_L2D2_NAME and parsedDomainHits[1][0] == S1_DOM_NAME:
                if float(parsedDomainHits[0][3]) > S1_L2D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > S1_L2D2_MIN_EVALUE:
                        return None
                if int(parsedDomainHits[0][1]) > S1_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + S1_L2_ENDS_WITHIN < len(sequence):
                        return None
                return S1_FAMILY_NAME + "_L2"
        S1_L1_MIN_EVALUE = 1e-50
        S1_L1_STARTS_WITHIN = 90
        S1_L1_ENDS_WITHIN = 80
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == S1_DOM_NAME:
                if float(parsedDomainHits[0][3]) > S1_L1_MIN_EVALUE:
                        return None
                if int(parsedDomainHits[0][1]) > S1_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + S1_L1_ENDS_WITHIN < len(sequence):
                        return None
                return S1_FAMILY_NAME + "_L1"

        # SHK
        SHK_DOM_NAME = "13_ShK-like_domain"
        SHK_FAMILY_NAME = "ShK-like"
        SHK_L3_MIN_EVALUE = 0.1
        SHK_L3_STARTS_WITHIN = 70
        SHK_L3_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 3 and parsedDomainHits[0][0] == SHK_DOM_NAME and parsedDomainHits[1][0] == SHK_DOM_NAME and parsedDomainHits[2][0] == SHK_DOM_NAME:
                if float(parsedDomainHits[0][3]) > SHK_L3_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > SHK_L3_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[2][3]) > SHK_L3_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > SHK_L3_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[2][2]) + SHK_L3_ENDS_WITHIN < len(sequence):
                        return None
                return SHK_FAMILY_NAME + "_L3"
        SHK_L2_MIN_EVALUE = 0.1
        SHK_L2_STARTS_WITHIN = 60
        SHK_L2_ENDS_WITHIN = 50
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == SHK_DOM_NAME and parsedDomainHits[1][0] == SHK_DOM_NAME:
                if float(parsedDomainHits[0][3]) > SHK_L2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > SHK_L2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > SHK_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + SHK_L2_ENDS_WITHIN < len(sequence):
                        return None
                return SHK_FAMILY_NAME + "_L2"
        SHK_L4D2_NAME = "44_Z6_Family"
        SHK_L4D1_MIN_EVALUE = 1e-3
        SHK_L4D2_MIN_EVALUE = 1e-50
        SHK_L4_STARTS_WITHIN = 90
        SHK_L4_ENDS_WITHIN = 15
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == SHK_DOM_NAME and parsedDomainHits[1][0] == SHK_L4D2_NAME:
                if float(parsedDomainHits[0][3]) > SHK_L4D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > SHK_L4D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > SHK_L4_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + SHK_L4_ENDS_WITHIN < len(sequence):
                        return None
                return SHK_FAMILY_NAME + "_L4"
        SHK_L1_MIN_EVALUE = 0.1
        SHK_L1_STARTS_WITHIN = 90
        SHK_L1_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == SHK_DOM_NAME:
                if float(parsedDomainHits[0][3]) > SHK_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > SHK_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + SHK_L1_ENDS_WITHIN < len(sequence):
                        return None
                return SHK_FAMILY_NAME + "_L1"

        # PEPTIDASE M12A
        M12A_DOM_NAME = "15_Peptidase_M12A_domain"
        M12A_FAMILY_NAME = "Peptidase_M12A"
        M12A_L7D2_NAME = "13_ShK-like_domain"
        M12A_L7D3_NAME = "13_ShK-like_domain"
        M12A_L7D1_MIN_EVALUE = 1e-60
        M12A_L7D2_MIN_EVALUE = 1e-1
        M12A_L7_STARTS_WITHIN = 160
        M12A_L7_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 3 and parsedDomainHits[0][0] == M12A_DOM_NAME and parsedDomainHits[1][0] == M12A_L7D2_NAME and parsedDomainHits[2][0] == M12A_L7D3_NAME:
                for i in range(1, len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != M12A_L7D2_NAME:
                                return None
                if float(parsedDomainHits[0][3]) > M12A_L7D1_MIN_EVALUE:
                        return None
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == M12A_L7D2_NAME) and float(parsedDomainHits[i][3]) <= M12A_L7D2_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L7_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + M12A_L7_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L7"
        M12A_L5D2_NAME = "52_Z15_Domain"
        M12A_L5D3_NAME = "49_Z11_Domain"
        M12A_L5D1_MIN_EVALUE = 1e-60
        M12A_L5D2_MIN_EVALUE = 1e-4
        M12A_L5D3_MIN_EVALUE = 1e-25
        M12A_L5_STARTS_WITHIN = 160
        M12A_L5_ENDS_WITHIN = 150
        if len(parsedDomainHits) >= 3 and parsedDomainHits[0][0] == M12A_DOM_NAME and parsedDomainHits[1][0] == M12A_L5D2_NAME and parsedDomainHits[2][0] == M12A_L5D3_NAME:
                for i in range(2, len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != M12A_L5D3_NAME:
                                return None
                if float(parsedDomainHits[0][3]) > M12A_L5D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > M12A_L5D2_MIN_EVALUE:
                        return None
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == M12A_L5D3_NAME) and float(parsedDomainHits[i][3]) <= M12A_L5D3_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L5_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + M12A_L5_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L5"
        M12A_L2D2_NAME = "49_Z12_Domain"
        M12A_L2D1_MIN_EVALUE = 1e-50
        M12A_L2D2_MIN_EVALUE = 1e-17
        M12A_L2_STARTS_WITHIN = 110
        M12A_L2_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == M12A_DOM_NAME and parsedDomainHits[1][0] == M12A_L2D2_NAME:
                if float(parsedDomainHits[0][3]) > M12A_L2D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > M12A_L2D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + M12A_L2_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L2"
        M12A_L3D2_NAME = "49_Z11_Domain"
        M12A_L3D1_MIN_EVALUE = 1e-60
        M12A_L3D2_MIN_EVALUE = 1e-4
        M12A_L3_STARTS_WITHIN = 160
        M12A_L3_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == M12A_DOM_NAME and parsedDomainHits[1][0] == M12A_L3D2_NAME:
                if float(parsedDomainHits[0][3]) > M12A_L3D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > M12A_L3D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L3_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + M12A_L3_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L3"
        M12A_L4D2_NAME = "52_Z15_Domain"
        M12A_L4D1_MIN_EVALUE = 1e-60
        M12A_L4D2_MIN_EVALUE = 1e-4
        M12A_L4_STARTS_WITHIN = 160
        M12A_L4_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == M12A_DOM_NAME and parsedDomainHits[1][0] == M12A_L4D2_NAME:
                if float(parsedDomainHits[0][3]) > M12A_L4D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > M12A_L4D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L4_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + M12A_L4_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L4"
        M12A_L6D2_NAME = "13_ShK-like_domain"
        M12A_L6D1_MIN_EVALUE = 1e-60
        M12A_L6D2_MIN_EVALUE = 1e-1
        M12A_L6_STARTS_WITHIN = 160
        M12A_L6_ENDS_WITHIN = 50
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == M12A_DOM_NAME and parsedDomainHits[1][0] == M12A_L6D2_NAME:
                if float(parsedDomainHits[0][3]) > M12A_L6D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > M12A_L6D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L6_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + M12A_L6_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L6"
        M12A_L8D1_NAME = "56_Z19_Domain"
        M12A_L8D1_MIN_EVALUE = 1e-30
        M12A_L8D2_MIN_EVALUE = 1e-60
        M12A_L8_STARTS_WITHIN = 80
        M12A_L8_ENDS_WITHIN = 120
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == M12A_L8D1_NAME and parsedDomainHits[1][0] == M12A_DOM_NAME:
                if float(parsedDomainHits[0][3]) > M12A_L8D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > M12A_L8D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L8_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + M12A_L8_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L8"
        M12A_L1_MIN_EVALUE = 1e-20
        M12A_L1_STARTS_WITHIN = 200
        M12A_L1_ENDS_WITHIN = 100
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == M12A_DOM_NAME:
                if float(parsedDomainHits[0][3]) > M12A_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > M12A_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + M12A_L1_ENDS_WITHIN < len(sequence):
                        return None
                return M12A_FAMILY_NAME + "_L1"

        # FACTORV
        FV_DOM_NAME = "16_FactorV-like_domain"
        FV_FAMILY_NAME = "FactorV-like"
        FV_L3D1_MIN_EVALUE = 1e-20
        FV_L3D2_MIN_EVALUE = 1e-20
        FV_L3D3_MIN_EVALUE = 0.1
        FV_L3_STARTS_WITHIN = 70
        FV_L3_ENDS_WITHIN = 100
        if len(parsedDomainHits) >= 3 and parsedDomainHits[0][0] == FV_DOM_NAME and parsedDomainHits[1][0] == FV_DOM_NAME and parsedDomainHits[2][0] == FV_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != FV_DOM_NAME:
                                return None
                if float(parsedDomainHits[0][3]) > FV_L3D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > FV_L3D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[2][3]) > FV_L3D3_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > FV_L3_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + FV_L3_ENDS_WITHIN < len(sequence):
                        return None
                return FV_FAMILY_NAME  + "_L3"
        FV_L2D1_MIN_EVALUE = 1e-20
        FV_L2D2_MIN_EVALUE = 0.1
        FV_L2_STARTS_WITHIN = 70
        FV_L2_ENDS_WITHIN = 100
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == FV_DOM_NAME and parsedDomainHits[1][0] == FV_DOM_NAME:
                if float(parsedDomainHits[0][3]) > FV_L2D1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[1][3]) > FV_L2D2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > FV_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + FV_L2_ENDS_WITHIN < len(sequence):
                        return None
                return FV_FAMILY_NAME + "_L2"
        FV_L1_MIN_EVALUE = 1e-20
        FV_L1_STARTS_WITHIN = 130
        FV_L1_ENDS_WITHIN = 250
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == FV_DOM_NAME:
                if float(parsedDomainHits[0][3]) > FV_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > FV_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + FV_L1_ENDS_WITHIN < len(sequence):
                        return None
                return FV_FAMILY_NAME + "_L1"
        
        # ACTINOPORIN
        APORIN_DOM_NAME = "2_Actinoporin_domain"
        APORIN_FAMILY_NAME = "Actinoporin"
        APORIN_MIN_EVALUE = 1e-30
        APORIN_STARTS_WITHIN = 60
        APORIN_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == APORIN_DOM_NAME:
                if float(parsedDomainHits[0][3]) > APORIN_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > APORIN_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + APORIN_ENDS_WITHIN < len(sequence):
                        return None
                return APORIN_FAMILY_NAME
        

        # BBH-LIKE
        BBH_DOM_NAME = "3_BBH-like_domain"
        BBH_FAMILY_NAME = "BBH-like"
        BBH_L3_MIN_EVALUE = 1e-5
        BBH_L3_STARTS_WITHIN = 50
        BBH_L3_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 3 and parsedDomainHits[0][0] == BBH_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != BBH_DOM_NAME:
                                return None
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == BBH_DOM_NAME) and float(parsedDomainHits[i][3]) <= BBH_L3_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > BBH_L3_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + BBH_L3_ENDS_WITHIN < len(sequence):
                        return None
                return BBH_FAMILY_NAME + "_L3"
        BBH_L2_MIN_EVALUE = 1e-4
        BBH_L2_STARTS_WITHIN = 80
        BBH_L2_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == BBH_DOM_NAME and parsedDomainHits[1][0] == BBH_DOM_NAME:
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == BBH_DOM_NAME) and float(parsedDomainHits[i][3]) <= BBH_L2_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > BBH_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + BBH_L2_ENDS_WITHIN < len(sequence):
                        return None
                return BBH_FAMILY_NAME + "_L2"
        BBH_L1_MIN_EVALUE = 1e-5
        BBH_L1_STARTS_WITHIN = 80
        BBH_L1_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == BBH_DOM_NAME:
                if float(parsedDomainHits[0][3]) > BBH_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > BBH_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + BBH_L1_ENDS_WITHIN < len(sequence):
                        return None
                return BBH_FAMILY_NAME + "_L1"
        
        # DEFENSIN
        DEF_DOM_NAME = "4_Defensin-like_domain"
        DEF_FAMILY_NAME = "Defensin"
        DEF_MIN_EVALUE = 1e-3
        DEF_STARTS_WITHIN = 100
        DEF_ENDS_WITHIN = 35
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == DEF_DOM_NAME:
                if float(parsedDomainHits[0][3]) > DEF_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > DEF_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + DEF_ENDS_WITHIN < len(sequence):
                        return None
                return DEF_FAMILY_NAME
        
        # IGFBP-LIKE
        IGFBP_DOM_NAME = "38_IGFBP-like_domain"
        IGFBP_FAMILY_NAME = "IGFBP-like"
        IGFBP_L3_MIN_EVALUE = 1e-7
        IGFBP_L3_STARTS_WITHIN = 120
        IGFBP_L3_ENDS_WITHIN = 700
        if len(parsedDomainHits) >= 3 and parsedDomainHits[0][0] == IGFBP_DOM_NAME and parsedDomainHits[1][0] == IGFBP_DOM_NAME and parsedDomainHits[2][0] == IGFBP_DOM_NAME:
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == IGFBP_DOM_NAME) and float(parsedDomainHits[i][3]) <= IGFBP_L3_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > IGFBP_L3_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + IGFBP_L3_ENDS_WITHIN < len(sequence):
                        return None
                return IGFBP_FAMILY_NAME + "_L3"
        IGFBP_L2_MIN_EVALUE = 1e-15
        IGFBP_L2_STARTS_WITHIN = 300
        IGFBP_L2_ENDS_WITHIN = 600
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == IGFBP_DOM_NAME and parsedDomainHits[1][0] == IGFBP_DOM_NAME:
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == IGFBP_DOM_NAME) and float(parsedDomainHits[i][3]) <= IGFBP_L2_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > IGFBP_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + IGFBP_L2_ENDS_WITHIN < len(sequence):
                        return None
                return IGFBP_FAMILY_NAME + "_L2"
        IGFBP_L1_MIN_EVALUE = 1e-15
        IGFBP_L1_STARTS_WITHIN = 160
        IGFBP_L1_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == IGFBP_DOM_NAME:
                if float(parsedDomainHits[0][3]) > IGFBP_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > IGFBP_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + IGFBP_L1_ENDS_WITHIN < len(sequence):
                        return None
                return IGFBP_FAMILY_NAME + "_L1"

        # ICK-LIKE
        ICK_DOM_NAME = "6_ICK-like_domain"
        ICK_FAMILY_NAME = "ICK-like"
        ICK_MIN_EVALUE = 1e-2
        ICK_STARTS_WITHIN = 110
        ICK_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == ICK_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != ICK_DOM_NAME:
                                return None
                if float(parsedDomainHits[0][3]) > ICK_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > ICK_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + ICK_ENDS_WITHIN < len(sequence):
                        return None
                return ICK_FAMILY_NAME

        # KAZAL-LIKE
        KAZAL_DOM_NAME = "7_Kazal-like_domain_domain"
        KAZAL_FAMILY_NAME = "Kazal-like"
        KAZAL_L2_MIN_EVALUE = 1e-12
        KAZAL_L2_STARTS_WITHIN = 90
        KAZAL_L2_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == KAZAL_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != KAZAL_DOM_NAME:
                                return None
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == KAZAL_DOM_NAME) and float(parsedDomainHits[i][3]) <= KAZAL_L2_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > KAZAL_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + KAZAL_L2_ENDS_WITHIN < len(sequence):
                        return None
                return KAZAL_FAMILY_NAME + "_L2"
        KAZAL_L1_MIN_EVALUE = 1e-12
        KAZAL_L1_STARTS_WITHIN = 90
        KAZAL_L1_ENDS_WITHIN = 30
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == KAZAL_DOM_NAME:
                if float(parsedDomainHits[0][3]) > KAZAL_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > KAZAL_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + KAZAL_L1_ENDS_WITHIN < len(sequence):
                        return None
                return KAZAL_FAMILY_NAME + "_L1"

        # KUNITZ-TYPE
        KUNITZ_DOM_NAME = "8_Kunitz-type_domain"
        KUNITZ_FAMILY_NAME = "Kunitz-type"
        KUNITZ_L3_MIN_EVALUE = 1e-7
        KUNITZ_L3_STARTS_WITHIN = 160
        KUNITZ_L3_ENDS_WITHIN = 40
        if len(parsedDomainHits) >= 3 and parsedDomainHits[0][0] == KUNITZ_DOM_NAME and parsedDomainHits[1][0] == KUNITZ_DOM_NAME and parsedDomainHits[2][0] == KUNITZ_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != KUNITZ_DOM_NAME:
                                return None
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == KUNITZ_DOM_NAME) and float(parsedDomainHits[i][3]) > KUNITZ_L3_MIN_EVALUE:
                                ALL_GOOD = False
                                return None
                if float(parsedDomainHits[0][1]) > KUNITZ_L3_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + KUNITZ_L3_ENDS_WITHIN < len(sequence):
                        return None
                return KUNITZ_FAMILY_NAME + "_L3"
        KUNITZ_L2_MIN_EVALUE = 1e-7
        KUNITZ_L2_STARTS_WITHIN = 160
        KUNITZ_L2_ENDS_WITHIN = 40
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == KUNITZ_DOM_NAME and parsedDomainHits[1][0] == KUNITZ_DOM_NAME:
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == KUNITZ_DOM_NAME) and float(parsedDomainHits[i][3]) <= KUNITZ_L2_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > KUNITZ_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[1][2]) + KUNITZ_L2_ENDS_WITHIN < len(sequence):
                        return None
                return KUNITZ_FAMILY_NAME + "_L2"
        KUNITZ_L1_MIN_EVALUE = 1e-8
        KUNITZ_L1_STARTS_WITHIN = 160
        KUNITZ_L1_ENDS_WITHIN = 120
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == KUNITZ_DOM_NAME:
                if float(parsedDomainHits[0][3]) > KUNITZ_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > KUNITZ_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + KUNITZ_L1_ENDS_WITHIN < len(sequence):
                        return None
                return KUNITZ_FAMILY_NAME + "_L1"

        # CREC
        CREC_DOM_NAME = "37_CREC_domain"
        CREC_FAMILY_NAME = "CREC"
        CREC_MIN_EVALUE = 1e-20
        CREC_STARTS_WITHIN = 70
        CREC_ENDS_WITHIN = 15
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == CREC_DOM_NAME:
                if float(parsedDomainHits[0][3]) > CREC_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > CREC_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + CREC_ENDS_WITHIN < len(sequence):
                        return None
                return CREC_FAMILY_NAME
        
        # EGF
        EGF_DOM_NAME = "5_EGF-like_domain"
        EGF_FAMILY_NAME = "EGF-like"
        EGF_MIN_EVALUE = 1e-8
        EGF_STARTS_WITHIN = 80
        EGF_ENDS_WITHIN = 45
        if len(parsedDomainHits) == 1 and parsedDomainHits[0][0] == EGF_DOM_NAME:
                if float(parsedDomainHits[0][3]) > EGF_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > EGF_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + EGF_ENDS_WITHIN < len(sequence):
                        return None
                return EGF_FAMILY_NAME

        # 17_U1
        U1_DOM_NAME = "17_U1_domain"
        U1_FAMILY_NAME = "U1"
        U1_MIN_EVALUE = 1e-12
        U1_STARTS_WITHIN = 30
        U1_ENDS_WITHIN = 100
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U1_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U1_ENDS_WITHIN < len(sequence):
                        return None
                return U1_FAMILY_NAME
        
        # 18_U2
        U2_DOM_NAME = "18_U2_domain"
        U2_FAMILY_NAME = "U2"
        U2_MIN_EVALUE = 1e-5
        U2_STARTS_WITHIN = 70
        U2_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U2_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U2_ENDS_WITHIN < len(sequence):
                        return None
                return U2_FAMILY_NAME
        
        # 20_U4
        U4_DOM_NAME = "20_U4_domain"
        U4_FAMILY_NAME = "U5"
        U4_MIN_EVALUE = 1e-25
        U4_STARTS_WITHIN = 70
        U4_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U4_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U4_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U4_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U4_ENDS_WITHIN < len(sequence):
                        return None
                return U4_FAMILY_NAME

        # 21_U5
        U5_DOM_NAME = "21_U5_domain"
        U5_FAMILY_NAME = "U5"
        U5_MIN_EVALUE = 1e-25
        U5_STARTS_WITHIN = 70
        U5_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U5_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U5_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U5_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U5_ENDS_WITHIN < len(sequence):
                        return None
                return U5_FAMILY_NAME

        # 23_U7
        U7_DOM_NAME = "23_U7_domain"
        U7_FAMILY_NAME = "U7"
        U7_MIN_EVALUE = 1e-45
        U7_STARTS_WITHIN = 20
        U7_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U7_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U7_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U7_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U7_ENDS_WITHIN < len(sequence):
                        return None
                return U7_FAMILY_NAME

        # 24_U8
        U8_DOM_NAME = "24_U8_domain"
        U8_FAMILY_NAME = "U8"
        U8_MIN_EVALUE = 1e-3
        U8_STARTS_WITHIN = 70
        U8_ENDS_WITHIN = 25
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U8_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U8_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U8_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U8_ENDS_WITHIN < len(sequence):
                        return None
                return U8_FAMILY_NAME

        # 25_U9
        U9_DOM_NAME = "25_U9_domain"
        U9_FAMILY_NAME = "U9"
        U9_MIN_EVALUE = 1e-35
        U9_STARTS_WITHIN = 60
        U9_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U9_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U9_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U9_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U9_ENDS_WITHIN < len(sequence):
                        return None
                return U9_FAMILY_NAME

        # 26_U10
        U10_DOM_NAME = "26_U10_domain"
        U10_FAMILY_NAME = "U10"
        U10_MIN_EVALUE = 1e-35
        U10_STARTS_WITHIN = 70
        U10_ENDS_WITHIN = 15
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U10_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U10_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U10_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U10_ENDS_WITHIN < len(sequence):
                        return None
                return U10_FAMILY_NAME

        # 27_U11
        U11_DOM_NAME = "27_U11_domain"
        U11_FAMILY_NAME = "U11"
        U11_L2_MIN_EVALUE = 1e-10
        U11_L2_STARTS_WITHIN = 45
        U11_L2_ENDS_WITHIN = 15
        if len(parsedDomainHits) >= 2 and parsedDomainHits[0][0] == U11_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != U11_DOM_NAME:
                                return None
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == U11_DOM_NAME) and float(parsedDomainHits[i][3]) <= U11_L2_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > U11_L2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + U11_L2_ENDS_WITHIN < len(sequence):
                        return None
                return U11_FAMILY_NAME + "_L2"
        U11_L1_MIN_EVALUE = 1e-15
        U11_L1_STARTS_WITHIN = 90
        U11_L1_ENDS_WITHIN = 35
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U11_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U11_L1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U11_L1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U11_L1_ENDS_WITHIN < len(sequence):
                        return None
                return U11_FAMILY_NAME + "_L1"

        # 28_U12
        U12_DOM_NAME = "28_U12_domain"
        U12_FAMILY_NAME = "U12"
        U12_MIN_EVALUE = 1e-45
        U12_STARTS_WITHIN = 90
        U12_ENDS_WITHIN = 50
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U12_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U12_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U12_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U12_ENDS_WITHIN < len(sequence):
                        return None
                return U12_FAMILY_NAME

        # 29_U13
        U13_DOM_NAME = "29_U13_domain"
        U13_FAMILY_NAME = "U13"
        U13_MIN_EVALUE = 1e-25
        U13_STARTS_WITHIN = 40
        U13_ENDS_WITHIN = 45
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U13_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U13_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U13_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U13_ENDS_WITHIN < len(sequence):
                        return None
                return U13_FAMILY_NAME

        # 30_U14
        U14_DOM_NAME = "30_U14_domain"
        U14_FAMILY_NAME = "U14"
        U14_MIN_EVALUE = 1e-50
        U14_STARTS_WITHIN = 25
        U14_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U14_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U14_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U14_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U14_ENDS_WITHIN < len(sequence):
                        return None
                return U14_FAMILY_NAME

        # 31_U15
        U15_DOM_NAME = "31_U15_domain"
        U15_FAMILY_NAME = "U15"
        U15_MIN_EVALUE = 1e-15
        U15_STARTS_WITHIN = 190
        U15_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U15_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != U15_DOM_NAME:
                                return None
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == U15_DOM_NAME) and float(parsedDomainHits[i][3]) <= U15_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > U15_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + U15_ENDS_WITHIN < len(sequence):
                        return None
                return U15_FAMILY_NAME

        # 32_U16
        U16_DOM_NAME = "32_U16_domain"
        U16_FAMILY_NAME = "U16"
        U16_MIN_EVALUE = 1e-10
        U16_STARTS_WITHIN = 120
        U16_ENDS_WITHIN = 50
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U16_DOM_NAME:
                for i in range(len(parsedDomainHits)):
                        if parsedDomainHits[i][0] != U16_DOM_NAME:
                                return None
                ONE_GOOD = False
                for i in range(len(parsedDomainHits)):
                        if (parsedDomainHits[i][0] == U16_DOM_NAME) and float(parsedDomainHits[i][3]) <= U16_MIN_EVALUE:
                                ONE_GOOD = True
                                break
                if ONE_GOOD == False:
                        return None
                if float(parsedDomainHits[0][1]) > U16_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[-1][2]) + U16_ENDS_WITHIN < len(sequence):
                        return None
                return U16_FAMILY_NAME

        # 35_U19
        U19_DOM_NAME = "35_U19_domain"
        U19_FAMILY_NAME = "U19"
        U19_MIN_EVALUE = 1e-60
        U19_STARTS_WITHIN = 180
        U19_ENDS_WITHIN = 40
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U19_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U19_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U19_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U19_ENDS_WITHIN < len(sequence):
                        return None
                return U19_FAMILY_NAME

        # 36_U20
        U20_DOM_NAME = "36_U20_domain"
        U20_FAMILY_NAME = "U20"
        U20_MIN_EVALUE = 1e-100
        U20_STARTS_WITHIN = 40
        U20_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == U20_DOM_NAME:
                if float(parsedDomainHits[0][3]) > U20_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > U20_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + U20_ENDS_WITHIN < len(sequence):
                        return None
                return U20_FAMILY_NAME

        # 39_Z1
        Z1_DOM_NAME = "39_Z1_Family"
        Z1_FAMILY_NAME = "Z1"
        Z1_MIN_EVALUE = 1e-45
        Z1_STARTS_WITHIN = 60
        Z1_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z1_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z1_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z1_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z1_ENDS_WITHIN < len(sequence):
                        return None
                return Z1_FAMILY_NAME

        # 40_Z2
        Z2_DOM_NAME = "40_Z2_Family"
        Z2_FAMILY_NAME = "Z2"
        Z2_MIN_EVALUE = 1e-70
        Z2_STARTS_WITHIN = 110
        Z2_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z2_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z2_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z2_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z2_ENDS_WITHIN < len(sequence):
                        return None
                return Z2_FAMILY_NAME

        # 41_Z3
        Z3_DOM_NAME = "41_Z3_Family"
        Z3_FAMILY_NAME = "Z3"
        Z3_MIN_EVALUE = 1e-30
        Z3_STARTS_WITHIN = 120
        Z3_ENDS_WITHIN = 40
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z3_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z3_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z3_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z3_ENDS_WITHIN < len(sequence):
                        return None
                return Z3_FAMILY_NAME

        # 42_Z4
        Z4_DOM_NAME = "42_Z4_Family"
        Z4_FAMILY_NAME = "Z4"
        Z4_MIN_EVALUE = 1e-90
        Z4_STARTS_WITHIN = 35
        Z4_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z4_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z4_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z4_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z4_ENDS_WITHIN < len(sequence):
                        return None
                return Z4_FAMILY_NAME

        # 43_Z5
        Z5_DOM_NAME = "43_Z5_Family"
        Z5_FAMILY_NAME = "Z5"
        Z5_MIN_EVALUE = 1e-30
        Z5_STARTS_WITHIN = 80
        Z5_ENDS_WITHIN = 90
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z5_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z5_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z5_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z5_ENDS_WITHIN < len(sequence):
                        return None
                return Z5_FAMILY_NAME

        # 44_Z6
        Z6_DOM_NAME = "44_Z6_Family"
        Z6_FAMILY_NAME = "Z6"
        Z6_MIN_EVALUE = 1e-40
        Z6_STARTS_WITHIN = 80
        Z6_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z6_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z6_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z6_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z6_ENDS_WITHIN < len(sequence):
                        return None
                return Z6_FAMILY_NAME

        # 45_Z7
        Z7_DOM_NAME = "45_Z7_Family"
        Z7_FAMILY_NAME = "Z7"
        Z7_MIN_EVALUE = 1e-140
        Z7_STARTS_WITHIN = 60
        Z7_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z7_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z7_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z7_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z7_ENDS_WITHIN < len(sequence):
                        return None
                return Z7_FAMILY_NAME

        # 46_Z8
        Z8_DOM_NAME = "46_Z8_Family"
        Z8_FAMILY_NAME = "Z8"
        Z8_MIN_EVALUE = 1e-40
        Z8_STARTS_WITHIN = 50
        Z8_ENDS_WITHIN = 20
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z8_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z8_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z8_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z8_ENDS_WITHIN < len(sequence):
                        return None
                return Z8_FAMILY_NAME

        # 47_Z9
        Z9_DOM_NAME = "47_Z9_Family"
        Z9_FAMILY_NAME = "Z9"
        Z9_MIN_EVALUE = 1e-65
        Z9_STARTS_WITHIN = 70
        Z9_ENDS_WITHIN = 10
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z9_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z9_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z9_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z9_ENDS_WITHIN < len(sequence):
                        return None
                return Z9_FAMILY_NAME

        # 48_Z10
        Z10_DOM_NAME = "48_Z10_Family"
        Z10_FAMILY_NAME = "Z10"
        Z10_MIN_EVALUE = 1e-120
        Z10_STARTS_WITHIN = 60
        Z10_ENDS_WITHIN = 400
        if len(parsedDomainHits) >= 1 and parsedDomainHits[0][0] == Z10_DOM_NAME:
                if float(parsedDomainHits[0][3]) > Z10_MIN_EVALUE:
                        return None
                if float(parsedDomainHits[0][1]) > Z10_STARTS_WITHIN:
                        return None
                if int(parsedDomainHits[0][2]) + Z10_ENDS_WITHIN < len(sequence):
                        return None
                return Z10_FAMILY_NAME

## Classify toxins with hard-coded rules
def classify_all_toxins(fastaDict, parsedDomtbloutDict):
        toxinDict = {}
        for key, value in parsedDomtbloutDict.items():
                if key in fastaDict:
                        seqid = key
                else:
                        seqid = effort_find_unique_seqid_matches(key, fastaDict)
                        if seqid == None:
                                print("This fasta file is jank AF. We can't use it, sorry.")
                                print("(sequence ID '{0}' wasn't found in the fasta file)".format(key))
                                return None
                toxinClass = toxin_classifier(str(fastaDict[seqid].seq), value)
                toxinDict[key] = toxinClass
        return toxinDict

def debug_helper(key, fastaDict, parsedDomtbloutDict):
        if key in fastaDict:
                seqid = key
        else:
                seqid = effort_find_unique_seqid_matches(key, fastaDict)
        return str(fastaDict[seqid].seq), parsedDomtbloutDict[key]

## Output function
def filename_generator(prefix, suffix):
        buffer = ""
        while os.path.exists(prefix + str(buffer) + suffix):
                if buffer == "":
                        buffer = 1
                else:
                        buffer += 1
        return prefix + str(buffer) + suffix

def output_func(toxinDict, fastaDict, outputDirectory):
        groupsDict = {}
        for key, value in toxinDict.items():
                if value == None: continue
                if value not in groupsDict:
                        groupsDict[value] = [key]
                else:
                        groupsDict[value].append(key)
        for groupKey, groupValue in groupsDict.items():
                fileName = filename_generator(os.path.join(outputDirectory, groupKey), ".fasta")
                with open(fileName, "w") as fileOut:
                        for value in groupValue:
                                seq = str(fastaDict[value].seq)
                                fileOut.write(">{0}\n{1}\n".format(value, seq))

# Fasta parser
def fasta_to_dict(inputFasta):
        try:
                records = SeqIO.to_dict(SeqIO.parse(open(inputFasta, "r"), "fasta"))
        except:
                records = {}
                with open(inputFasta, "r") as fileIn:
                        for line in fileIn:
                                if line.startswith(">"):
                                        seqid = line.strip(">\r\n")
                                        records[seqid] = SeqIO.SeqRecord(Seq(""), id=seqid)
                                        prevSeqid = seqid
                                else:
                                        records[prevSeqid].seq += line.strip("\r\n ")
        return records

def effort_find_unique_seqid_matches(shortID, fastaDict):
        idMatches = [seqid for seqid in fastaDict.keys() if seqid.startswith(shortID)]
        if idMatches == [] or len(idMatches) > 1:
                return None
        else:
                return idMatches[0]

#### USER INPUT SECTION
usage = """%(prog)s reads a parsed HMMER domtblout file, as created by parse_domtblout_segregate.py,
and classifies sequences into putative toxin families on the basis of their domain architecture.
Certain hard-coded rules are enforced to remove sequences that do not conform to previous identified
family structures.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "--input", dest="inputParse",
                   help="Input parsed HMMER3 domtblout file")
p.add_argument("-f", "--fasta", dest="inputFasta",
                   help="Input FASTA file of queried sequences")
p.add_argument("-o", "--ouputDir", dest="outputDirectory", default = "",
                   help="Optionally specify the directory to write files to")

args = p.parse_args()
args = validate_args(args)

# Read parsed HMMER file
parseDict = read_parsed_domtblout(args.inputParse)

# Obtain FASTA dictionary
fastaDict = fasta_to_dict(args.inputFasta)

# Classify toxins with hard-coded rules
toxinDict = classify_all_toxins(fastaDict, parseDict)

# Generate output
output_func(toxinDict, fastaDict, args.outputDirectory)

# All done!
print('Program completed successfully!')
