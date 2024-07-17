'''
This script was created by Zhang Yujian on January 22nd, 2024.
This script is aimed to convert different sequence infomation involved in bioinformation.
'''
from collections import OrderedDict

class ConvKit:
    # codon
    dna2ProDict: dict[str, str] = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', \
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', \
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', \
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', \
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', \
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', \
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', \
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', \
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', \
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', \
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', \
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', \
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', \
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', \
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', \
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }


    rna2ProDict: dict[str, str] = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', \
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', \
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', \
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', \
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', \
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', \
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', \
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', \
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', \
        'ACU': 'U', 'ACC': 'U', 'ACA': 'U', 'ACG': 'U', \
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', \
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', \
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', \
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', \
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', \
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }


    # T/U
    dna2RnaDict: dict[str, str] = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'U'}

    rna2DnaDict: dict[str ,str] = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'T'}

    # DNA Substitution Matrix
    # Unitary matrix
    dnaUnitaryMatrix: dict[bool, int] = {True: 1, False: 0}

    # Transition-transversion matrix
    dnaTTMatrix: dict[str, int] = {True: 1, 'transition': -1, 'transversion': -5}

    # BLAST matrix
    dnaBlastMatrix: dict[bool, int] = {True: 5, False: -4}

    # Protein Substitution Matrix
    # Unitary matrix
    protUnitaryMatrix: dict[bool, int] = {True: 1, False: 0}

    # PAM-250 matrix
    pam250Matrix: dict[str, int] = {
        'CC': 12, \
        'CS': 0, 'CT': -2, 'CP': -3, 'CA': -2, 'CG': -3, \
        'SC': 0, 'TC': -2, 'PC': -3, 'AC': -2, 'GC': -3, \
        'CN': -4, 'CD': -5, 'CE': -5, 'CQ': -5, \
        'NC': -4, 'DC': -5, 'EC': -5, 'QC': -5, \
        'CH': -3, 'CR': -4, 'CK': -5, \
        'HC': -3, 'RC': -4, 'KC': -5, \
        'CM': -5, 'CI': -2, 'CL': -6, 'CV': -2, \
        'MC': -5, 'IC': -2, 'LC': -6, 'VC': -2, \
        'CF': -4, 'CY': 0, 'CW': -8, \
        'FC': -4, 'YC': 0, 'WC': -8, \
        'SS': 2, 'ST': 1, 'SP': 1, 'SA': 1, 'SG': 1, \
        'TS': 1, 'PS': 1, 'AS': 1, 'GS': 1, \
        'SN': 1, 'SD': 0, 'SE': 0, 'SQ': -1, \
        'NS': 1, 'DS': 0, 'ES': 0, 'QS': -1, \
        'SH': -1, 'SR': 0, 'SK': 0, \
        'HS': -1, 'RS': 0, 'KS': 0, \
        'SM': -2, 'SI': -1, 'SL': -3, 'SV': -1, \
        'MS': -2, 'IS': -1, 'LS': -3, 'VS': -1, \
        'SF': -3, 'SY': -3, 'SW': -2, \
        'FS': -3, 'YS': -3, 'WS': -2, \
        'TT': 3, 'TP': 0, 'TA': 1, 'TG': 0, \
        'PT': 0, 'AT': 1, 'GT': 0, \
        'TN': 0, 'TD': 0, 'TE': 0, 'TQ': -1, \
        'NT': 0, 'DT': 0, 'ET': 0, 'QT': -1, \
        'TH': -1, 'TR': -1, 'TK': 0, \
        'HT': -1, 'RT': -1, 'KT': 0, \
        'TM': -1, 'TI': 0, 'TL': -2, 'TV': 0, \
        'MT': -1, 'IT': 0, 'LT': -2, 'VT': 0, \
        'TF': -3, 'TY': -3, 'TW': -5, \
        'FT': -3, 'YT': -3, 'WT': -5, \
        'PP': 6, 'PA': 1, 'PG': -1, \
        'AP': 1, 'GP': -1, \
        'PN': -1, 'PD': -1, 'PE': -1, 'PQ': 0, \
        'NP': -1, 'DP': -1, 'EP': -1, 'QP': 0, \
        'PH': 0, 'PR': 0, 'PK': -1, \
        'HP': 0, 'RP': 0, 'KP': -1, \
        'PM': -2, 'PI': -2, 'PL': -3, 'PV': -1, \
        'MP': -2, 'IP': -2, 'LP': -3, 'VP': -1, \
        'PF': -5, 'PY': -5, 'PW': -6, \
        'FP': -5, 'YP': -5, 'WP': -6, \
        'AA': 2, 'AG': 1, \
        'GA': 1, \
        'AN': 0, 'AD': 0, 'AE': 0, 'AQ': 0, \
        'NA': 0, 'DA': 0, 'EA': 0, 'QA': 0, \
        'AH': -1, 'AR': -2, 'AK': -1, \
        'HA': -1, 'RA': -2, 'KA': -1, \
        'AM': -1, 'AI': -1, 'AL': -2, 'AV': 0, \
        'MA': -1, 'IA': -1, 'LA': -2, 'VA': 0, \
        'AF': -4, 'AY': -3, 'AW': -6, \
        'FA': -4, 'YA': -3, 'WA': -6, \
        'GG': 5, \
        'GN': 0, 'GD': 1, 'GE': 0, 'GQ': -1, \
        'NG': 0, 'DG': 1, 'EG': 0, 'QG': -1, \
        'GH': -2, 'GR': -3, 'GK': -2, \
        'HG': -2, 'RG': -3, 'KG': -2, \
        'GM': -3, 'GI': -3, 'GL': -4, 'GV': -1, \
        'MG': -3, 'IG': -3, 'LG': -4, 'VG': -1, \
        'GF': -5, 'GY': -7, 'GW': -7, \
        'FG': -5, 'YG': -7, 'WG': -7, \
        'NN': 2, 'ND': 2, 'NE': 1, 'NQ': 1, \
        'DN': 2, 'EN': 1, 'QN': 1, \
        'NH': 2, 'NR': 0, 'NK': 1, \
        'HN': 2, 'RN': 0, 'KN': 1, \
        'NM': -2, 'NI': -2, 'NL': -3, 'NV': -2, \
        'MN': -2, 'IN': -2, 'LN': -3, 'VN': -2, \
        'NF': -2, 'NY': -2, 'NW': -4, \
        'FN': -2, 'YN': -2, 'WN': -4, \
        'DD': 4, 'DE': 3, 'DQ': 2, \
        'ED': 3, 'QD': 2, \
        'DH': 1, 'DR': -1, 'DK': 0, \
        'HD': 1, 'RD': -1, 'KD': 0, \
        'DM': -3, 'DI': -2, 'DL': -4, 'DV': -2, \
        'MD': -3, 'ID': -2, 'LD': -4, 'VD': -2, \
        'DF': -6, 'DY': -4, 'DW': -7, \
        'FD': -6, 'YD': -4, 'WD': -7, \
        'EE': 4, 'EQ': 2, \
        'QE': 2, \
        'EH': 1, 'ER': -1, 'EK': 0, \
        'HE': 1, 'RE': -1, 'KE': 0, \
        'EM': -2, 'EI': -2, 'EL': -3, 'EV': -2, \
        'ME': -2, 'IE': -2, 'LE': -3, 'VE': -2, \
        'EF': -5, 'EY': -4, 'EW': -7, \
        'FE': -5, 'YE': -4, 'WE': -7, \
        'QQ': 4, \
        'QH': 3, 'QR': 1, 'QK': 1, \
        'HQ': 3, 'RQ': 1, 'KQ': 1, \
        'QM': -1, 'QI': -2, 'QL': -2, 'QM': -2, \
        'MQ': -1, 'IQ': -2, 'LQ': -2, 'MQ': -2, \
        'QF': -5, 'QY': -4, 'QW': -5, \
        'FQ': -5, 'YQ': -4, 'WQ': -5, \
        'HH': 6, 'HR': 2, 'HK': 0, \
        'RH': 2, 'KH': 0, \
        'HM': -2, 'HI': -2, 'HL': -2, 'HV': -2, \
        'MH': -2, 'IH': -2, 'LH': -2, 'VH': -2, \
        'HF': -2, 'HY': 0, 'HW': -3, \
        'FH': -2, 'YH': 0, 'WH': -3, \
        'RR': 6, 'RK': 3, \
        'KR': 3, \
        'RM': 0, 'RI': -2, 'RL': -3, 'RV': -2, \
        'MR': 0, 'IR': -2, 'LR': -3, 'VR': -2, \
        'RF': -4, 'RY': -4, 'RW': 2, \
        'FR': -4, 'YR': -4, 'WR': 2, \
        'KK': 5, \
        'KM': 0, 'KI': -2, 'KL': -3, 'KV': -2, \
        'MK': 0, 'IK': -2, 'LK': -3, 'VK': -2, \
        'KF': -5, 'KY': -4, 'KW': -3, \
        'FK': -5, 'YK': -4, 'WK': -3, \
        'MM': 6, 'MI': 2, 'ML': 4, 'MV': 2, \
        'IM': 2, 'LM': 4, 'VM': 2, \
        'MF': 0, 'MY': -2, 'MW': -4, \
        'FM': 0, 'YM': -2, 'WM': -4, \
        'II': 5, 'IL': 2, 'IV': 4, \
        'LI': 2, 'VI': 4, \
        'IF': 1, 'IY': -1, 'IW': -5, \
        'FI': 1, 'YI': -1, 'WI': -5, \
        'LL': 6, 'LV': 2, \
        'VL': 2, \
        'LF': 2, 'LY': -1, 'LW': -2, \
        'FL': 2, 'YL': -1, 'WL': -2, \
        'VV': 4, \
        'VF': -1, 'VY': -2, 'VW': -6, \
        'FV': -1, 'YV': -2, 'WV': -6, \
        'FF': 9, 'FY': 7, 'FW': 0, \
        'YF': 7, 'WF': 0, \
        'YY': 10, 'YW': 0, \
        'WY': 0
    }

    # BLOSUM-62 matrix
    blosum62Matrix: dict[str, int] = {
        'AA': 4, 'AR': -1, 'AN': -2, 'AD': -2, 'AC': 0, \
        'RA': -1, 'NA': -2, 'DA': -2, 'CA': 0, \
        'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -1, \
        'QA': -1, 'EA': -1, 'GA': 0, 'HA': -2, 'IA': -1, \
        'AL': -1, 'AK': -1, 'AM': -1, 'AF': -2, 'AP': -1, \
        'LA': -1, 'KA': -1, 'MA': -1, 'FA': -2, 'PA': -1, \
        'AS': 1, 'AT': 0, 'AW': -3, 'AY': -2, 'AV': 0, \
        'SA': 1, 'TA': 0, 'WA': -3, 'YA': -2, 'VA': 0, \
        'RR': 5, 'RN': 0, 'RD': -2, 'RC': -3, 'RQ': 1, \
        'NR': 0, 'DR': -2, 'CR': -3, 'QR': 1, \
        'RE': 0, 'RG': -2, 'RH': 0, 'RI': -3, 'RL': -2, \
        'ER': 0, 'GR': -2, 'HR': 0, 'IR': -3, 'LR': -2, \
        'RK': 2, 'RM': -1, 'RF': -3, 'RP': -2, 'RS': -1, \
        'KR': 2, 'MR': -1, 'FR': -3, 'PR': -2, 'SR': -1, \
        'RT': -1, 'RW': -3, 'RY': -2, 'RV': -3, 'NN': 6, \
        'TR': -1, 'WR': -3, 'YR': -2, 'VR': -3, \
        'ND': 1, 'NC': -3, 'NQ': 0, 'NE': 0, 'NG': 0, \
        'DN': 1, 'CN': -3, 'QN': 0, 'EN': 0, 'GN': 0, \
        'NH': 1, 'NI': -3, 'NL': -3, 'NK': 0, 'NM': -2, \
        'HN': 1, 'IN': -3, 'LN': -3, 'KN': 0, 'MN': -2, \
        'NF': -3, 'NP': -2, 'NS': 1, 'NT': 0, 'NW': -4, \
        'FN': -3, 'PN': -2, 'SN': 1, 'TN': 0, 'WN': -4, \
        'NY': -2, 'NV': -3, 'DD': 6, 'DC': -3, 'DQ': 0, \
        'YN': -2, 'VN': -3, 'CD': -3, 'QD': 0, \
        'DE': 2, 'DG': -1, 'DH': -1, 'DI': -3, 'DL': -4, \
        'ED': 2, 'GD': -1, 'HD': -1, 'ID': -3, 'LD': -4, \
        'DK': -1, 'DM': -3, 'DF': -3, 'DP': -1, 'DS': 0, \
        'KD': -1, 'MD': -3, 'FD': -3, 'PD': -1, 'SD': 0, \
        'DT': -1, 'DW': -4, 'DY': -3, 'DV': -3, 'CC': 9, \
        'TD': -1, 'WD': -4, 'YD': -3, 'VD': -3, \
        'CQ': -3, 'CE': -4, 'CG': -3, 'CH': -3, 'CI': -1, \
        'QC': -3, 'EC': -4, 'GC': -3, 'HC': -3, 'IC': -1, \
        'CL': -1, 'CK': -3, 'CM': -1, 'CF': -2, 'CP': -3, \
        'LC': -1, 'KC': -3, 'MC': -1, 'FC': -2, 'PC': -3, \
        'CS': -1, 'CT': -1, 'CW': -2, 'CY': -2, 'CV': -1, \
        'SC': -1, 'TC': -1, 'WC': -2, 'YC': -2, 'VC': -1, \
        'QQ': 5, 'QE': 2, 'QG': -2, 'QH': 0, 'QI': -3, \
        'EQ': 2, 'GQ': -2, 'HQ': 0, 'IQ': -3, \
        'QL': -2, 'QK': 1, 'QM': 0, 'QF': -3, 'QP': -1, \
        'LQ': -2, 'KQ': 1, 'MQ': 0, 'FQ': -3, 'PQ': -1, \
        'QS': 0, 'QT': -1, 'QW': -2, 'QY': -1, 'QV': -2, \
        'SQ': 0, 'TQ': -1, 'WQ': -2, 'YQ': -1, 'VQ': -2, \
        'EE': 5, 'EG': -2, 'EH': 0, 'EI': -3, 'EL': -3, \
        'GE': -2, 'HE': 0, 'IE': -3, 'LE': -3, \
        'EK': 1, 'EM': -2, 'EF': -3, 'EP': -1, 'ES': 0, \
        'KE': 1, 'ME': -2, 'FE': -3, 'PE': -1, 'SE': 0, \
        'ET': -1, 'EW': -3, 'EY': -2, 'EV': -2, 'GG': 6, \
        'TE': -1, 'WE': -3, 'YE': -2, 'VE': -2, \
        'GH': -2, 'GI': -4, 'GL': -4, 'GK': -2, 'GM': -3, \
        'HG': -2, 'IG': -4, 'LG': -4, 'KG': -2, 'MG': -3, \
        'GF': -3, 'GP': -2, 'GS': 0, 'GT': -2, 'GW': -2, \
        'FG': -3, 'PG': -2, 'SG': 0, 'TG': -2, 'WG': -2, \
        'GY': -3, 'GV': -3, 'HH': 8, 'HI': -3, 'HL': -3, \
        'YG': -3, 'VG': -3, 'IH': -3, 'LH': -3, \
        'HK': -1, 'HM': -2, 'HF': -1, 'HP': -2, 'HS': -1, \
        'KH': -1, 'MH': -2, 'FH': -1, 'PH': -2, 'SH': -1, \
        'HT': -2, 'HW': -2, 'HY': 2, 'HV': -3, 'II': 4, \
        'TH': -2, 'WH': -2, 'YH': 2, 'VH': -3, \
        'IL': 2, 'IK': -3, 'IM': 1, 'IF': 0, 'IP': -3, \
        'LI': 2, 'KI': -3, 'MI': 1, 'FI': 0, 'PI': -3, \
        'IS': -2, 'IT': -1, 'IW': -3, 'IY': -1, 'IV': 3, \
        'SI': -2, 'TI': -1, 'WI': -3, 'YI': -1, 'VI': 3, \
        'LL': 4, 'LK': -2, 'LM': 2, 'LF': 0, 'LP': -3, \
        'KL': -2, 'ML': 2, 'FL': 0, 'PL': -3, \
        'LS': -2, 'LT': -1, 'LW': -2, 'LY': -1, 'LV': 1, \
        'SL': -2, 'TL': -1, 'WL': -2, 'YL': -1, 'VL': 1, \
        'KK': 5, 'KM': -1, 'KF': -3, 'KP': -1, 'KS': 0, \
        'MK': -1, 'FK': -3, 'PK': -1, 'SK': 0, \
        'KT': -1, 'KW': -3, 'KY': -2, 'KV': -2, 'MM': 5, \
        'TK': -1, 'WK': -3, 'YK': -2, 'VK': -2, \
        'MF': 0, 'MP': -2, 'MS': -1, 'MT': -1, 'MW': -1, \
        'FM': 0, 'PM': -2, 'SM': -1, 'TM': -1, 'WM': -1, \
        'MY': -1, 'MV': 1, 'FF': 6, 'FP': -4, 'FS': -2, \
        'YM': -1, 'VM': 1, 'PF': -4, 'SF': -2, \
        'FT': -2, 'FW': 1, 'FY': 3, 'FV': -1, 'PP': 7, \
        'TF': -2, 'WF': 1, 'YF': 3, 'VF': -1, \
        'PS': -1, 'PT': -1, 'PW': -4, 'PY': -3, 'PV': -2, \
        'SP': -1, 'TP': -1, 'WP': -4, 'YP': -3, 'VP': -2, \
        'SS': 4, 'ST': 1, 'SW': -3, 'SY': -2, 'SV': -2, \
        'TS': 1, 'WS': -3, 'YS': -2, 'VS': -2, \
        'TT': 5, 'TW': -2, 'TY': -2, 'TV': 0, 'WW': 11, \
        'WT': -2, 'YT': -2, 'VT': 0, \
        'WY': 2, 'WV': -3, 'YY': 7, 'YV': -1, 'VV': 4, \
        'YW': 2, 'VW': -3, 'VY': -1
    }

    def __init__(self):
        pass


# Read fasta files
def readFasta(fastaFilePath: str) -> OrderedDict[str, str]:
    """To read FASTA file and get its identifiers and sequences. The output is a ordered dictionary whose keys are identifiers, values are sequences

    Args:
        fastaFilePath (str): Path of FASTA file.

    Returns:
        OrderedDict[str, str]: Keys are identifiers. And values are sequences.

    """

    outputDict: OrderedDict[str, str] = OrderedDict()
    with open(file=fastaFilePath, mode='r') as f:
        line: str = f.readline()
        tax: str = ''
        seq: str = ''
        if (line.startswith('>')):
            tax = line.lstrip('>').rstrip('\n')
        while (line):
            line = f.readline()
            if (line.startswith('>')):
                if (tax != ''):
                    outputDict[tax] = seq
                tax = line.lstrip('>').rstrip('\n')
                seq = ''
            else:
                seq += line.rstrip('\n')
        outputDict[tax] = seq
    return (outputDict)


# transcription
def dna2Rna(dnaSeq: str) -> str:
    """To transcribe DNA into RNA.

    Args:
        dnaSeq (str): DNA sequence. Capitalization does not matter.

    Returns:
        str: RNA sequence. All letters are uppercase.
    """

    convKit: ConvKit = ConvKit()
    output = dnaSeq.upper().translate(str.maketrans(convKit.dna2RnaDict))
    return (output)


# reverse transcription
def rna2Dna(rnaSeq: str) -> str:
    """To reverse transcribe RNA sequence into DNA sequence.

    Args:
        rnaSeq (str): RNA sequence. Capitalization does not matter.

    Returns:
        str: DNA sequence. All letters are uppercase.
    """

    convKit: ConvKit = ConvKit()
    output = rnaSeq.upper().translate(str.maketrans(convKit.rna2DnaDict))
    return (output)


# DNA translation
def dna2Pro(dnaSeq: str, start: int = 0, end: int = -1) -> str:
    """To translate DNA sequence into protein sequence directly.

    Args:
        dnaSeq (str): DNA sequence. Capitalization does not matter.
        start (int): Beginning site of the translation.
        end (int): Ending site of the translation.

    Returns:
        str: Protein sequence. All letters are uppercase.
    """

    if (end == -1):
        end = len(dnaSeq) - 1
    assert (end > start), 'Invalid arguments!'
    convKit: ConvKit = ConvKit()
    i: int = int(start)
    outputList: list[str] = []
    while (i <= end - 2):
        outputList.append(convKit.dna2ProDict[dnaSeq[i: i+3].upper()])
        i += 3
    output: str = ''.join(outputList)
    return (output)


# RNA translation
def rna2Pro(rnaSeq: str, start: int = 0, end: int = -1) -> str:
    """To translate RNA sequence into protein sequence.

    Args:
        rnaSeq (str): RNA sequence. Capitalization does not matter.
        start (int): Beginning site of the translation.
        end (int): Ending site of the translation.

    Returns:
        str: Protein sequence. All letters are uppercase.
    """

    if (end == -1):
        end = len(rnaSeq) - 1
    assert (end > start), 'Invalid arguments!'
    convKit: ConvKit = ConvKit()
    i: int = int(start)
    outputList: list[str] = []
    while (i <= end - 2):
        outputList.append(convKit.rna2ProDict[rnaSeq[i: i+3].upper()])
        i += 3
    output: str = ''.join(outputList)
    return (output)


# Alignment
def pairwiseDnaAlign(fasta: str = '', seqA: str = 'ACGT', seqB: str = 'ACGTA', matrix: str = 'unitary', gapOpen: float = -10, gapExtend: float = -0.5) -> None:
    """Pairwise DNA Alignment. Results will be output in console.

    Args:
        fasta (str): FASTA file with two sequences.
        seqA (str): If no FASTA file is given. Sequence can be put into argument seqA.
        seqB (str): If no FASTA file is given. Sequence can be put into argument seqB.
        matrix (str): score matrix. Can be 'unitary', 'blast' or 'tt'.
        gapOpen (float): Punishment for opening a gap.
        gapExtend (float): Punishment for extending a gap.
    """

    if (fasta != ''):
        sequences: OrderedDict[str, str] = readFasta(fastaFilePath=fasta)
        seqA: str = list(sequences.values())[0]
        seqB: str = list(sequences.values())[1]
    seqA = seqA.upper()
    seqB = seqB.upper()
    # sequence cleaning (remove elements excluding [A, C, G, T])
    def dnaSeqClean(seq: str) -> str:
        seqOut: str = ''.join(list(map(lambda x:x if (x in ['A', 'C', 'G', 'T']) else '', seq)))
        return (seqOut)

    seqA = dnaSeqClean(seq=seqA)
    seqB = dnaSeqClean(seq=seqB)
    assert (len(seqA) * len(seqB) != 0), 'empty string or not DNA sequence'

    convKit: ConvKit = ConvKit()
    scoreMatrix: list[list[float]] = []
    for i in range(len(seqB)+1):
        scoreMatrix.append([])
        for j in range(len(seqA)+1):
            scoreMatrix[i].append(0)
    # Needleman-Wunsch (NW)

    # Matrix initiation
    for rowIdx in range(1, len(seqB)+1):
        scoreMatrix[rowIdx][0] = gapOpen + gapExtend * (rowIdx - 1)
    for columnIdx in range(len(seqA)+1):
        scoreMatrix[0][columnIdx] = gapOpen + gapExtend * (columnIdx - 1)

    actionArray: list[list[bool]] = [[False, True] for columnIdx in range(len(seqA)+2)]
    actionArray[0] = [True, False]

    # Matrix scanning
    match (matrix):
        case ('unitary'):
            for rowIdx in range(1, len(seqB)+1):
                buffer: list[list[bool]] = [[True, False]]
                for columnIdx in range(1, len(seqA)+1):
                    # thisR
                    if (actionArray[0][1]):
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapExtend
                    else:
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapOpen

                    # thisD
                    if (actionArray[columnIdx+1][0]):
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapExtend
                    else:
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapOpen

                    # thisC
                    thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.dnaUnitaryMatrix[(seqA[columnIdx-1]==seqB[rowIdx-1])]

                    # Action decision
                    scoreMatrix[rowIdx][columnIdx] = max(thisR, thisD, thisC)
                    if (thisC == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, False]
                    elif (thisR == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, True]
                    else:
                        actionArray[0] = [True, False]
                    buffer.append(actionArray[0])

                # buffer added into actionArray
                actionArray = [[True, False], [True, False]]
                actionArray.extend(buffer)
                del(buffer)

        case ('blast'):
            for rowIdx in range(1, len(seqB)+1):
                buffer: list[list[bool]] = [[True, False]]
                for columnIdx in range(1, len(seqA)+1):
                    # thisR
                    if (actionArray[0][1]):
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapExtend
                    else:
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapOpen

                    # thisD
                    if (actionArray[columnIdx+1][0]):
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapExtend
                    else:
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapOpen

                    # thisC
                    thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.dnaBlastMatrix[(seqA[columnIdx-1]==seqB[rowIdx-1])]

                    # Action decision
                    scoreMatrix[rowIdx][columnIdx] = max(thisR, thisD, thisC)
                    if (thisC == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, False]
                    elif (thisR == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, True]
                    else:
                        actionArray[0] = [True, False]
                    buffer.append(actionArray[0])

                # buffer added into actionArray
                actionArray = [[True, False], [True, False]]
                actionArray.extend(buffer)
                del(buffer)

        case ('tt'):
            purine:list[str] = ['A', 'G']
            pyrimidine:list[str] = ['C', 'T']
            for rowIdx in range(1, len(seqB)+1):
                buffer: list[list[bool]] = [[True, False]]
                for columnIdx in range(1, len(seqA)+1):
                    # thisR
                    if (actionArray[0][1]):
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapExtend
                    else:
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapOpen

                    # thisD
                    if (actionArray[columnIdx+1][0]):
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapExtend
                    else:
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapOpen

                    # thisC
                    if (seqB[rowIdx-1] == seqA[columnIdx-1]):
                        thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.dnaTTMatrix[True]
                    elif (((seqB[rowIdx-1] in purine) and (seqA[columnIdx-1] in purine)) or ((seqB[rowIdx-1] in pyrimidine) and seqA[columnIdx-1] in pyrimidine)):
                        thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.dnaTTMatrix['transition']
                    else:
                        thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.dnaTTMatrix['transversion']

                    # Action decision
                    scoreMatrix[rowIdx][columnIdx] = max(thisR, thisD, thisC)
                    if (thisC == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, False]
                    elif (thisR == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, True]
                    else:
                        actionArray[0] = [True, False]
                    buffer.append(actionArray[0])

                # buffer added into actionArray
                actionArray = [[True, False], [True, False]]
                actionArray.extend(buffer)
                del(buffer)

        case (_):
            print('Unknown matrix. \'unitary\', \'blast\' and \'tt\' are supported')
            return (1)

    finalScore: float = scoreMatrix[-1][-1]
    # backforward
    steps: list[str] = []
    rowIdx = len(seqB)
    columnIdx = len(seqA)
    while ((rowIdx != 0) or (columnIdx != 0)):

        if ((rowIdx != 0) and (columnIdx != 0)):
            bestPathScore = max(scoreMatrix[rowIdx-1][columnIdx], \
                                scoreMatrix[rowIdx][columnIdx-1], \
                                    scoreMatrix[rowIdx-1][columnIdx-1])
            # up-leftward
            if (bestPathScore == scoreMatrix[rowIdx-1][columnIdx-1]):
                steps.insert(0, 'c')
                rowIdx -= 1
                columnIdx -= 1
            # upward
            elif (bestPathScore == scoreMatrix[rowIdx-1][columnIdx]):
                steps.insert(0, 'd')
                rowIdx -= 1
            # leftward
            else:
                steps.insert(0, 'r')
                columnIdx -= 1

        elif ((rowIdx == 0) and (columnIdx != 0)):
            steps.insert(0, 'r')
            columnIdx -= 1

        else:
            steps.insert(0, 'd')
            rowIdx -= 1

    seqAIdx: int = 0
    seqBIdx: int = 0
    seqAOut: str = ''
    seqBOut: str = ''
    for step in steps:
        if (step == 'r'):
            seqAOut += seqA[seqAIdx]
            seqBOut += '-'
            seqAIdx += 1
        elif (step == 'd'):
            seqAOut += '-'
            seqBOut += seqB[seqBIdx]
            seqBIdx += 1
        else:
            seqAOut += seqA[seqAIdx]
            seqBOut += seqB[seqBIdx]
            seqBIdx += 1
            seqAIdx += 1


    seqAlign: str = ''
    for idx in range(len(seqAOut)):
        if (seqAOut[idx] == seqBOut[idx]):
            seqAlign += '|'
        else:
            seqAlign += ' '

    consoleWidth: int = os.get_terminal_size().columns
    seqAOut = [seqAOut[i:i+consoleWidth] for i in range(0, len(seqAOut), consoleWidth)]
    seqAlign = [seqAlign[i:i+consoleWidth] for i in range(0, len(seqAlign), consoleWidth)]
    seqBOut = [seqBOut[i:i+consoleWidth] for i in range(0, len(seqBOut), consoleWidth)]
    for idx in range(len(seqAOut)):
        print(seqAOut[idx])
        print(seqAlign[idx])
        print(seqBOut[idx])
    print(finalScore)


def pairwiseProtAlign(fasta: str = '', seqA: str = 'ACGT', seqB: str = 'ACGTT', matrix: str = 'unitary', gapOpen: float = -10, gapExtend: float = -0.5, consoleWidth: int = 50) -> None:
    """Pairwise Protein Alignment. Results will be output in console.

    Args:
        fasta (str): FASTA file with two sequences.
        seqA (str): If no FASTA file is given. Sequence can be put into argument seqA.
        seqB (str): If no FASTA file is given. Sequence can be put into argument seqB.
        matrix (str): score matrix. Can be 'unitary', 'blosum62' or 'pam250'.
        gapOpen (float): Punishment for opening a gap.
        gapExtend (float): Punishment for extending a gap.
        consoleWidth (int): This determine how many bases shown in a single line.
    """
    if (fasta != ''):
        sequences: dict[str, str] = readFasta(fastaFilePath=fasta)
    seqA: str = list(sequences.values())[0]
    seqB: str = list(sequences.values())[1]

    seqA = seqA.upper()
    seqB = seqB.upper()

    convKit: ConvKit = ConvKit()
    scoreMatrix: list[list[float]] = []
    for i in range(len(seqB)+1):
        scoreMatrix.append([])
        for j in range(len(seqA)+1):
            scoreMatrix[i].append(0)
    # Needleman-Wunsch (NW)

    # Matrix scanning
    match (matrix):
        case ('unitary'):

            # Matrix initiation
            for rowIdx in range(1, len(seqB)+1):
                scoreMatrix[rowIdx][0] = gapOpen + gapExtend * (rowIdx - 1)
            for columnIdx in range(len(seqA)+1):
                scoreMatrix[0][columnIdx] = gapOpen + gapExtend * (columnIdx - 1)

            actionArray: list[list[bool]] = [[False, True] for columnIdx in range(len(seqA)+2)]
            actionArray[0] = [True, False]

            for rowIdx in range(1, len(seqB)+1):
                buffer: list[list[bool]] = [[True, False]]
                for columnIdx in range(1, len(seqA)+1):
                    # thisR
                    if (actionArray[0][1]):
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapExtend
                    else:
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapOpen

                    # thisD
                    if (actionArray[columnIdx+1][0]):
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapExtend
                    else:
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapOpen

                    # thisC
                    thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.protUnitaryMatrix[(seqA[columnIdx-1]==seqB[rowIdx-1])]

                    # Action decision
                    scoreMatrix[rowIdx][columnIdx] = max(thisR, thisD, thisC)
                    if (thisC == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, False]
                    elif (thisR == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, True]
                    else:
                        actionArray[0] = [True, False]
                    buffer.append(actionArray[0])

                # buffer added into actionArray
                actionArray = [[True, False], [True, False]]
                actionArray.extend(buffer)
                del(buffer)

        case ('pam250'):

            # sequence cleaning (remove elements excluding [A, C, G, T])
            def dnaSeqClean(seq: str) -> str:
                seqOut: str = ''.join(list(map(lambda x:x if (x in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']) else '', seq)))
                return (seqOut)

            seqA = dnaSeqClean(seq=seqA)
            seqB = dnaSeqClean(seq=seqB)
            assert (len(seqA) * len(seqB) != 0), 'pam250 can only be used for 20 common residues'

            # Matrix initiation
            for rowIdx in range(1, len(seqB)+1):
                scoreMatrix[rowIdx][0] = gapOpen + gapExtend * (rowIdx - 1)
            for columnIdx in range(len(seqA)+1):
                scoreMatrix[0][columnIdx] = gapOpen + gapExtend * (columnIdx - 1)

            actionArray: list[list[bool]] = [[False, True] for columnIdx in range(len(seqA)+2)]
            actionArray[0] = [True, False]

            for rowIdx in range(1, len(seqB)+1):
                buffer: list[list[bool]] = [[True, False]]
                for columnIdx in range(1, len(seqA)+1):
                    # thisR
                    if (actionArray[0][1]):
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapExtend
                    else:
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapOpen

                    # thisD
                    if (actionArray[columnIdx+1][0]):
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapExtend
                    else:
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapOpen

                    # thisC
                    thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.pam250Matrix[seqA[columnIdx-1]+seqB[rowIdx-1]]

                    # Action decision
                    scoreMatrix[rowIdx][columnIdx] = max(thisR, thisD, thisC)
                    if (thisC == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, False]
                    elif (thisR == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, True]
                    else:
                        actionArray[0] = [True, False]
                    buffer.append(actionArray[0])

                # buffer added into actionArray
                actionArray = [[True, False], [True, False]]
                actionArray.extend(buffer)
                del(buffer)

        case ('blosum62'):

            # sequence cleaning (remove elements excluding [A, C, G, T])
            def dnaSeqClean(seq: str) -> str:
                seqOut: str = ''.join(list(map(lambda x:x if (x in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']) else '', seq)))
                return (seqOut)

            seqA = dnaSeqClean(seq=seqA)
            seqB = dnaSeqClean(seq=seqB)
            assert (len(seqA) * len(seqB) != 0), 'blosum62 can only be used for 20 common residues'

            # Matrix initiation
            for rowIdx in range(1, len(seqB)+1):
                scoreMatrix[rowIdx][0] = gapOpen + gapExtend * (rowIdx - 1)
            for columnIdx in range(len(seqA)+1):
                scoreMatrix[0][columnIdx] = gapOpen + gapExtend * (columnIdx - 1)

            actionArray: list[list[bool]] = [[False, True] for columnIdx in range(len(seqA)+2)]
            actionArray[0] = [True, False]

            for rowIdx in range(1, len(seqB)+1):
                buffer: list[list[bool]] = [[True, False]]
                for columnIdx in range(1, len(seqA)+1):
                    # thisR
                    if (actionArray[0][1]):
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapExtend
                    else:
                        thisR: float = scoreMatrix[rowIdx][columnIdx-1] + gapOpen

                    # thisD
                    if (actionArray[columnIdx+1][0]):
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapExtend
                    else:
                        thisD: float = scoreMatrix[rowIdx-1][columnIdx] + gapOpen

                    # thisC
                    thisC: float = scoreMatrix[rowIdx-1][columnIdx-1] + convKit.blosum62Matrix[seqA[columnIdx-1]+seqB[rowIdx-1]]

                    # Action decision
                    scoreMatrix[rowIdx][columnIdx] = max(thisR, thisD, thisC)
                    if (thisC == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, False]
                    elif (thisR == max(thisR, thisD, thisC)):
                        actionArray[0] = [False, True]
                    else:
                        actionArray[0] = [True, False]
                    buffer.append(actionArray[0])

                # buffer added into actionArray
                actionArray = [[True, False], [True, False]]
                actionArray.extend(buffer)
                del(buffer)

        case (_):
            print('Unknown matrix. \'unitary\', \'pam250\' and \'blosum62\' are supported')
            return (1)

    finalScore: float = scoreMatrix[-1][-1]
    # backforward
    steps: list[str] = []
    rowIdx = len(seqB)
    columnIdx = len(seqA)
    while ((rowIdx != 0) or (columnIdx != 0)):

        if ((rowIdx != 0) and (columnIdx != 0)):
            bestPathScore = max(scoreMatrix[rowIdx-1][columnIdx], \
                                scoreMatrix[rowIdx][columnIdx-1], \
                                    scoreMatrix[rowIdx-1][columnIdx-1])
            # up-leftward
            if (bestPathScore == scoreMatrix[rowIdx-1][columnIdx-1]):
                steps.insert(0, 'c')
                rowIdx -= 1
                columnIdx -= 1
            # upward
            elif (bestPathScore == scoreMatrix[rowIdx-1][columnIdx]):
                steps.insert(0, 'd')
                rowIdx -= 1
            # leftward
            else:
                steps.insert(0, 'r')
                columnIdx -= 1

        elif ((rowIdx == 0) and (columnIdx != 0)):
            steps.insert(0, 'r')
            columnIdx -= 1

        else:
            steps.insert(0, 'd')
            rowIdx -= 1

    seqAIdx: int = 0
    seqBIdx: int = 0
    seqAOut: str = ''
    seqBOut: str = ''
    for step in steps:
        if (step == 'r'):
            seqAOut += seqA[seqAIdx]
            seqBOut += '-'
            seqAIdx += 1
        elif (step == 'd'):
            seqAOut += '-'
            seqBOut += seqB[seqBIdx]
            seqBIdx += 1
        else:
            seqAOut += seqA[seqAIdx]
            seqBOut += seqB[seqBIdx]
            seqBIdx += 1
            seqAIdx += 1


    seqAlign: str = ''
    for idx in range(len(seqAOut)):
        if (seqAOut[idx] == seqBOut[idx]):
            seqAlign += '|'
        else:
            seqAlign += ' '


    seqAOut = [seqAOut[i:i+consoleWidth] for i in range(0, len(seqAOut), consoleWidth)]
    seqAlign = [seqAlign[i:i+consoleWidth] for i in range(0, len(seqAlign), consoleWidth)]
    seqBOut = [seqBOut[i:i+consoleWidth] for i in range(0, len(seqBOut), consoleWidth)]
    for idx in range(len(seqAOut)):
        print(seqAOut[idx])
        print(seqAlign[idx])
        print(seqBOut[idx])
    print(finalScore)






