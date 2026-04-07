"""Amino acid data, Grantham matrix, and classification helpers."""

import pandas as pd

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

AA_PROPERTIES = {
    'A': {'size': 'small',  'charge': 'neutral',  'hydrophobic': True},
    'R': {'size': 'large',  'charge': 'positive', 'hydrophobic': False},
    'N': {'size': 'medium', 'charge': 'neutral',  'hydrophobic': False},
    'D': {'size': 'medium', 'charge': 'negative', 'hydrophobic': False},
    'C': {'size': 'small',  'charge': 'neutral',  'hydrophobic': True},
    'E': {'size': 'medium', 'charge': 'negative', 'hydrophobic': False},
    'Q': {'size': 'medium', 'charge': 'neutral',  'hydrophobic': False},
    'G': {'size': 'small',  'charge': 'neutral',  'hydrophobic': False},
    'H': {'size': 'medium', 'charge': 'positive', 'hydrophobic': False},
    'I': {'size': 'medium', 'charge': 'neutral',  'hydrophobic': True},
    'L': {'size': 'medium', 'charge': 'neutral',  'hydrophobic': True},
    'K': {'size': 'large',  'charge': 'positive', 'hydrophobic': False},
    'M': {'size': 'medium', 'charge': 'neutral',  'hydrophobic': True},
    'F': {'size': 'large',  'charge': 'neutral',  'hydrophobic': True},
    'P': {'size': 'small',  'charge': 'neutral',  'hydrophobic': False},
    'S': {'size': 'small',  'charge': 'neutral',  'hydrophobic': False},
    'T': {'size': 'small',  'charge': 'neutral',  'hydrophobic': False},
    'W': {'size': 'large',  'charge': 'neutral',  'hydrophobic': True},
    'Y': {'size': 'large',  'charge': 'neutral',  'hydrophobic': False},
    'V': {'size': 'small',  'charge': 'neutral',  'hydrophobic': True},
}

# Max SASA (Tien et al 2013, theoretical Gly-X-Gly)
MAX_SASA = {
    'A': 129, 'R': 274, 'N': 195, 'D': 193, 'C': 167, 'E': 223,
    'Q': 225, 'G': 104, 'H': 224, 'I': 197, 'L': 201, 'K': 236,
    'M': 224, 'F': 240, 'P': 159, 'S': 155, 'T': 172, 'V': 174,
    'W': 285, 'Y': 263,
}

# Full Grantham matrix (upper triangle)
_GRANTHAM_RAW = {
    ('A','R'):112,('A','N'):111,('A','D'):126,('A','C'):195,('A','Q'):91,
    ('A','E'):107,('A','G'):60,('A','H'):86,('A','I'):94,('A','L'):96,
    ('A','K'):106,('A','M'):84,('A','F'):113,('A','P'):27,('A','S'):99,
    ('A','T'):58,('A','W'):148,('A','Y'):112,('A','V'):64,
    ('R','N'):86,('R','D'):96,('R','C'):180,('R','Q'):43,('R','E'):54,
    ('R','G'):125,('R','H'):29,('R','I'):97,('R','L'):102,('R','K'):26,
    ('R','M'):91,('R','F'):97,('R','P'):103,('R','S'):110,('R','T'):71,
    ('R','W'):101,('R','Y'):77,('R','V'):96,
    ('N','D'):23,('N','C'):139,('N','Q'):46,('N','E'):42,('N','G'):80,
    ('N','H'):68,('N','I'):149,('N','L'):153,('N','K'):94,('N','M'):142,
    ('N','F'):158,('N','P'):91,('N','S'):46,('N','T'):65,('N','W'):174,
    ('N','Y'):143,('N','V'):133,
    ('D','C'):154,('D','Q'):61,('D','E'):45,('D','G'):94,('D','H'):81,
    ('D','I'):168,('D','L'):172,('D','K'):101,('D','M'):160,('D','F'):177,
    ('D','P'):108,('D','S'):65,('D','T'):85,('D','W'):181,('D','Y'):160,
    ('D','V'):152,
    ('C','Q'):154,('C','E'):170,('C','G'):159,('C','H'):174,('C','I'):198,
    ('C','L'):198,('C','K'):202,('C','M'):196,('C','F'):205,('C','P'):169,
    ('C','S'):112,('C','T'):149,('C','W'):215,('C','Y'):194,('C','V'):192,
    ('Q','E'):29,('Q','G'):87,('Q','H'):24,('Q','I'):109,('Q','L'):113,
    ('Q','K'):53,('Q','M'):101,('Q','F'):116,('Q','P'):76,('Q','S'):68,
    ('Q','T'):42,('Q','W'):130,('Q','Y'):99,('Q','V'):96,
    ('E','G'):98,('E','H'):40,('E','I'):134,('E','L'):138,('E','K'):56,
    ('E','M'):126,('E','F'):140,('E','P'):93,('E','S'):80,('E','T'):65,
    ('E','W'):152,('E','Y'):122,('E','V'):121,
    ('G','H'):98,('G','I'):135,('G','L'):138,('G','K'):127,('G','M'):127,
    ('G','F'):153,('G','P'):42,('G','S'):56,('G','T'):59,('G','W'):184,
    ('G','Y'):147,('G','V'):109,
    ('H','I'):94,('H','L'):99,('H','K'):32,('H','M'):87,('H','F'):100,
    ('H','P'):77,('H','S'):89,('H','T'):47,('H','W'):115,('H','Y'):83,
    ('H','V'):84,
    ('I','L'):5,('I','K'):102,('I','M'):10,('I','F'):21,('I','P'):95,
    ('I','S'):142,('I','T'):89,('I','W'):61,('I','Y'):33,('I','V'):29,
    ('L','K'):107,('L','M'):15,('L','F'):22,('L','P'):98,('L','S'):145,
    ('L','T'):92,('L','W'):61,('L','Y'):36,('L','V'):32,
    ('K','M'):95,('K','F'):102,('K','P'):103,('K','S'):121,('K','T'):78,
    ('K','W'):110,('K','Y'):85,('K','V'):97,
    ('M','F'):28,('M','P'):87,('M','S'):135,('M','T'):81,('M','W'):67,
    ('M','Y'):36,('M','V'):21,
    ('F','P'):114,('F','S'):155,('F','T'):103,('F','W'):40,('F','Y'):22,
    ('F','V'):50,
    ('P','S'):74,('P','T'):38,('P','W'):147,('P','Y'):110,('P','V'):68,
    ('S','T'):58,('S','W'):177,('S','Y'):144,('S','V'):124,
    ('T','W'):128,('T','Y'):92,('T','V'):69,
    ('W','Y'):37,('W','V'):88,
    ('Y','V'):55,
}

# Build symmetric lookup
GRANTHAM = {}
for (a, b), v in _GRANTHAM_RAW.items():
    GRANTHAM[(a, b)] = v
    GRANTHAM[(b, a)] = v

MAX_GRANTHAM = 215  # C-W


def get_grantham(ref_aa, alt_aa):
    """Get Grantham distance. Returns 0 for identical substitutions."""
    if ref_aa == alt_aa:
        return 0
    return GRANTHAM.get((ref_aa, alt_aa), None)


def classify_grantham(distance):
    if distance is None:
        return None
    if distance <= 60:
        return 'conservative'
    elif distance <= 100:
        return 'moderately_conservative'
    elif distance <= 150:
        return 'moderately_radical'
    return 'radical'


def grantham_severity(distance):
    if distance is None:
        return 0.0
    return distance / MAX_GRANTHAM


def get_property_changes(ref_aa, alt_aa):
    """Describe directional property changes."""
    if ref_aa not in AA_PROPERTIES or alt_aa not in AA_PROPERTIES:
        return ''
    r, a = AA_PROPERTIES[ref_aa], AA_PROPERTIES[alt_aa]
    changes = []
    if r['size'] != a['size']:
        changes.append(f"size:{r['size']}→{a['size']}")
    if r['charge'] != a['charge']:
        changes.append(f"charge:{r['charge']}→{a['charge']}")
    if r['hydrophobic'] != a['hydrophobic']:
        changes.append(f"hydro:{'yes' if r['hydrophobic'] else 'no'}→{'yes' if a['hydrophobic'] else 'no'}")
    return ';'.join(changes)


# ─── Classification helpers ───────────────────────────────────────────────────

def classify_plddt(p, thresholds=None):
    if p is None or (isinstance(p, float) and pd.isna(p)):
        return None
    if thresholds is None:
        thresholds = {'very_high': 90, 'high': 70, 'medium': 50}
    p = float(p)
    if p >= thresholds['very_high']:
        return 'very_high'
    elif p >= thresholds['high']:
        return 'high'
    elif p >= thresholds['medium']:
        return 'medium'
    return 'low'


def classify_contacts(c):
    if c is None:
        return None
    c = int(c)
    if c >= 10:
        return 'high_contact'
    elif c >= 4:
        return 'medium_contact'
    elif c >= 1:
        return 'low_contact'
    return 'no_contacts'


def classify_burial(a):
    if a is None or (isinstance(a, float) and pd.isna(a)):
        return 'unknown'
    a = float(a)
    if a < 0.05:
        return 'buried_core'
    elif a < 0.25:
        return 'partially_buried'
    return 'surface_exposed'


def classify_ddg(v):
    """Classify a DDG value into a semantic category."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    v = float(v)
    if v > 2.0:
        return 'highly_destabilizing'
    elif v > 1.0:
        return 'destabilizing'
    elif v > 0.5:
        return 'mildly_destabilizing'
    elif v > -0.5:
        return 'neutral'
    elif v > -1.0:
        return 'mildly_stabilizing'
    elif v > -2.0:
        return 'stabilizing'
    return 'highly_stabilizing'


# Safe accessors
def sf(v, default=0.0):
    """Safe float."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return default
    try:
        return float(v)
    except (ValueError, TypeError):
        return default


def ss(v):
    """Safe string."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ''
    return str(v).strip()


def sb(v):
    """Safe bool."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return False
    if isinstance(v, bool):
        return v
    return str(v).lower() in ('true', '1', 'yes')


BURIAL_RANK = {'unknown': 0, 'surface_exposed': 1, 'partially_buried': 2, 'buried_core': 3}
RANK_TO_BURIAL = {v: k for k, v in BURIAL_RANK.items()}
