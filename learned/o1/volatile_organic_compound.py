"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound (VOC)
Definition: Any organic compound having an initial boiling point less than or equal to 250°C measured at a standard atmospheric pressure of 101.3 kPa.
"""
from rdkit import Chem

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    Estimates the boiling point using an improved Joback method.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        return False, "Invalid SMILES string"

    # Define Joback group contributions for boiling point estimation
    joback_groups = {
        # Alkanes
        "CH3": {"smarts": "[C;H3]", "delta_tb": -23.58},
        "CH2": {"smarts": "[C;H2]", "delta_tb": -7.24},
        "CH": {"smarts": "[C;H1]", "delta_tb": 9.17},
        "C_alkane": {"smarts": "[C;H0]", "delta_tb": 21.87},
        
        # Alkenes
        "C=C": {"smarts": "[C]=[C]", "delta_tb": 14.70},
        "CH2=CH2": {"smarts": "[C;H2]=[C;H2]", "delta_tb": -10.35},
        "CH=CH": {"smarts": "[C;H1]=[C;H1]", "delta_tb": 5.40},
        "C=CH2": {"smarts": "[C;H1]=[C;H2]", "delta_tb": -10.35},

        # Alkynes
        "C#C": {"smarts": "[C]#[C]", "delta_tb": 5.15},

        # Aromatics
        "Carom_H": {"smarts": "[cH]", "delta_tb": -12.51},
        "Carom": {"smarts": "[c;H0]", "delta_tb": 4.01},

        # Alcohols
        "OH(alcohol)": {"smarts": "[OX2H][C]", "delta_tb": 26.86},
        "OH(phenol)": {"smarts": "[OX2H][c]", "delta_tb": 9.27},

        # Ethers
        "Ether(aliphatic)": {"smarts": "[OD2]([C])[C]", "delta_tb": 13.60},
        "Ether(aromatic)": {"smarts": "[OD2R][c]", "delta_tb": 13.60},

        # Aldehydes and Ketones
        "Carbonyl(aldehyde)": {"smarts": "[CX3H1](=O)[#6]", "delta_tb": 11.70},
        "Carbonyl(ketone)": {"smarts": "[#6][CX3](=O)[#6]", "delta_tb": 7.78},

        # Carboxylic Acids and Esters
        "COOH": {"smarts": "[CX3](=O)[OX2H1]", "delta_tb": 20.00},
        "Ester": {"smarts": "[CX3](=O)[OX2][#6]", "delta_tb": 7.43},

        # Amines
        "PrimaryAmine": {"smarts": "[NX3;H2][#6]", "delta_tb": 48.50},
        "SecondaryAmine": {"smarts": "[NX3;H1]([#6])[#6]", "delta_tb": 9.92},
        "TertiaryAmine": {"smarts": "[NX3;H0]([#6])([#6])[#6]", "delta_tb": 21.60},

        # Halides
        "Fluorine": {"smarts": "[F]", "delta_tb": -18.95},
        "Chlorine": {"smarts": "[Cl]", "delta_tb": 26.86},
        "Bromine": {"smarts": "[Br]", "delta_tb": 68.10},
        "Iodine": {"smarts": "[I]", "delta_tb": 98.72},

        # Nitro
        "Nitro": {"smarts": "[NX3](=O)=O", "delta_tb": 33.60},

        # Thiols
        "Thiol": {"smarts": "[SX2H]", "delta_tb": 21.00},

        # Sulfides
        "Sulfide": {"smarts": "[#16X2]([#6])[#6]", "delta_tb": 21.00},

        # Rings (approximate contribution)
        "Ring": {"smarts": "[R]", "delta_tb": 19.81},
        
        # Others can be added as needed
    }

    # Initialize delta_tb sum
    delta_tb_sum = 0.0
    groups_found = []

    # Keep track of matched atoms
    matched_atoms = set()

    # Loop over Joback groups to find matches and sum delta_tb
    for group_name, group_info in joback_groups.items():
        pattern = Chem.MolFromSmarts(group_info["smarts"])
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern, uniquify=True)
        num_matches = 0
        for match in matches:
            # Check if atoms have been matched already
            if not any(atom_idx in matched_atoms for atom_idx in match):
                num_matches += 1
                matched_atoms.update(match)
        if num_matches > 0:
            delta_tb = group_info["delta_tb"]
            delta_tb_sum += delta_tb * num_matches
            groups_found.append((group_name, num_matches, delta_tb))

    # Account for any unmatched carbon atoms (as CH, CH2, or CH3)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in matched_atoms:
            hydrogen_count = atom.GetTotalNumHs()
            if hydrogen_count == 3:
                group_name = "CH3"
                delta_tb = -23.58
            elif hydrogen_count == 2:
                group_name = "CH2"
                delta_tb = -7.24
            elif hydrogen_count == 1:
                group_name = "CH"
                delta_tb = 9.17
            else:
                group_name = "C_alkane"
                delta_tb = 21.87
            delta_tb_sum += delta_tb
            groups_found.append((group_name, 1, delta_tb))
            matched_atoms.add(atom.GetIdx())

    # Base boiling point
    Tb_K = 198.0 + delta_tb_sum  # in Kelvin

    # Convert to Celsius
    Tb_C = Tb_K - 273.15

    # Classify based on boiling point
    if Tb_C <= 250:
        result = True
        reason = f"Estimated boiling point is {Tb_C:.2f}°C, which is ≤ 250°C"
    else:
        result = False
        reason = f"Estimated boiling point is {Tb_C:.2f}°C, which is > 250°C"

    return result, reason