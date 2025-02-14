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
    Estimates the boiling point using the Joback method.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define Joback group contributions for boiling point estimation
    joback_groups = {
        "C_C": {"smarts": "[CX4H3]", "delta_tb": -23.58},
        "C_CH": {"smarts": "[CX4H2][CX4]", "delta_tb": -7.24},
        "C_C=C": {"smarts": "[CX3]=[CX3]", "delta_tb": 14.73},
        "C_C#C": {"smarts": "[CX2]#[CX2]", "delta_tb": 5.15},
        "C_Cl": {"smarts": "[ClX1][CX4]", "delta_tb": 26.86},
        "C_Br": {"smarts": "[BrX1][CX4]", "delta_tb": 68.10},
        "C_I": {"smarts": "[IX1][CX4]", "delta_tb": 100.00},
        "O_H": {"smarts": "[OX2H]", "delta_tb": 33.60},
        "O_C": {"smarts": "[OX2H][CX4]", "delta_tb": 21.78},
        "N_H2": {"smarts": "[NX3H2]", "delta_tb": 35.00},
        # Add more functional groups as needed
    }
    
    # Initialize delta_tb sum
    delta_tb_sum = 0.0
    groups_found = []
    
    # Loop over Joback groups to find matches and sum delta_tb
    for group_name, group_info in joback_groups.items():
        smarts = group_info["smarts"]
        delta_tb = group_info["delta_tb"]
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        num_matches = len(matches)
        if num_matches > 0:
            delta_tb_sum += delta_tb * num_matches
            groups_found.append((group_name, num_matches, delta_tb))
    
    # Base boiling point
    Tb = 198.0 + delta_tb_sum  # in Kelvin
    
    # Convert to Celsius
    Tb_C = Tb - 273.15
    
    # Classify based on boiling point
    if Tb_C <= 250:
        result = True
        reason = f"Estimated boiling point is {Tb_C:.2f}°C, which is ≤ 250°C"
    else:
        result = False
        reason = f"Estimated boiling point is {Tb_C:.2f}°C, which is > 250°C"
    
    return result, reason