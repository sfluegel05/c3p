"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: volatile organic compound
"""
from rdkit import Chem

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is any organic compound having an initial boiling point less than or equal to 250°C
    measured at a standard atmospheric pressure of 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import rdMolDescriptors

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains carbon atoms (organic compound)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not an organic compound (contains no carbon atoms)"

    # Define Joback group contributions for boiling point estimation
    joback_groups = [
        {'name': 'CH3', 'smarts': '[CX4H3]', 'deltaTb': 23.58},
        {'name': 'CH2', 'smarts': '[CX4H2]', 'deltaTb': 22.88},
        {'name': 'CH', 'smarts': '[CX4H1]', 'deltaTb': 21.74},
        {'name': 'C', 'smarts': '[CX4H0]', 'deltaTb': 21.16},
        {'name': '=CH2', 'smarts': '[CX3H2][!#1]', 'deltaTb': 18.89},
        {'name': '=CH-', 'smarts': '[CX3H1][!#1]', 'deltaTb': 14.32},
        {'name': '=C', 'smarts': '[CX3H0][!#1]', 'deltaTb': 14.70},
        {'name': '≡CH', 'smarts': '[CX2H1]', 'deltaTb': 10.80},
        {'name': '≡C', 'smarts': '[CX2H0]', 'deltaTb': 10.60},
        {'name': 'Ring CH2', 'smarts': '[R][CH2][R]', 'deltaTb': 21.34},
        {'name': 'Ring CH', 'smarts': '[R][CH][R]', 'deltaTb': 20.00},
        {'name': 'Ring C', 'smarts': '[R][C][R]', 'deltaTb': 17.15},
        {'name': 'OH (alcohol)', 'smarts': '[OX2H][CX4]', 'deltaTb': 34.40},
        {'name': 'O (ether)', 'smarts': '[OX2][CX4]', 'deltaTb': 16.23},
        {'name': 'C=O (carbonyl)', 'smarts': '[CX3]=[OX1]', 'deltaTb': 19.00},
        {'name': 'COOH (carboxylic acid)', 'smarts': 'C(=O)[OH]', 'deltaTb': 58.50},
        {'name': 'NH2 (primary amine)', 'smarts': '[NX3H2]', 'deltaTb': 33.60},
        {'name': 'NH (secondary amine)', 'smarts': '[NX3H1][#6]', 'deltaTb': 22.88},
        {'name': 'N (tertiary amine)', 'smarts': '[NX3]([#6])([#6])[#6]', 'deltaTb': 14.03},
        {'name': 'Cl', 'smarts': '[Cl]', 'deltaTb': 8.76},
        {'name': 'Br', 'smarts': '[Br]', 'deltaTb': 8.51},
        {'name': 'F', 'smarts': '[F]', 'deltaTb': 6.44},
        {'name': 'Phenyl', 'smarts': 'c1ccccc1', 'deltaTb': 19.87},
        {'name': 'NO2 (nitro)', 'smarts': 'N(=O)=O', 'deltaTb': 36.63},
        {'name': 'Sulfoxide', 'smarts': '[#16X3]=[OX1]', 'deltaTb': 30.48},
        {'name': 'Sulfone', 'smarts': '[#16X4](=O)(=O)', 'deltaTb': 52.26},
        {'name': 'Ester', 'smarts': 'C(=O)O[CX4]', 'deltaTb': 31.51},
        {'name': 'Ether', 'smarts': '[OD2]([#6])[#6]', 'deltaTb': 16.23},
        {'name': 'Aldehyde', 'smarts': '[CX3H1]=[O]', 'deltaTb': 19.00},
        {'name': 'Ketone', 'smarts': '[CX3](=O)[#6]', 'deltaTb': 19.00},
        {'name': 'Thiol', 'smarts': '[SX2H]', 'deltaTb': 34.10},
        {'name': 'Thioether', 'smarts': '[#16X2][#6]', 'deltaTb': 21.00},
        # Add more groups if necessary
    ]

    # Initialize total group contribution
    total_deltaTb = 0
    # Keep track of atoms assigned to groups to avoid double counting
    assigned_atoms = set()

    for group in joback_groups:
        pattern = group['smarts']
        smarts = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(smarts)
        count = 0
        for match in matches:
            # Check if any atom in the match has been assigned already
            if any(atom_idx in assigned_atoms for atom_idx in match):
                continue  # Skip if atoms have been assigned
            assigned_atoms.update(match)
            count += 1
        # Sum the contributions
        deltaTb = count * group['deltaTb']
        total_deltaTb += deltaTb
        # Uncomment the following line to see debugging info
        # print(f"Group: {group['name']}, Count: {count}, ΔTb: {deltaTb}")

    # Estimate boiling point
    Tb_K = 198.0 + total_deltaTb  # Boiling point in Kelvin
    Tb_C = Tb_K - 273.15  # Convert to Celsius

    # Compare with 250°C
    if Tb_C <= 250.0:
        return True, f"Estimated boiling point is {Tb_C:.2f}°C, which is ≤ 250°C"
    else:
        return False, f"Estimated boiling point is {Tb_C:.2f}°C, which is > 250°C"