"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is a mineral that is an inorganic nutrient essential for the human body.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """

    # Essential mineral elements required by the human body
    essential_minerals = {
        9,   # F - Fluorine
        11,  # Na - Sodium
        12,  # Mg - Magnesium
        13,  # Al - Aluminum
        15,  # P - Phosphorus
        16,  # S - Sulfur
        17,  # Cl - Chlorine
        19,  # K - Potassium
        20,  # Ca - Calcium
        26,  # Fe - Iron
        29,  # Cu - Copper
        30,  # Zn - Zinc
        51,  # Sb - Antimony
        53,  # I - Iodine
        55,  # Cs - Cesium
        56,  # Ba - Barium
        57,  # La - Lanthanum
    }

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get set of atomic numbers present in the molecule
    atomic_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}

    # Check for presence of essential mineral elements
    minerals_in_mol = atomic_nums & essential_minerals
    if not minerals_in_mol:
        return False, "Does not contain essential mineral elements required by the human body"

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # If number of carbon atoms is high, likely an organic molecule
    if num_carbons > 6:
        return False, f"Contains too many carbon atoms ({num_carbons}), likely an organic molecule"

    # Check for aromatic rings, characteristic of organic molecules
    if mol.GetAromaticAtomCount() > 0:
        return False, "Contains aromatic rings, likely an organic molecule"

    # Check for presence of functional groups characteristic of organic molecules
    # Define SMARTS patterns for common organic functional groups
    organic_functional_groups = [
        '[#6][#7]',  # Carbon attached to nitrogen
        '[#6][#8]',  # Carbon attached to oxygen
        '[#6]=[#6]', # Carbon-carbon double bond
        '[#6]#[#6]', # Carbon-carbon triple bond
        '[#6]=[O]',  # Carbonyl group
        '[#6]-[#6]-[#6]',  # Carbon chains of length 3 or more
    ]

    for pattern in organic_functional_groups:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            return False, "Contains functional groups characteristic of organic molecules"

    # Passed all checks
    return True, "Molecule is an inorganic compound containing essential mineral nutrients"