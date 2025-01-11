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
    A mineral nutrient is a compound containing essential minerals required by the human body.

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

    # Toxic or non-essential elements to exclude
    toxic_elements = {
        33, # As - Arsenic
        74, # W - Tungsten
        80, # Hg - Mercury
        82, # Pb - Lead
    }

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get set of atomic numbers present in the molecule
    atomic_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}

    # Check for presence of toxic elements
    if atomic_nums & toxic_elements:
        return False, "Contains toxic or non-essential elements"

    # Check for presence of essential minerals
    minerals_in_mol = atomic_nums & essential_minerals
    if not minerals_in_mol:
        return False, "Does not contain essential mineral elements required by the human body"

    # Passed all checks
    return True, "Molecule contains essential mineral nutrients"