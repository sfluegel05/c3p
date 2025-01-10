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
    A mineral nutrient is an inorganic compound containing essential minerals required by the human body.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    # Essential mineral elements required by the human body
    essential_minerals = {
        11,  # Na - Sodium
        12,  # Mg - Magnesium
        13,  # Al - Aluminium
        19,  # K - Potassium
        20,  # Ca - Calcium
        25,  # Mn - Manganese
        26,  # Fe - Iron
        29,  # Cu - Copper
        30,  # Zn - Zinc
        34,  # Se - Selenium
        38,  # Sr - Strontium
        56,  # Ba - Barium
        74,  # W - Tungsten
        79,  # Au - Gold
        80,  # Hg - Mercury
        82,  # Pb - Lead
        83,  # Bi - Bismuth
    }

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is inorganic (contains no carbon atoms)
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if c_atoms:
        return False, "Molecule is organic"

    # Get set of atomic numbers present in the molecule
    atomic_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}

    # Check for presence of essential minerals
    minerals_in_mol = atomic_nums & essential_minerals
    if not minerals_in_mol:
        return False, "Does not contain essential mineral elements required by the human body"

    # Passed all checks
    return True, "Molecule is an inorganic compound containing essential mineral nutrients"