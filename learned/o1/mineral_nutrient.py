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

    A mineral nutrient is a mineral that is an inorganic nutrient which must be ingested and absorbed in adequate amounts to satisfy a wide range of essential metabolic and/or structural functions in the human body.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # List of essential mineral elements
    essential_minerals = {
        "Ca",  # Calcium
        "P",   # Phosphorus
        "K",   # Potassium
        "Na",  # Sodium
        "Cl",  # Chlorine
        "Mg",  # Magnesium
        "S",   # Sulfur
        "Fe",  # Iron
        "Mn",  # Manganese
        "Cu",  # Copper
        "I",   # Iodine
        "Zn",  # Zinc
        "Se",  # Selenium
        "Mo",  # Molybdenum
        "Cr",  # Chromium
        "Co",  # Cobalt
        "F"    # Fluorine
    }

    # Get set of elements present in the molecule
    elements_in_mol = set()
    for atom in mol.GetAtoms():
        elements_in_mol.add(atom.GetSymbol())

    # Check if molecule contains at least one essential mineral element
    if elements_in_mol.intersection(essential_minerals):
        return True, "Contains essential mineral elements"

    return False, "Does not contain any essential mineral elements"