"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    This function checks for the presence of metal ions and specific anionic groups
    typically found in mineral nutrients.

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

    # List of metal ions typically found in mineral nutrients
    metal_ions = {'Li', 'Na', 'K', 'Mg', 'Ca', 'Fe', 'Cu', 'Zn', 'Ba', 'La', 'Cs', 'Al', 'Mn', 'Co', 'Ni'}

    # Correct SMARTS patterns for common anions
    common_anions_patterns = [
        "[OX2-]",                 # Simplified example for an oxide ion
        "OS(=O)(=O)[O-]",         # Sulfate ion
        "O[P](=O)([O-])[O-]",     # Phosphate ion
        "[N+](=O)[O-]",           # Nitrate ion
        "[C](=O)([O-])[O-]",      # Carbonate ion
        "[Cl-]", "[F-]", "[Br-]", "[I-]" # Halides
    ]

    # Check for presence of metal ions
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if not any(metal in atoms for metal in metal_ions):
        return False, "No essential metal ion typically found in mineral nutrients"

    # Check for presence of specific anionic groups
    any_match = False
    for anion_smarts in common_anions_patterns:
        anion_mol = Chem.MolFromSmarts(anion_smarts)
        if anion_mol is not None and mol.HasSubstructMatch(anion_mol):
            any_match = True
            break

    if any_match:
        return True, "Contains metal ion and specific anion groups typical of mineral nutrients"

    return False, "No specific anionic groups found typical of mineral nutrients"