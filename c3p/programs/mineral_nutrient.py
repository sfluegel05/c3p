"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is an inorganic compound essential for the human body, typically containing a metal ion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define common metal ion pattern
    metal_ions_pattern = Chem.MolFromSmarts("[Al, Ba, Ca, Cs, Fe, K, La, Mg, Na, Pd, Sb, Zn]")

    # Define common counterion pattern including sulfates, chlorides, etc.
    counterions_pattern = Chem.MolFromSmarts("[Cl, F, O, P, S, N]")

    # Check for presence of metal ions
    if not mol.HasSubstructMatch(metal_ions_pattern):
        return False, "No metal ions found that are commonly essential nutrients"

    # Check for presence of common counterions
    if not mol.HasSubstructMatch(counterions_pattern):
        return False, "No recognizable counterions found"

    # Additional logic can include specific checks for known essential nutrient patterns.
    # For simplicity, here we check that at least one metal ion and one counterion is present.
    return True, "Contains metal ions with suitable anionic components typical of mineral nutrients"