"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a mineral, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of metal cations (basic metal types: Alkali, Alkaline earth, transition metals)
    metals = ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Zn', 'Cu', 'Ba', 'Cs', 'La', 'Sb', 'Ni', 'Al']
    metal_count = 0

    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            metal_count += 1

    if metal_count < 1:
        return False, "No metal ions detected, unlikely to be a mineral"

    # Check if there are counter ions or complex anions often found in minerals
    anion_patterns = [
        Chem.MolFromSmarts('[O-]', asSmarts=True),  # Common oxide or simple anionic pattern
        Chem.MolFromSmarts('P([O-])(=O)', asSmarts=True),  # Phosphate type pattern
        Chem.MolFromSmarts('S([O-])(=O)', asSmarts=True)  # Sulfate type pattern
    ]
    anion_count = 0

    for pattern in anion_patterns:
        if mol.HasSubstructMatch(pattern):
            anion_count += 1

    if anion_count > 0:
        return True, "Contains characteristic anions and metal ions typical of minerals"

    # Check for hydrate structures (e.g., multiple water molecules present)
    water_pattern = Chem.MolFromSmarts("O")
    water_matches = mol.GetSubstructMatches(water_pattern)
    if len(water_matches) >= 4:
        return True, "Contains multiple water molecules similar to hydrates found in minerals"

    return False, "Does not meet typical mineral structure criteria"