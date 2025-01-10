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

    # Check for presence of metal cations or major mineral components
    elements = ['Na', 'K', 'Ca', 'Mg', 'Fe', 'Zn', 'Cu', 'Ba', 'Cs', 'La', 'Sb', 'Ni', 'Al', 'Si']
    element_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in elements:
            element_count += 1

    if element_count < 1:
        return False, "No typical mineral components detected"

    # Check for mineral cationic or anionic patterns
    anion_patterns = [
        Chem.MolFromSmarts('[O-]'),  # Oxide pattern
        Chem.MolFromSmarts('P(=O)([O-])[O-]'),  # Phosphate pattern
        Chem.MolFromSmarts('S(=O)([O-])[O-]'),  # Sulfate pattern
        Chem.MolFromSmarts('[C](=O)([O-])'),  # Carbonate pattern
    ]

    for pattern in anion_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains characteristic anions and elements typical of minerals"

    # Check for hydrate or water patterns suggesting geological or hydrated forms
    water_pattern = Chem.MolFromSmarts('O')
    water_matches = mol.GetSubstructMatches(water_pattern)
    if len(water_matches) >= 3:
        return True, "Contains multiple water molecules characteristic of hydrated minerals"

    return False, "Does not meet typical mineral structure criteria"