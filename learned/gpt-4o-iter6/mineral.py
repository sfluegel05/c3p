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
    
    # Define relevant metals and non-metals frequently found in minerals
    elements = {'Na', 'K', 'Ca', 'Mg', 'Fe', 'Cu', 'Ba', 'Cs', 'La', 'Al', 'Si', 'Sb', 'Ni', 'Zn', 'P', 'Cl', 'F', 'O', 'S'}
    element_count = sum(atom.GetSymbol() in elements for atom in mol.GetAtoms())

    # If at least one major mineral element is found, proceed with checks
    if element_count < 1:
        return False, "No typical mineral components detected"
    
    # Define anion patterns common in minerals
    anion_patterns = [
        Chem.MolFromSmarts('[O-]'),  # Oxide pattern
        Chem.MolFromSmarts('P(=O)([O-])[O-]'),  # Phosphate pattern
        Chem.MolFromSmarts('S(=O)([O-])[O-]'),  # Sulfate pattern
        Chem.MolFromSmarts('[C](=O)([O-])'),  # Carbonate pattern
        Chem.MolFromSmarts('[F-]'),  # Fluoride pattern
        Chem.MolFromSmarts('[Cl-]')  # Chloride pattern
    ]
    
    for pattern in anion_patterns:
        if mol.HasSubstructMatch(pattern) and element_count >= 1:
            return True, "Contains characteristic anions and elements typical of minerals"
    
    # Water molecule pattern for hydrates
    water_pattern = Chem.MolFromSmarts('O')
    water_matches = len(mol.GetSubstructMatches(water_pattern))
    if water_matches >= 3:
        return True, "Contains multiple water molecules characteristic of hydrated minerals"
    
    return False, "Does not meet typical mineral structure criteria"