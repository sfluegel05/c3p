"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:35755 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid has a single carboxylic acid group with a ketone group
    at the alpha position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if len(carboxylic_matches) == 0:
        return False, "No carboxylic acid group found"
    elif len(carboxylic_matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Check for 2-oxo pattern (ketone adjacent to carboxylic acid)
    # Pattern: O=C-C(=O)OH
    oxo_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(oxo_acid_pattern):
        return False, "No 2-oxo group adjacent to carboxylic acid"
    
    # Additional check to ensure the ketone is not part of another functional group
    # like an ester or another carboxylic acid
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "The 2-oxo group is part of another functional group"

    return True, "Contains a carboxylic acid group with adjacent 2-oxo substituent"