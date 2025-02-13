"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("a")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"
    
    # Look for a primary alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[CX4;H2][OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No primary alcohol (C-OH) group found"

    return True, "Contains an aromatic ring and a primary alcohol group"

# Testing example
# is_aromatic_primary_alcohol("C1=CC(=CN=C1)CO")  # 3-pyridinemethanol, expected True.