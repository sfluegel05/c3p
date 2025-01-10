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

    # Look for a primary alcohol group attached directly to an aromatic ring
    primary_alcohol_aromatic_pattern = Chem.MolFromSmarts("a[CX4H2][OX2H]")
    
    if not mol.HasSubstructMatch(primary_alcohol_aromatic_pattern):
        return False, "No primary alcohol group directly attached to an aromatic ring found"

    return True, "Contains a primary alcohol group directly attached to an aromatic ring"