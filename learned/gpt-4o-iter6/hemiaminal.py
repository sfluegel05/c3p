"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is characterized by an amino group and a hydroxy group attached to the same carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an enhanced SMARTS pattern to capture more potential hemiaminal configurations
    # Focus on flexible hybridization and multiple equivalent representations
    hemiaminal_patterns = [
        Chem.MolFromSmarts("[CX4](O)(N)"),  # Typical carbon with single bonds
        Chem.MolFromSmarts("[CX3](O)(N)"),  # Allow sp2 hybridized carbon
        Chem.MolFromSmarts("[C](O)(N)(*)"), # Handle attachments to additional groups
    ]

    # Check for matches with any defined pattern
    for pattern in hemiaminal_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains hemiaminal structure: amino and hydroxy groups attached to the same carbon"
    
    return False, "No hemiaminal structure found: missing amino and hydroxy groups on the same carbon"