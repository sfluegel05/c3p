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

    # Enhanced SMARTS patterns for hemiaminal detection
    hemiaminal_patterns = [
        Chem.MolFromSmarts("[CX4;R0](O)(N)"),  # Accommodating both open-chain and ring hemiaminals
        Chem.MolFromSmarts("[CX3;R0](O)(N)"),  # Permit sp2 carbon within non-ring systems
        Chem.MolFromSmarts("[R][CX4](O)(N)"),  # Cyclic with an aliphatic carbon
        Chem.MolFromSmarts("[R][CX3](O)(N)"),  # Cyclic with a planar carbon
    ]

    # Check for matches with any defined pattern
    for pattern in hemiaminal_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains hemiaminal structure: amino and hydroxy groups attached to the same carbon"
    
    return False, "No hemiaminal structure found: missing amino and hydroxy groups on the same carbon"