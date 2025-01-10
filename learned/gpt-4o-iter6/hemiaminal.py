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

    # Define enhanced hemiaminal pattern to match more configurations
    # This pattern captures carbons attached to a nitrogen and oxygen, regardless of hybridization
    hemiaminal_patterns = [
        Chem.MolFromSmarts("[C;X4,X3,X2](O)(N)"),
        # Consider patterns where carbon is adjacent to nitrogen and has oxygen,
        # allowing for structural diversity in hemiaminals
        Chem.MolFromSmarts("[C;X4,X3,X2](N)([OH])"),
    ]

    # Check for matches with any defined pattern
    for pattern in hemiaminal_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains hemiaminal structure: amino and hydroxy groups attached to the same carbon"
    
    return False, "No hemiaminal structure found: missing amino and hydroxy groups on the same carbon"