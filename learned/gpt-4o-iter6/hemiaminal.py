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
    
    # Expand SMARTS pattern for hemiaminal to include both cyclic and acyclic cases
    # Match carbon with hydroxy and amino group while being flexible for cyclic systems and additional groups
    hemiaminal_pattern_linear = Chem.MolFromSmarts("[CX4](O)(N)")  # For linear hemiaminals
    hemiaminal_pattern_cyclic = Chem.MolFromSmarts("[C](O)(N)")  # Accounts for cyclic cases

    # Check for matches with either pattern
    if mol.HasSubstructMatch(hemiaminal_pattern_linear) or mol.HasSubstructMatch(hemiaminal_pattern_cyclic):
        return True, "Contains hemiaminal structure: amino and hydroxy groups attached to the same carbon"
    
    return False, "No hemiaminal structure found: missing amino and hydroxy groups on the same carbon"