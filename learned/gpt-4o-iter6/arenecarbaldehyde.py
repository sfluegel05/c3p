"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as any aldehyde with the carbonyl group attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an aldehyde pattern where the carbonyl carbon is attached to an aromatic carbon
    aldehyde_aromatic_attachment = Chem.MolFromSmarts("[C](=[O])[c]")

    if mol.HasSubstructMatch(aldehyde_aromatic_attachment):
        return True, "Aldehyde group attached to an aromatic moiety"
    
    return False, "Aldehyde group not attached to aromatic moiety"