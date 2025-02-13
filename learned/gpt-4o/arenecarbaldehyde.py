"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde in which the carbonyl group is attached to an aromatic moiety.
    
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

    # Define aldehyde pattern
    aldehyde_pattern = Chem.MolFromSmarts("CO")
    
    # Check if the molecule contains an aldehyde group
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    # Define aromatic carbon pattern with attached aldehyde
    aromatic_aldehyde_pattern = Chem.MolFromSmarts("c[CH]=O")
    
    # Check for the specific structural pattern
    if mol.HasSubstructMatch(aromatic_aldehyde_pattern):
        return True, "Aromatic ring with attached aldehyde group found"
    else:
        return False, "No aromatic ring directly bound to an aldehyde group found"