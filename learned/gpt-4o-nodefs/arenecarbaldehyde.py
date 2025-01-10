"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is characterized by an aldehyde group directly attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for aromatic ring with aldehyde group directly attached
    arenecarbaldehyde_pattern = Chem.MolFromSmarts("[$([C;R]=O)]")
    
    # Check if the pattern is present in the molecule
    if mol.HasSubstructMatch(arenecarbaldehyde_pattern):
        return True, "Contains aromatic ring with an aldehyde group directly attached"
    else:
        return False, "No arenecarbaldehyde structure found"