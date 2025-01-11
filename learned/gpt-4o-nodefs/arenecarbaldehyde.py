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
    
    # SMARTS pattern for an aldehyde group directly attached to an aromatic carbon
    aldehyde_attached_aromatic_pattern = Chem.MolFromSmarts("[#6]=O")
    
    # Check if there's an aromatic carbon bounded to an aldehyde group
    aromatic_and_aldehyde_pattern = Chem.MolFromSmarts("[$([c]C=O)]")
    
    # Check if molecule has aromatic ring matched with aldehyde group directly
    if mol.HasSubstructMatch(aromatic_and_aldehyde_pattern):
        return True, "Contains aromatic ring with an aldehyde group directly attached"
    else:
        return False, "No arenecarbaldehyde structure found"