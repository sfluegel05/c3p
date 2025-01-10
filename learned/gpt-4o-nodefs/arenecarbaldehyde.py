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
    # The aldehyde carbon [C](=O) should be connected to an aromatic carbon [c]
    aldehyde_attached_aromatic_pattern = Chem.MolFromSmarts("[c][C](=O)")
    
    # Verify the molecule contains the necessary substructure
    if mol.HasSubstructMatch(aldehyde_attached_aromatic_pattern):
        return True, "Contains aromatic ring with an aldehyde group directly attached"
    else:
        return False, "No arenecarbaldehyde structure found"