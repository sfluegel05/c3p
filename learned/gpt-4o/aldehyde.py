"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a compound contains an aldehyde group (RC(=O)H) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains an aldehyde group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for an aldehyde group: carbonyl group (C=O) directly bonded to hydrogen (C-H)
    # This is a more general pattern to allow for variability in surrounding atoms, but focuses on C=O with C-H
    aldehyde_pattern = Chem.MolFromSmarts("[CX3]=O")
    hydrogen_pattern = Chem.MolFromSmarts("[CX2H]")  # Must be bonded to hydrogen as well

    # Check if the molecule contains the aldehyde functional group
    if mol.HasSubstructMatch(aldehyde_pattern) and mol.HasSubstructMatch(hydrogen_pattern):
        return True, "Contains an aldehyde group (RC(=O)H)"
    
    return False, "No aldehyde group found"