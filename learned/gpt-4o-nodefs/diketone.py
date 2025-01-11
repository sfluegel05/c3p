"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as having exactly two distinct ketone groups (C=O), 
    where each ketone carbon is only bonded to carbon (C) atoms on the sides,
    thus ensuring the functional structure of typical diketones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a ketone pattern with more strict neighbor conditions
    # This pattern: carbon double-bonded to oxygen, bonded to two distinct carbon atoms
    ketone_pattern = Chem.MolFromSmarts("C(=O)(C)C")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Count the number of distinct ketone groups
    ketone_count = len(ketone_matches)
    
    if ketone_count == 2:
        return True, "Contains exactly two ketone groups"
    else:
        return False, f"Found {ketone_count} ketone groups, need exactly 2"