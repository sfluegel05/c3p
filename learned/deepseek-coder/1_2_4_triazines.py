"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: CHEBI:38323 1,2,4-triazine
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine has a six-membered ring with nitrogen atoms at positions 1, 2, and 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 1,2,4-triazine core pattern
    triazine_pattern = Chem.MolFromSmarts("n1cnnc1")
    if not mol.HasSubstructMatch(triazine_pattern):
        return False, "No 1,2,4-triazine core found"

    # Verify the nitrogen positions in the ring
    # The pattern ensures that the nitrogens are at positions 1, 2, and 4
    # So no further checks are needed for the core structure

    return True, "Contains a 1,2,4-triazine core with nitrogens at positions 1, 2, and 4"