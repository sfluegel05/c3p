"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is identified as a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for porphyrin-like structures
    # General porphyrin has four pyrrole-like units connected by methine bridges
    porphyrin_pattern = Chem.MolFromSmarts("[n]1c2cc[nH]c2ccc3c4cc[nH]c4ccc5c6cc[nH]c6ccc1c35")

    # Check if the molecule matches the porphyrin core pattern
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "Matches porphyrin macrocyclic structure"

    return False, "Does not match porphyrin structure"