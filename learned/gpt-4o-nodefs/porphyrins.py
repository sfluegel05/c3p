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
    
    # Define a SMARTS pattern for core porphyrin structure
    # A broader pattern is used to account for variations across real examples
    # This pattern is intended to match typical porphyrin macrocycles with or without metals
    porphyrin_pattern = Chem.MolFromSmarts('[n&H]1ccc(c2nc(c3ccc([n&H]c4ccc([n&H]c1c34)c2)=C)=C)')

    # Check if the molecule matches the porphyrin core pattern
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "Matches porphyrin macrocyclic structure"

    return False, "Does not match porphyrin structure"