"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxy group at the beta- or 3-position
    from the carboxylic acid and is characterized by a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive 3-hydroxy fatty acid pattern
    # This pattern looks for a carboxylic acid followed by two carbons and an OH group
    pattern = Chem.MolFromSmarts("C(=O)O[C](C)[C@H](O)")

    # Check if the molecule matches the 3-hydroxy fatty acid pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a hydroxy group at the 3-position"
    
    # Also consider non-chiral variants
    pattern_non_chiral = Chem.MolFromSmarts("C(=O)OCC(O)")
    if mol.HasSubstructMatch(pattern_non_chiral):
        return True, "Contains a hydroxy group at the 3-position (non-chiral)"

    return False, "No hydroxy group at the 3-position found"