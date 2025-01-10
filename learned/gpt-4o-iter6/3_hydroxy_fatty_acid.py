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
    # C(=O)O represents the carboxylic acid.
    # [C;H2][C;H](O) represents a methine group with a hydroxy group at 3-position.
    # The long chain is implicit and these patterns should suffice to wrap around it.
    pattern = Chem.MolFromSmarts("C(=O)O[C;H2][C;H](O)")

    if mol.HasSubstructMatch(pattern):
        return True, "Contains a hydroxy group at the 3-position"
    
    return False, "No hydroxy group at the 3-position found"