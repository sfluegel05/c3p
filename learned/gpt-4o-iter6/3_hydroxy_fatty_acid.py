"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxy functional group at the beta- or 3-position
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

    # SMARTS: C(=O)O (the carboxylic acid part) - C - C(O)
    pattern = Chem.MolFromSmarts("C(=O)O-C-C(O)")

    if mol.HasSubstructMatch(pattern):
        return True, "Contains a hydroxy group at the 3-position"
    
    return False, "No hydroxy group at the 3-position found"