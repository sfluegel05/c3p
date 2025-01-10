"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.

    A hemiaminal contains both a hydroxyl group (-OH) and an amine group (-NH2 or -NH)
    attached to the same carbon atom, including consideration of stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined pattern for hemiaminal:
    # Ensure patterns accurately reflect plausible hemiaminal configurations
    hemiaminal_patterns = [
        "[C;!R]([OH])[NH2]",     # Non-cyclic, open-chain hemiaminal
        "[C;!R]([OH])[NH]",      # Non-cyclic secondary amine variant
        "[C&R]([OH])[NH2]",      # Cyclic, primary amine variant
        "[C&R]([OH])[NH]",       # Cyclic, secondary amine variant
    ]

    # Iterate through the refined hemiaminal patterns to find a match
    for pattern in hemiaminal_patterns:
        hemiaminal_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(hemiaminal_pattern):
            return True, f"Contains hemiaminal substructure pattern matching {pattern}"

    return False, "Does not contain detectable hemiaminal substructure"