"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.

    A hemiaminal contains both a hydroxyl group (-OH) and an amine group (-NH2 or -NH)
    attached to the same carbon atom, including cyclic variations.

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

    # Improved pattern for hemiaminal:
    # Allow for both primary and secondary amines, and ring structures
    hemiaminal_patterns = [
        "[C;R0]([OH])([NH2])",   # Standard hemiaminal structure
        "[C;R0]([OH])([NH])",    # Secondary amine variant
        "[C&!R1]([OH])N",        # Exclude direct C-N bonds that form rings, using NSMARTS
        "[$(C([OH])[NH,NH2])]"] # Check for broader environment, may find both 

    # Iterate through hemiaminal_patterns to find a match
    for pattern in hemiaminal_patterns:
        hemiaminal_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(hemiaminal_pattern):
            return True, f"Contains hemiaminal substructure pattern matching {pattern}"

    return False, "Does not contain detectable hemiaminal substructure"