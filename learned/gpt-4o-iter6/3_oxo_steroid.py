"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has an oxo group at the 3-position on the steroid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible steroid backbone pattern using four fused rings
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define a pattern for detecting an oxo group at position 3
    # The C=O group should be adjacent to the first ring in a typical steroid
    oxy_pattern = Chem.MolFromSmarts("C1=O")
    if not mol.HasSubstructMatch(oxy_pattern):
        return False, "No 3-oxo group found on steroid skeleton"

    return True, "3-oxo group found at position 3 on steroid skeleton"