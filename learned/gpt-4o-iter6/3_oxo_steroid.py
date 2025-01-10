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

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # A general pattern for the steroid backbone: four fused rings (flexible with aromatic possibilities)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CC3C2CCC4C3CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Pattern to identify an oxo group (C=O) at the 3 position in the context of a steroid
    oxo_3_pattern = Chem.MolFromSmarts("C1(=O)[C;R2]2[C;R1]3C[C;R1]4[C;R1][C;R1][C;R1][C;R1]12")
    if not mol.HasSubstructMatch(oxo_3_pattern):
        return False, "No 3-oxo group found on steroid skeleton"

    return True, "3-oxo group found at position 3 on steroid skeleton"