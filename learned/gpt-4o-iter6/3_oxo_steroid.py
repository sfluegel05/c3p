"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Look for the steroid backbone structure (4-ring system typical of steroids)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C(=O)CCC4=C3C2CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Find and verify the presence of an oxo group (=O) at position 3
    oxo_pattern = Chem.MolFromSmarts("C1=CC(=O)[C@@H]2CC[C@H]3C=C4C(=O)CCC4CCC3=C12")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found on steroid skeleton"

    return True, "3-oxo group found at position 3 on steroid skeleton"