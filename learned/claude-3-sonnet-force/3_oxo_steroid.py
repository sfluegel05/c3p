"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:35483 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is a steroid compound with an oxo substituent at the 3rd position.

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

    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]12[C@@]3([C@H]([C@@H]4[C@@]([C@]([C@]5([C@@]6([C@@H]([C@@H]([C@H]([C@@H]([C@H]([C@H](3)CC[C@@H](2)O)C5)C)C)C(C[C@]46C)(=O)C)C)(C)C)(C)C)C)C(=O)CC1"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for oxo group at position 3
    oxo_pattern = Chem.MolFromSmarts("[C@]12[C@@]3([C@H]([C@@H]4[C@@]([C@]([C@]5([C@@]6([C@@H]([C@@H]([C@H]([C@@H]([C@H]([C@H](3)CC[C@@H](2)O)C5)C)C)C(C[C@]46C)(=O)C)C)(C)C)(C)C)C)C(=O)CC1=O"
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group found at position 3"

    return True, "Contains a steroid backbone with an oxo substituent at position 3"