"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:35713 3-oxo steroid
Definition: Any oxo steroid where an oxo substituent is located at position 3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.

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
    
    # Look for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2[C@@]3([C@H](C[C@@H]4CC=C5C[C@@H](CC[C@]5(C)[C@H]4C3)C2)C1)C"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for oxo group at position 3
    oxo_pattern = Chem.MolFromSmarts("[C](=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H](CC3)[C@@H]12")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group at position 3"
    
    # Additional checks
    mol_wt = Chem.Descriptors.MolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, "Molecular weight outside typical range for 3-oxo steroids"
    
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Fewer than 3 rings, not a steroid"
    
    return True, "Contains steroid backbone with oxo group at position 3"