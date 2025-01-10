"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is characterized as a cyclic hemiacetal, featuring a carbon atom
    bonded to both a hydroxyl group and an ether oxygen within a ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Refined lactol pattern
    # Ensures a carbon is bonded to both an hydroxyl group and an ether oxygen, and both oxygens are part of the ring
    lactol_smart_pattern = '[C;R]([O;R;H1])[O;R]'

    lactol_pattern = Chem.MolFromSmarts(lactol_smart_pattern)
    
    if lactol_pattern and mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains cyclic hemiacetal (lactol) structure"
    else:
        return False, "Does not contain cyclic hemiacetal (lactol) structure"

# Example usage:
# result, reason = is_lactol("O1C2=C(C(C[C@@]1(C3=CC=...=C2)O)O")
# print(result, reason)